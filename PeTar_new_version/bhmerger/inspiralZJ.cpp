#include <Python.h>
#include <numpy/arrayobject.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#ifndef PETAR_PRECESSION_SOURCE_DIR
#define PETAR_PRECESSION_SOURCE_DIR ""
#endif
#ifndef PETAR_PRECESSION_INSTALL_DIR
#define PETAR_PRECESSION_INSTALL_DIR ""
#endif

#include "inspiralZJ.hpp"

namespace {

// Only cancellation at the circular limit may be rounded to zero.
bool eccentricityFromEnergyTerm(double term, double& eccentricity) {
    if (!std::isfinite(term)) return false;
    const double e2=1.0+term;
    const double tolerance=32.0*std::numeric_limits<double>::epsilon()
                         *std::max(1.0,std::abs(term));
    if (e2 < -tolerance) return false;
    eccentricity=std::sqrt(std::max(0.0,e2));
    return true;
}

bool configurePrecessionPath() {
    const char* override_path=std::getenv("PETAR_PRECESSION_PATH");
    const char* candidates[]={override_path, PETAR_PRECESSION_SOURCE_DIR,
                              PETAR_PRECESSION_INSTALL_DIR};
    for (const char* candidate: candidates) {
        if (!candidate || !candidate[0]) continue;
        std::error_code error;
        const auto path=std::filesystem::absolute(candidate,error);
        if (!error && std::ifstream(path/"precession"/"__init__.py").good()) {
            PyObject* python_path=PySys_GetObject("path"); // borrowed
            PyObject* entry=PyUnicode_DecodeFSDefault(path.c_str());
            if (!entry) return false;
            const bool ok=python_path && PyList_Insert(python_path,0,entry)==0;
            Py_DECREF(entry);
            return ok;
        }
        if (candidate==override_path) {
            std::fprintf(stderr,"Invalid PETAR_PRECESSION_PATH: %s\n",candidate);
            return false; // An explicit override must not silently select another package.
        }
    }
    std::fprintf(stderr,"Cannot locate bundled precession; set PETAR_PRECESSION_PATH\n");
    return false;
}

bool setDoubleKeyword(PyObject* dict, const char* key, double value) {
    if (!dict || PyErr_Occurred()) return false;
    PyObject* number = PyFloat_FromDouble(value);
    if (!number) return false;
    const int result = PyDict_SetItemString(dict, key, number);
    Py_DECREF(number);
    return result == 0;
}

PyObject* callKeywords(PyObject* function, PyObject* kwargs) {
    if (!kwargs || PyErr_Occurred()) return nullptr;
    PyObject* args = PyTuple_New(0);
    if (!args) return nullptr;
    PyObject* result = PyObject_Call(function, args, kwargs);
    Py_DECREF(args);
    return result;
}

bool validArray(PyObject* object, npy_intp size) {
    if (!object || !PyArray_Check(object)) return false;
    auto* array = reinterpret_cast<PyArrayObject*>(object);
    if (PyArray_TYPE(array) != NPY_DOUBLE || !PyArray_ISCARRAY_RO(array)
        || !PyArray_ISNOTSWAPPED(array) || PyArray_SIZE(array) != size) return false;
    const double* data = static_cast<const double*>(PyArray_DATA(array));
    for (npy_intp i=0; i<size; ++i) if (!std::isfinite(data[i])) return false;
    return true;
}

double arrayValue(PyObject* object, bool two_dimensional, bool circular_e=false) {
    if (!object || !PyArray_Check(object)) return NAN;
    auto* array = reinterpret_cast<PyArrayObject*>(object);
    if (PyArray_NDIM(array) != (two_dimensional ? 2 : 1)
        || PyArray_DIM(array, 0) < 1
        || (two_dimensional && PyArray_DIM(array, 1) < (circular_e ? 1 : 2))) return NAN;
    const bool constant_e = two_dimensional && PyArray_DIM(array, 1)==1;
    char* address = static_cast<char*>(two_dimensional
                                      ? PyArray_GETPTR2(array,0,PyArray_DIM(array,1)-1)
                                                       : PyArray_GETPTR1(array,0));
    PyObject* item = PyArray_GETITEM(array, address);
    if (!item) return NAN;
    const double value = PyFloat_AsDouble(item);
    Py_DECREF(item);
    return constant_e && value!=0.0 ? NAN : value;
}

bool setStringListKeyword(PyObject* dict, const char* key,
                          const char* const values[], Py_ssize_t count) {
    if (!dict || PyErr_Occurred()) return false;
    PyObject* list = PyList_New(count);
    if (!list) return false;
    for (Py_ssize_t i=0; i<count; ++i) {
        PyObject* value = PyUnicode_FromString(values[i]);
        if (!value) {
            Py_DECREF(list);
            return false;
        }
        PyList_SET_ITEM(list, i, value); // steals value
    }
    const int result = PyDict_SetItemString(dict, key, list);
    Py_DECREF(list);
    return result == 0;
}

struct InspiralParameters {
    double theta1;
    double theta2;
    double deltaphi;
    double deltachi;
    double kappa;
    double a;
    double e;
    double uc;
    double u;
    double chieff;
    double q;
    double chi1;
    double chi2;
    double m;
};
thread_local InspiralParameters Inspiraloutputs{};
thread_local BHMerger::RemnantParameters Remnantoutputs{};
bool bhmerger_python_owned = false;
std::atomic<bool> bhmerger_api_ready{false};
PyThreadState* bhmerger_main_thread_state = nullptr;

// Recover the adiabatically evolved constants without requesting a random
// nutation phase.  If the two physical deltachi roots coalesce, reconstruct
// the angles at the fitting boundary from the resonant limiting solution.
bool recoverResonantInspiral(PyObject* module,
                             const InspiralParameters& encounter,
                             InspiralParameters& recovered) {
    if (!module) return false;

    PyObject* inspiral = PyObject_GetAttrString(module, "inspiral_precav");
    if (!inspiral || !PyCallable_Check(inspiral)) {
        Py_XDECREF(inspiral);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }

    npy_intp dims[1] = {2};
    double a_data[2] = {encounter.a, 10.0};
    PyObject* a_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, a_data);
    PyObject* kwargs = PyDict_New();
    setDoubleKeyword(kwargs, "theta1", encounter.theta1);
    setDoubleKeyword(kwargs, "theta2", encounter.theta2);
    setDoubleKeyword(kwargs, "deltaphi", encounter.deltaphi);
    if (kwargs && a_array) PyDict_SetItemString(kwargs, "a", a_array);
    setDoubleKeyword(kwargs, "e", encounter.e);
    setDoubleKeyword(kwargs, "q", encounter.q);
    setDoubleKeyword(kwargs, "chi1", encounter.chi1);
    setDoubleKeyword(kwargs, "chi2", encounter.chi2);
    const char* outputs[] = {
        "kappa", "a", "e", "uc", "u", "deltachiminus",
        "deltachiplus", "deltachi3", "chieff", "q", "chi1", "chi2"
    };
    setStringListKeyword(kwargs, "requested_outputs", outputs,
                         static_cast<Py_ssize_t>(sizeof(outputs)/sizeof(outputs[0])));

    PyObject* result = callKeywords(inspiral, kwargs);
    Py_XDECREF(kwargs);
    Py_XDECREF(a_array);
    Py_DECREF(inspiral);
    if (!result || !PyDict_Check(result)) {
        Py_XDECREF(result);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }

    InspiralParameters candidate = encounter;
    candidate.kappa = arrayValue(PyDict_GetItemString(result, "kappa"), true);
    candidate.a = arrayValue(PyDict_GetItemString(result, "a"), true);
    candidate.e = arrayValue(PyDict_GetItemString(result, "e"), true, encounter.e==0.0);
    candidate.uc = arrayValue(PyDict_GetItemString(result, "uc"), true);
    candidate.u = arrayValue(PyDict_GetItemString(result, "u"), true);
    const double root_minus = arrayValue(
        PyDict_GetItemString(result, "deltachiminus"), true);
    const double root_plus = arrayValue(
        PyDict_GetItemString(result, "deltachiplus"), true);
    const double root_three = arrayValue(
        PyDict_GetItemString(result, "deltachi3"), true);
    candidate.chieff = arrayValue(PyDict_GetItemString(result, "chieff"), false);
    candidate.q = arrayValue(PyDict_GetItemString(result, "q"), false);
    candidate.chi1 = arrayValue(PyDict_GetItemString(result, "chi1"), false);
    candidate.chi2 = arrayValue(PyDict_GetItemString(result, "chi2"), false);
    Py_DECREF(result);
    if (PyErr_Occurred()) PyErr_Clear();

    const bool constants_valid =
        std::isfinite(candidate.kappa) && std::isfinite(candidate.a)
        && std::isfinite(candidate.e) && std::isfinite(candidate.uc)
        && std::isfinite(candidate.u) && std::isfinite(candidate.chieff)
        && std::isfinite(candidate.q) && std::isfinite(candidate.chi1)
        && std::isfinite(candidate.chi2) && std::isfinite(root_minus)
        && std::isfinite(root_plus) && std::isfinite(root_three)
        && candidate.a > 0.0 && candidate.e >= 0.0 && candidate.e < 1.0
        && candidate.q > 0.0 && candidate.q <= 1.0
        && candidate.chi1 >= 0.0 && candidate.chi1 <= 1.0
        && candidate.chi2 >= 0.0 && candidate.chi2 <= 1.0;
    if (!constants_valid) return false;

    // root_three is not one of the two physical turning points and can become
    // very large near equal masses; it must not relax this coalescence test.
    const double root_scale = std::max(
        {1.0, std::fabs(root_minus), std::fabs(root_plus)});
    const double root_gap = std::fabs(root_plus-root_minus);
    const double resonance_tolerance =
        64.0*std::sqrt(std::numeric_limits<double>::epsilon())*root_scale;
    if (root_gap > resonance_tolerance) return false;

    PyObject* resonance = PyObject_GetAttrString(module, "deltachiresonance");
    if (!resonance || !PyCallable_Check(resonance)) {
        Py_XDECREF(resonance);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }
    kwargs = PyDict_New();
    setDoubleKeyword(kwargs, "kappa", candidate.kappa);
    setDoubleKeyword(kwargs, "a", candidate.a);
    setDoubleKeyword(kwargs, "e", candidate.e);
    setDoubleKeyword(kwargs, "chieff", candidate.chieff);
    setDoubleKeyword(kwargs, "q", candidate.q);
    setDoubleKeyword(kwargs, "chi1", candidate.chi1);
    setDoubleKeyword(kwargs, "chi2", candidate.chi2);
    result = callKeywords(resonance, kwargs);
    Py_XDECREF(kwargs);
    Py_DECREF(resonance);
    if (!validArray(result, 1)) {
        Py_XDECREF(result);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }
    candidate.deltachi = *static_cast<double*>(
        PyArray_DATA(reinterpret_cast<PyArrayObject*>(result)));
    Py_DECREF(result);
    const double root_mean = 0.5*(root_minus+root_plus);
    if (!std::isfinite(candidate.deltachi)
        || std::fabs(candidate.deltachi-root_mean) > 16.0*resonance_tolerance)
        return false;

    PyObject* to_angles = PyObject_GetAttrString(module, "conserved_to_angles");
    if (!to_angles || !PyCallable_Check(to_angles)) {
        Py_XDECREF(to_angles);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }
    kwargs = PyDict_New();
    setDoubleKeyword(kwargs, "deltachi", candidate.deltachi);
    setDoubleKeyword(kwargs, "kappa", candidate.kappa);
    setDoubleKeyword(kwargs, "a", candidate.a);
    setDoubleKeyword(kwargs, "e", candidate.e);
    setDoubleKeyword(kwargs, "chieff", candidate.chieff);
    setDoubleKeyword(kwargs, "q", candidate.q);
    setDoubleKeyword(kwargs, "chi1", candidate.chi1);
    setDoubleKeyword(kwargs, "chi2", candidate.chi2);
    setDoubleKeyword(kwargs, "cyclesign", encounter.deltaphi < 0.0 ? -1.0 : 1.0);
    result = callKeywords(to_angles, kwargs);
    Py_XDECREF(kwargs);
    Py_DECREF(to_angles);
    if (!validArray(result, 3)) {
        Py_XDECREF(result);
        if (PyErr_Occurred()) PyErr_Clear();
        return false;
    }
    const double* angles = static_cast<const double*>(
        PyArray_DATA(reinterpret_cast<PyArrayObject*>(result)));
    candidate.theta1 = angles[0];
    candidate.theta2 = angles[1];
    candidate.deltaphi = angles[2];
    Py_DECREF(result);

    const double pi = std::acos(-1.0);
    if (!std::isfinite(candidate.theta1) || !std::isfinite(candidate.theta2)
        || !std::isfinite(candidate.deltaphi)
        || candidate.theta1 < 0.0 || candidate.theta1 > pi
        || candidate.theta2 < 0.0 || candidate.theta2 > pi
        || candidate.deltaphi < -pi || candidate.deltaphi > pi)
        return false;

    recovered = candidate;
    std::fprintf(stderr,
                 "Recovered near-resonance inspiral at a=%.17g, e=%.17g "
                 "(deltachi root gap %.3e)\n",
                 candidate.a, candidate.e, root_gap);
    return true;
}

} // namespace

// ====================== BHMERGER-ZJ: remnant output accessor ======================
// Return a value-copy of the latest remnant result.
BHMerger::RemnantParameters BHMerger::getRemnant() {
    return Remnantoutputs;
}
// ==================== BHMERGER-ZJ: remnant output accessor end ====================

// ======================= BHMERGER-ZJ: Python lifecycle =======================
// Initialize the embedded interpreter on the application's main thread before
// any OpenMP worker can enter calculateRemnant(). Returns 0 on success.
int BHMerger::initializePython() {
    if (Py_IsInitialized() && bhmerger_api_ready.load()) return 0;
    const bool own_interpreter = !Py_IsInitialized();
    if (own_interpreter) Py_Initialize();
    if (!Py_IsInitialized()) {
        std::fprintf(stderr, "Failed to initialize embedded Python\n");
        return -1;
    }
    const PyGILState_STATE gil = PyGILState_Ensure();
    // Another caller may have completed bridge setup while we waited for GIL.
    if (!bhmerger_api_ready.load()) {
        if (_import_array() < 0 || !configurePrecessionPath()) {
            if (PyErr_Occurred()) PyErr_Print();
            std::fprintf(stderr, "Failed to initialize BH Python bridge\n");
            PyGILState_Release(gil);
            if (own_interpreter) Py_FinalizeEx();
            return -1;
        }
        bhmerger_python_owned = own_interpreter;
        bhmerger_api_ready.store(true);
    }
    PyGILState_Release(gil);
    if (own_interpreter) {
        // Only release our own main thread state. The host retains its GIL policy.
        bhmerger_main_thread_state = PyEval_SaveThread();
    }
    return 0;
}

// Finalize only after all integration and OpenMP work has completed. Calling
// this when no interpreter was initialized is intentionally a no-op.
int BHMerger::finalizePython() {
    if (!Py_IsInitialized() || !bhmerger_python_owned) return 0;

    // Reacquire the main interpreter thread state before finalization.
    if (bhmerger_main_thread_state != nullptr) {
        PyEval_RestoreThread(bhmerger_main_thread_state);
        bhmerger_main_thread_state = nullptr;
    }

    bhmerger_api_ready.store(false);
    const int status = Py_FinalizeEx();
    bhmerger_python_owned = false;
    return status;
}
// ===================== BHMERGER-ZJ: Python lifecycle end =====================

namespace {

double norm(const double v[3]) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Produce a reproducible pseudo-random unit vector from the merger phase-space
// values. This avoids global RNG state and keeps OpenMP/MPI results repeatable.
void deterministicUnitVector(const double values[], const int count,
                             const std::uint64_t stream, double direction[3]) {
    std::uint64_t state = 0x9e3779b97f4a7c15ULL ^ stream;
    for (int i = 0; i < count; ++i) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &values[i], sizeof(bits));
        state ^= bits + 0x9e3779b97f4a7c15ULL + (state << 6) + (state >> 2);
    }
    auto uniform = [&state]() {
        state += 0x9e3779b97f4a7c15ULL;
        std::uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
        return (static_cast<double>(z >> 11) + 0.5)
             * (1.0 / 9007199254740992.0);
    };
    const double z = 2.0 * uniform() - 1.0;
    const double phi = 2.0 * std::acos(-1.0) * uniform();
    const double rho = std::sqrt(std::max(0.0, 1.0 - z*z));
    direction[0] = rho * std::cos(phi);
    direction[1] = rho * std::sin(phi);
    direction[2] = z;
}

// 向量归一化
void normalize(double v[3]) {
    double n = norm(v);
    if (n > 0) {
        v[0] /= n; v[1] /= n; v[2] /= n;
    }
}

// 向量叉积: c = a × b
void cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// 向量点积
double dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// 局部坐标系到全局坐标系的转换
void localToGlobal(const double L[3], const double S[3], const double v_local[3], double v_global[3]) {
    double z_axis[3], x_axis[3], y_axis[3];

    // 1. z轴方向就是角动量方向
    std::memcpy(z_axis, L, 3*sizeof(double));
    normalize(z_axis);

    // 2. x轴方向通过施密特正交化得到
    // S在z轴上的分量
    double S_dot_z = dot(S, z_axis);
    double S_perp[3] = { S[0] - S_dot_z*z_axis[0],
                         S[1] - S_dot_z*z_axis[1],
                         S[2] - S_dot_z*z_axis[2] };
    // Zero or aligned spin does not define an azimuth. Pick a deterministic
    // reference direction without perturbing the physical spin magnitude.
    if (norm(S_perp) <= 1.0e-14) {
        const double reference[3] = {
            std::fabs(z_axis[0]) < 0.9 ? 1.0 : 0.0,
            std::fabs(z_axis[0]) < 0.9 ? 0.0 : 1.0,
            0.0
        };
        cross(reference, z_axis, S_perp);
    }
    normalize(S_perp);  // 得到 x 轴
    std::memcpy(x_axis, S_perp, 3*sizeof(double));

    // 3. y轴通过右手系叉积得到
    cross(z_axis, x_axis, y_axis);
    normalize(y_axis);

    // 4. 构造旋转矩阵 R = [x_axis, y_axis, z_axis]
    // 将局部向量转换到全局向量
    // v_global = R * v_local
    v_global[0] = x_axis[0]*v_local[0] + y_axis[0]*v_local[1] + z_axis[0]*v_local[2];
    v_global[1] = x_axis[1]*v_local[0] + y_axis[1]*v_local[1] + z_axis[1]*v_local[2];
    v_global[2] = x_axis[2]*v_local[0] + y_axis[2]*v_local[1] + z_axis[2]*v_local[2];
}

// === 帮助函数：加载 Python 模块（只加载一次） ===
PyObject* getPythonModule() {
    static PyObject* pModule = nullptr;
    if (pModule == nullptr) {
        PyObject *pName = PyUnicode_DecodeFSDefault("precession.eccentricity");
        pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (!pModule) {
            PyErr_Print();
            std::fprintf(stderr, "[Error] Failed to import Python module precession.eccentricity\n");
            return nullptr;
        }
    }
    return pModule;
}

void setDirectMergerSpinDirection(double theta1, double theta2,
                                  double deltaphi, double q,
                                  double chi1, double chi2,
                                  double direction[3]) {
    // Match precession.remnantspindirection() at plunge separation r=10M.
    const double one_plus_q = 1.0 + q;
    const double mass_scale = one_plus_q * one_plus_q;
    const double orbital_angular_momentum = q * std::sqrt(10.0) / mass_scale;
    const double spin1_angular_momentum = chi1 / mass_scale;
    const double spin2_angular_momentum = chi2 * q * q / mass_scale;

    direction[0] =
        spin1_angular_momentum * std::sin(theta1)
        + spin2_angular_momentum * std::sin(theta2) * std::cos(deltaphi);
    direction[1] =
        spin2_angular_momentum * std::sin(theta2) * std::sin(deltaphi);
    direction[2] =
        orbital_angular_momentum
        + spin1_angular_momentum * std::cos(theta1)
        + spin2_angular_momentum * std::cos(theta2);
    normalize(direction);
}

} // namespace


// === 主函数 ===
bool BHMerger::calculateRemnant(double m1,double m2,double dx,double dy,double dz,
                                double vrx,double vry,double vrz,
                                double spin1x,double spin1y,double spin1z,
                                double spin2x,double spin2y,double spin2z,
                                double G,double C)
{
    // ---- Initialize lazily only as a fallback for non-PeTar callers. ----
    // PeTar initializes Python explicitly on its main thread before integration.
    if (initializePython() != 0) return false;
    PyGILState_STATE gil_state = PyGILState_Ensure();
#define CAL_REMNANT_RETURN(value) do { PyGILState_Release(gil_state); return (value); } while (0)

    // ---- 基本物理计算 ----
    double q, m;
    double s1[3], s2[3];
    double r=std::sqrt(dx*dx+dy*dy+dz*dz);
    const double perturbation_seed[] = {
        m1, m2, dx, dy, dz, vrx, vry, vrz,
        spin1x, spin1y, spin1z, spin2x, spin2y, spin2z
    };
    constexpr double spin_perturbation = 1.0e-6;

    if (m1 < m2) {
        q = m1 / m2; m = m1 + m2;
        s1[0]=spin2x;
        s1[1]=spin2y;
        s1[2]=spin2z;
        s2[0]=spin1x;
        s2[1]=spin1y;
        s2[2]=spin1z;
       
    } else {
        q = m2 / m1; m = m1 + m2;
        s1[0]=spin1x;
        s1[1]=spin1y;
        s1[2]=spin1z;
        s2[0]=spin2x;
        s2[1]=spin2y;
        s2[2]=spin2z;
    }
    double tots1 = norm(s1);
    double tots2 = norm(s2);
    if (tots1 == 0.0) {
        deterministicUnitVector(perturbation_seed, 14, 1, s1);
        for (int k = 0; k < 3; ++k) s1[k] *= spin_perturbation;
        tots1 = spin_perturbation;
    }
    if (tots2 == 0.0) {
        deterministicUnitVector(perturbation_seed, 14, 2, s2);
        for (int k = 0; k < 3; ++k) s2[k] *= spin_perturbation;
        tots2 = spin_perturbation;
    }
    // Some fitting expressions are singular exactly at equal mass. Keep the
    // perturbation local to the fitting input; the particle masses are unchanged.
    if (q == 1.0) q = 1.0 - 1.0e-10;

    double meff = m1 * m2 / (m1 + m2);
    //已经翻转了方向，与外部定义一致
    double Lx = meff*(vry*dz - dy*vrz);
    double Ly = meff*(dx*vrz - vrx*dz);
    double Lz = meff*(vrx*dy - dx*vry);
    double Lvec[3]={Lx,Ly,Lz};
    double L  = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    if (!std::isfinite(L)) CAL_REMNANT_RETURN(false);
    if (L == 0.0) {
        deterministicUnitVector(perturbation_seed, 14, 3, Lvec);
        const double angular_momentum_scale =
            std::max(1.0, meff * r
                            * std::sqrt(vrx*vrx + vry*vry + vrz*vrz));
        for (int k = 0; k < 3; ++k) Lvec[k] *= 1.0e-12 * angular_momentum_scale;
        Lx = Lvec[0];
        Ly = Lvec[1];
        Lz = Lvec[2];
        L = norm(Lvec);
    }

    double E = 0.5*meff*(vrx*vrx+vry*vry+vrz*vrz)
             - G*m1*m2/r;

    double clight = C;
    double a0 = G*m1*m2/(2*std::fabs(E))*clight*clight/(G*m);
    double e;
    if (!eccentricityFromEnergyTerm(
            2*E*L*L / (meff*G*G*m1*m1*m2*m2), e)) {
        std::fprintf(stderr,"Invalid eccentricity squared in BH remnant input\n");
        CAL_REMNANT_RETURN(false);
    }
    const bool is_bound_inspiral = (E < 0.0 && e >= 0.0 && e < 1.0);

    double chi1 = tots1;
    double chi2 = tots2;

    // --- 夹角计算（防止 NaN） ---
    #define SAFE_ACOS(x) (std::acos(std::clamp((x), -1.0, 1.0)))
    double theta1 = tots1 > 0.0
                  ? SAFE_ACOS((Lx*s1[0]+Ly*s1[1]+Lz*s1[2])/(L*tots1)) : 0.0;
    double theta2 = tots2 > 0.0
                  ? SAFE_ACOS((Lx*s2[0]+Ly*s2[1]+Lz*s2[2])/(L*tots2)) : 0.0;
    double s1_cross_l[3], s2_cross_l[3];
    cross(s1, Lvec, s1_cross_l);
    cross(s2, Lvec, s2_cross_l);
    const double projected_norm1 = norm(s1_cross_l);
    const double projected_norm2 = norm(s2_cross_l);
    double projected_cross[3];
    cross(s1_cross_l, s2_cross_l, projected_cross);
    // precession.eval_cyclesign uses S1 dot (S2 cross L).
    double deltaphi = (projected_norm1 > 0.0 && projected_norm2 > 0.0 && L > 0.0)
                    ? std::atan2(dot(projected_cross, Lvec)/L,
                                 dot(s1_cross_l, s2_cross_l)) : 0.0;
    //printf("\n~~~~ZJ  debug chi1=%f, chi2=%f, theta1=%f,theta2=%f,deltaphi=%f, a0=%f, e=%f,q=%f\n",chi1,chi2,theta1,theta2,deltaphi,a0,e,q );


    // --- 初始化输出 ---
    Inspiraloutputs = {};
    Remnantoutputs = {};
    Inspiraloutputs.theta1 = theta1;
    Inspiraloutputs.theta2 = theta2;
    Inspiraloutputs.deltaphi = deltaphi;
    Inspiraloutputs.a = a0;
    Inspiraloutputs.e = e;
    Inspiraloutputs.q = q;
    Inspiraloutputs.chi1 = chi1;
    Inspiraloutputs.chi2 = chi2;
    Inspiraloutputs.m = m;
    const InspiralParameters encounter_outputs = Inspiraloutputs;
    bool inspiral_evolution_valid = false;

    // ---- 导入模块 ----
    PyObject *pModule = getPythonModule();
    if (!pModule) CAL_REMNANT_RETURN(false);

    // =====================================================
    // 调用 Python 函数 inspiral_precav
    // =====================================================
    // Peters may already have reached the 10 GM/c^2 fitting boundary.
    // Never evolve the Python inspiral backwards from a smaller separation.
    if (is_bound_inspiral && a0 > 10.0*(1.0+1.0e-8)) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "inspiral_precav");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find inspiral_precav\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        npy_intp dims[1] = {2};
        double a_arr_data[2] = {a0, 10.0};
        PyObject *a_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, a_arr_data);


        PyObject *kwargs = PyDict_New();
        setDoubleKeyword(kwargs, "theta1", theta1);
        setDoubleKeyword(kwargs, "theta2", theta2);
        setDoubleKeyword(kwargs, "deltaphi", deltaphi);
        if (kwargs && a_array) PyDict_SetItemString(kwargs, "a", a_array);
        setDoubleKeyword(kwargs, "e", e);
        setDoubleKeyword(kwargs, "q", q);
        setDoubleKeyword(kwargs, "chi1", chi1);
        setDoubleKeyword(kwargs, "chi2", chi2);

        PyObject *result = callKeywords(pFunc, kwargs);
        Py_XDECREF(kwargs);
        Py_DECREF(pFunc);
        Py_XDECREF(a_array);

        if (!result || !PyDict_Check(result)) {
            Py_XDECREF(result);
            if (PyErr_Occurred()) PyErr_Clear();
            Inspiraloutputs = encounter_outputs;
            inspiral_evolution_valid = recoverResonantInspiral(
                pModule, encounter_outputs, Inspiraloutputs);
            if (!inspiral_evolution_valid)
                std::fprintf(stderr,
                             "Call inspiral_precav failed; using encounter-time "
                             "spin configuration for remnant fit\n");
        }
        else {
        // 从结果字典中提取需要的值
        PyObject *theta1_arr = PyDict_GetItemString(result, "theta1");
        PyObject *theta2_arr = PyDict_GetItemString(result, "theta2"); 
        PyObject *deltaphi_arr = PyDict_GetItemString(result, "deltaphi"); 
        PyObject *deltachi_arr = PyDict_GetItemString(result, "deltachi"); 
        PyObject *kappa_arr = PyDict_GetItemString(result, "kappa"); 
        PyObject *a_arr = PyDict_GetItemString(result, "a"); 
        PyObject *e_arr = PyDict_GetItemString(result, "e"); 
        PyObject *uc_arr = PyDict_GetItemString(result, "uc"); 
        PyObject *u_arr = PyDict_GetItemString(result, "u"); 

        PyObject *chieff_arr = PyDict_GetItemString(result, "chieff");
        PyObject *q_arr      = PyDict_GetItemString(result, "q");
        PyObject *chi1_arr   = PyDict_GetItemString(result, "chi1");
        PyObject *chi2_arr   = PyDict_GetItemString(result, "chi2");
        Inspiraloutputs.theta1 = arrayValue(theta1_arr, true); 
        Inspiraloutputs.theta2 = arrayValue(theta2_arr, true); 
        Inspiraloutputs.deltaphi= arrayValue(deltaphi_arr, true);
        Inspiraloutputs.deltachi= arrayValue(deltachi_arr, true); 
        Inspiraloutputs.kappa = arrayValue(kappa_arr, true); 
        Inspiraloutputs.a = arrayValue(a_arr, true); 
        Inspiraloutputs.e = arrayValue(e_arr, true, e==0.0);
        Inspiraloutputs.uc = arrayValue(uc_arr, true); 
        Inspiraloutputs.u = arrayValue(u_arr, true); 
        // 这几个是一维数组，取第一个元素 
        
        Inspiraloutputs.chieff = arrayValue(chieff_arr, false);
        Inspiraloutputs.q      = arrayValue(q_arr, false);
        Inspiraloutputs.chi1   = arrayValue(chi1_arr, false);
        Inspiraloutputs.chi2   = arrayValue(chi2_arr, false);
        Inspiraloutputs.m      = m;

        Py_DECREF(result);
        if (PyErr_Occurred()) PyErr_Clear();
        inspiral_evolution_valid =
            std::isfinite(Inspiraloutputs.theta1)
            && std::isfinite(Inspiraloutputs.theta2)
            && std::isfinite(Inspiraloutputs.deltaphi)
            && std::isfinite(Inspiraloutputs.a)
            && std::isfinite(Inspiraloutputs.e)
            && std::isfinite(Inspiraloutputs.q)
            && std::isfinite(Inspiraloutputs.chi1)
            && std::isfinite(Inspiraloutputs.chi2);
        if (!inspiral_evolution_valid) {
            Inspiraloutputs = encounter_outputs;
            inspiral_evolution_valid = recoverResonantInspiral(
                pModule, encounter_outputs, Inspiraloutputs);
            if (!inspiral_evolution_valid)
                std::fprintf(stderr,
                             "Non-finite inspiral output; using encounter-time "
                             "spin configuration for remnant fit\n");
        }
        }
    }
    if (!is_bound_inspiral || !inspiral_evolution_valid) {
        // ================= BHMERGER-ZJ: hyperbolic/direct merger path =================
        // inspiral_precav() supports only 0 <= e < 1. Keep the encounter-time
        // spin angles and proceed directly to the remnant fits for e >= 1.
        setDirectMergerSpinDirection(Inspiraloutputs.theta1,
                                     Inspiraloutputs.theta2,
                                     Inspiraloutputs.deltaphi,
                                     Inspiraloutputs.q,
                                     Inspiraloutputs.chi1,
                                     Inspiraloutputs.chi2,
                                     Remnantoutputs.spindirection);
        // =============== BHMERGER-ZJ: hyperbolic/direct merger path end ===============
    }

    // =====================================================
    // remnantmass
    // =====================================================
    {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "remnantmass");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find remnantmass\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        PyObject *kwargs = PyDict_New();
        setDoubleKeyword(kwargs, "theta1", Inspiraloutputs.theta1);
        setDoubleKeyword(kwargs, "theta2", Inspiraloutputs.theta2);
        setDoubleKeyword(kwargs, "q", Inspiraloutputs.q);
        setDoubleKeyword(kwargs, "chi1", Inspiraloutputs.chi1);
        setDoubleKeyword(kwargs, "chi2", Inspiraloutputs.chi2);

        PyObject *result = callKeywords(pFunc, kwargs);
        Py_XDECREF(kwargs);
        Py_DECREF(pFunc);

        if (validArray(result, 1)) {
            PyArrayObject *arr = (PyArrayObject*)result;
            Remnantoutputs.mfin = *((double*)PyArray_DATA(arr));
        } else {
            PyErr_Print(); fprintf(stderr, "remnantmass call failed\n");
            Py_XDECREF(result);
            CAL_REMNANT_RETURN(false);
        }
        Py_XDECREF(result);
    }

    // =====================================================
    // remnantspin
    // =====================================================
    {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "remnantspin");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find remnantspin\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        PyObject *kwargs = PyDict_New();
        setDoubleKeyword(kwargs, "theta1", Inspiraloutputs.theta1);
        setDoubleKeyword(kwargs, "theta2", Inspiraloutputs.theta2);
        setDoubleKeyword(kwargs, "deltaphi", Inspiraloutputs.deltaphi);
        setDoubleKeyword(kwargs, "q", Inspiraloutputs.q);
        setDoubleKeyword(kwargs, "chi1", Inspiraloutputs.chi1);
        setDoubleKeyword(kwargs, "chi2", Inspiraloutputs.chi2);

        PyObject *result = callKeywords(pFunc, kwargs);
        Py_XDECREF(kwargs);
        Py_DECREF(pFunc);

        if (validArray(result, 1)) {
            PyArrayObject *arr = (PyArrayObject*)result;
            Remnantoutputs.spin = *((double*)PyArray_DATA(arr));
            if (!std::isfinite(Remnantoutputs.spin)
                || Remnantoutputs.spin < 0.0
                || Remnantoutputs.spin > 1.0) {
                std::fprintf(stderr, "Invalid remnant spin magnitude: %.17g\n",
                             Remnantoutputs.spin);
                Py_XDECREF(result);
                CAL_REMNANT_RETURN(false);
            }
        } else {
            PyErr_Print(); fprintf(stderr, "remnantspin call failed\n");
            Py_XDECREF(result);
            CAL_REMNANT_RETURN(false);
        }
        Py_XDECREF(result);
    }

    // =====================================================
    // remnantspindirection
    // =====================================================
    if (is_bound_inspiral && inspiral_evolution_valid) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "remnantspindirection");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find remnantspindirection\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        PyObject *kwargs = PyDict_New();
        setDoubleKeyword(kwargs, "theta1", Inspiraloutputs.theta1);
        setDoubleKeyword(kwargs, "theta2", Inspiraloutputs.theta2);
        setDoubleKeyword(kwargs, "deltaphi", Inspiraloutputs.deltaphi);
        setDoubleKeyword(kwargs, "a", Inspiraloutputs.a);
        setDoubleKeyword(kwargs, "e", Inspiraloutputs.e);
        setDoubleKeyword(kwargs, "q", Inspiraloutputs.q);
        setDoubleKeyword(kwargs, "chi1", Inspiraloutputs.chi1);
        setDoubleKeyword(kwargs, "chi2", Inspiraloutputs.chi2);

        PyObject *result = callKeywords(pFunc, kwargs);
        Py_XDECREF(kwargs);
        Py_DECREF(pFunc);

        if (validArray(result, 3)) {
            PyArrayObject *arr = (PyArrayObject*)result;
            double *data = (double*)PyArray_DATA(arr);
            Remnantoutputs.spindirection[0] = data[0];
            Remnantoutputs.spindirection[1] = data[1];
            Remnantoutputs.spindirection[2] = data[2];
            const double direction_norm = norm(Remnantoutputs.spindirection);
            if (!std::isfinite(direction_norm) || direction_norm == 0.0) {
                std::fprintf(stderr, "Invalid remnant spin direction\n");
                Py_XDECREF(result);
                CAL_REMNANT_RETURN(false);
            }
        } else {
            PyErr_Print(); fprintf(stderr, "remnantspindirection failed\n");
            Py_XDECREF(result);
            CAL_REMNANT_RETURN(false);
        }
        Py_XDECREF(result);
    }

    {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "remnantkick");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find remnantkick\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        PyObject *kwargs = PyDict_New();
        setDoubleKeyword(kwargs, "theta1", Inspiraloutputs.theta1);
        setDoubleKeyword(kwargs, "theta2", Inspiraloutputs.theta2);
        setDoubleKeyword(kwargs, "deltaphi", Inspiraloutputs.deltaphi);

        setDoubleKeyword(kwargs, "q", Inspiraloutputs.q);
        setDoubleKeyword(kwargs, "chi1", Inspiraloutputs.chi1);
        setDoubleKeyword(kwargs, "chi2", Inspiraloutputs.chi2);

        PyObject *result = callKeywords(pFunc, kwargs);
        Py_XDECREF(kwargs);
        Py_DECREF(pFunc);

        if (validArray(result, 4)) {
            PyArrayObject *arr = (PyArrayObject*)result;
            double *data = (double*)PyArray_DATA(arr);
            Remnantoutputs.vkick[0] = data[1];
            Remnantoutputs.vkick[1] = data[2];
            Remnantoutputs.vkick[2] = data[3];
        } else {
            PyErr_Print(); fprintf(stderr, "remnantkick failed\n");
            Py_XDECREF(result);
            CAL_REMNANT_RETURN(false);
        }
        Py_XDECREF(result);
    }
    localToGlobal(Lvec, s1, Remnantoutputs.vkick, Remnantoutputs.vkick_nor);
    localToGlobal(Lvec, s1, Remnantoutputs.spindirection, Remnantoutputs.spindirection_nor);

    // Python remains alive for later mergers and is finalized after the main loop.
    PyGILState_Release(gil_state);
#undef CAL_REMNANT_RETURN
    return true;
}





