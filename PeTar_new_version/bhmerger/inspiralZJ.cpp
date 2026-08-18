#include <Python.h>
#include <numpy/arrayobject.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "inspiralZJ.hpp"

namespace {

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
PyThreadState* bhmerger_main_thread_state = nullptr;

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
    if (Py_IsInitialized()) return 0;

    Py_Initialize();
    if (!Py_IsInitialized()) {
        std::fprintf(stderr, "Failed to initialize embedded Python\n");
        return -1;
    }
    bhmerger_python_owned = true;

    if (_import_array() < 0) {
        PyErr_Print();
        std::fprintf(stderr, "Failed to initialize NumPy C API\n");
        Py_FinalizeEx();
        bhmerger_python_owned = false;
        return -1;
    }

    if (PyRun_SimpleString(
            "import sys; "
            "sys.path.insert(0, './bhmerger/precession-master'); "
            "sys.path.insert(0, 'bhmerger/precession-master')") != 0) {
        PyErr_Print();
        std::fprintf(stderr, "Failed to configure embedded Python module path\n");
        Py_FinalizeEx();
        bhmerger_python_owned = false;
        return -1;
    }

    // Release the GIL so OpenMP worker threads can acquire it around Python calls.
    bhmerger_main_thread_state = PyEval_SaveThread();
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
    double e = std::sqrt(1 + 2*E*L*L / (meff*G*G*m1*m1*m2*m2));
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
    double deltaphi = (projected_norm1 > 0.0 && projected_norm2 > 0.0)
                    ? SAFE_ACOS(dot(s1_cross_l, s2_cross_l)
                                /(projected_norm1*projected_norm2)) : 0.0;
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

    // ---- 导入模块 ----
    PyObject *pModule = getPythonModule();
    if (!pModule) CAL_REMNANT_RETURN(false);

    // =====================================================
    // 调用 Python 函数 inspiral_precav
    // =====================================================
    if (is_bound_inspiral) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "inspiral_precav");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find inspiral_precav\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        npy_intp dims[1] = {2};
        double a_arr_data[2] = {a0, 10.0};
        PyObject *a_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, a_arr_data);
        Py_INCREF(a_array);

        PyObject *kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "theta1", PyFloat_FromDouble(theta1));
        PyDict_SetItemString(kwargs, "theta2", PyFloat_FromDouble(theta2));
        PyDict_SetItemString(kwargs, "deltaphi", PyFloat_FromDouble(deltaphi));
        PyDict_SetItemString(kwargs, "a", a_array);
        PyDict_SetItemString(kwargs, "e", PyFloat_FromDouble(e));
        PyDict_SetItemString(kwargs, "q", PyFloat_FromDouble(q));
        PyDict_SetItemString(kwargs, "chi1", PyFloat_FromDouble(chi1));
        PyDict_SetItemString(kwargs, "chi2", PyFloat_FromDouble(chi2));

        PyObject *result = PyObject_Call(pFunc, PyTuple_New(0), kwargs);
        Py_DECREF(kwargs);
        Py_DECREF(pFunc);
        Py_DECREF(a_array);

        if (!result) {
            PyErr_Print(); fprintf(stderr, "Call inspiral_precav failed\n");
            CAL_REMNANT_RETURN(false);
        }

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
        Inspiraloutputs.theta1 = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)theta1_arr, PyArray_GETPTR2((PyArrayObject*)theta1_arr,0,1))); 
        Inspiraloutputs.theta2 = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)theta2_arr, PyArray_GETPTR2((PyArrayObject*)theta2_arr,0,1))); 
        Inspiraloutputs.deltaphi= PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)deltaphi_arr,PyArray_GETPTR2((PyArrayObject*)deltaphi_arr,0,1)));
        Inspiraloutputs.deltachi= PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)deltachi_arr,PyArray_GETPTR2((PyArrayObject*)deltachi_arr,0,1))); 
        Inspiraloutputs.kappa = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)kappa_arr, PyArray_GETPTR2((PyArrayObject*)kappa_arr,0,1))); 
        Inspiraloutputs.a = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)a_arr, PyArray_GETPTR2((PyArrayObject*)a_arr,0,1))); 
        Inspiraloutputs.e = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)e_arr, PyArray_GETPTR2((PyArrayObject*)e_arr,0,1))); 
        Inspiraloutputs.uc = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)uc_arr, PyArray_GETPTR2((PyArrayObject*)uc_arr,0,1))); 
        Inspiraloutputs.u = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)u_arr, PyArray_GETPTR2((PyArrayObject*)u_arr,0,1))); 
        // 这几个是一维数组，取第一个元素 
        
        Inspiraloutputs.chieff = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)chieff_arr, PyArray_GETPTR1((PyArrayObject*)chieff_arr,0)));
        Inspiraloutputs.q      = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)q_arr, PyArray_GETPTR1((PyArrayObject*)q_arr,0)));
        Inspiraloutputs.chi1   = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)chi1_arr, PyArray_GETPTR1((PyArrayObject*)chi1_arr,0)));
        Inspiraloutputs.chi2   = PyFloat_AsDouble(PyArray_GETITEM((PyArrayObject*)chi2_arr, PyArray_GETPTR1((PyArrayObject*)chi2_arr,0)));
        Inspiraloutputs.m      = m;

        Py_DECREF(result);
    }
    else {
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
        PyDict_SetItemString(kwargs, "theta1", PyFloat_FromDouble(Inspiraloutputs.theta1));
        PyDict_SetItemString(kwargs, "theta2", PyFloat_FromDouble(Inspiraloutputs.theta2));
        PyDict_SetItemString(kwargs, "q", PyFloat_FromDouble(Inspiraloutputs.q));
        PyDict_SetItemString(kwargs, "chi1", PyFloat_FromDouble(Inspiraloutputs.chi1));
        PyDict_SetItemString(kwargs, "chi2", PyFloat_FromDouble(Inspiraloutputs.chi2));

        PyObject *result = PyObject_Call(pFunc, PyTuple_New(0), kwargs);
        Py_DECREF(kwargs);
        Py_DECREF(pFunc);

        if (PyArray_Check(result)) {
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
        PyDict_SetItemString(kwargs, "theta1", PyFloat_FromDouble(Inspiraloutputs.theta1));
        PyDict_SetItemString(kwargs, "theta2", PyFloat_FromDouble(Inspiraloutputs.theta2));
        PyDict_SetItemString(kwargs, "deltaphi", PyFloat_FromDouble(Inspiraloutputs.deltaphi));
        PyDict_SetItemString(kwargs, "q", PyFloat_FromDouble(Inspiraloutputs.q));
        PyDict_SetItemString(kwargs, "chi1", PyFloat_FromDouble(Inspiraloutputs.chi1));
        PyDict_SetItemString(kwargs, "chi2", PyFloat_FromDouble(Inspiraloutputs.chi2));

        PyObject *result = PyObject_Call(pFunc, PyTuple_New(0), kwargs);
        Py_DECREF(kwargs);
        Py_DECREF(pFunc);

        if (PyArray_Check(result)) {
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
    if (is_bound_inspiral) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "remnantspindirection");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Print(); fprintf(stderr, "Cannot find remnantspindirection\n");
            Py_XDECREF(pFunc); CAL_REMNANT_RETURN(false);
        }

        PyObject *kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "theta1", PyFloat_FromDouble(Inspiraloutputs.theta1));
        PyDict_SetItemString(kwargs, "theta2", PyFloat_FromDouble(Inspiraloutputs.theta2));
        PyDict_SetItemString(kwargs, "deltaphi", PyFloat_FromDouble(Inspiraloutputs.deltaphi));
        PyDict_SetItemString(kwargs, "a", PyFloat_FromDouble(Inspiraloutputs.a));
        PyDict_SetItemString(kwargs, "e", PyFloat_FromDouble(Inspiraloutputs.e));
        PyDict_SetItemString(kwargs, "q", PyFloat_FromDouble(Inspiraloutputs.q));
        PyDict_SetItemString(kwargs, "chi1", PyFloat_FromDouble(Inspiraloutputs.chi1));
        PyDict_SetItemString(kwargs, "chi2", PyFloat_FromDouble(Inspiraloutputs.chi2));

        PyObject *result = PyObject_Call(pFunc, PyTuple_New(0), kwargs);
        Py_DECREF(kwargs);
        Py_DECREF(pFunc);

        if (PyArray_Check(result)) {
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
        PyDict_SetItemString(kwargs, "theta1", PyFloat_FromDouble(Inspiraloutputs.theta1));
        PyDict_SetItemString(kwargs, "theta2", PyFloat_FromDouble(Inspiraloutputs.theta2));
        PyDict_SetItemString(kwargs, "deltaphi", PyFloat_FromDouble(Inspiraloutputs.deltaphi));

        PyDict_SetItemString(kwargs, "q", PyFloat_FromDouble(Inspiraloutputs.q));
        PyDict_SetItemString(kwargs, "chi1", PyFloat_FromDouble(Inspiraloutputs.chi1));
        PyDict_SetItemString(kwargs, "chi2", PyFloat_FromDouble(Inspiraloutputs.chi2));

        PyObject *result = PyObject_Call(pFunc, PyTuple_New(0), kwargs);
        Py_DECREF(kwargs);
        Py_DECREF(pFunc);

        if (PyArray_Check(result)) {
            PyArrayObject *arr = (PyArrayObject*)result;
            double *data = (double*)PyArray_DATA(arr);
            Remnantoutputs.vkick[0] = data[1];
            Remnantoutputs.vkick[1] = data[2];
            Remnantoutputs.vkick[2] = data[3];
        } else {
            PyErr_Print(); fprintf(stderr, "remnantspindirection failed\n");
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





