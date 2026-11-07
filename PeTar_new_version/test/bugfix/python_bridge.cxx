// Include the implementation to exercise its actual Python ownership helpers.
#include "../../bhmerger/inspiralZJ.cpp"
#include <cassert>
#include <iostream>

void python(const char* code) {
    const auto gil=PyGILState_Ensure();
    assert(PyRun_SimpleString(code)==0);
    PyGILState_Release(gil);
}
bool calculate(double spin_y=0.5, bool bound=true) {
    return BHMerger::calculateRemnant(2.,1.,10.,0.,0.,0.,bound ? -0.5 : -2.,0.,
                                     0.5,0.,0.,0.,spin_y,0.,1.,10.);
}
int main(int argc, char** argv) {
    const bool host=argc>1 && std::string(argv[1])=="--host";
    if(host) Py_Initialize();
    double e;
    assert(eccentricityFromEnergyTerm(-1.-std::numeric_limits<double>::epsilon(),e) && e==0.);
    assert(!eccentricityFromEnergyTerm(-1.01,e));
    assert(!eccentricityFromEnergyTerm(INFINITY,e));
    assert(BHMerger::initializePython()==0);
    if(argc>1 && std::string(argv[1])=="--import-only") {
        python("import precession.eccentricity as real_module; print(real_module.__file__)");
        assert(BHMerger::finalizePython()==0);
        return 0;
    }
    python(R"PY(
import sys, types, numpy as np, gc, tracemalloc
pkg=types.ModuleType('precession'); mod=types.ModuleType('precession.eccentricity')
sys.modules['precession']=pkg; sys.modules['precession.eccentricity']=mod
mod.fail=''; mod.bad=''; mod.last_phi=0.; mod.circular=False
mod.fail_inspiral=False; mod.near_resonance=True; mod.resonance_calls=0
def fit(name,size):
    def call(**kw):
        mod.last_phi=kw.get('deltaphi',mod.last_phi)
        if mod.fail==name: raise RuntimeError('injected '+name)
        if mod.bad==name: return np.array([],dtype=float)
        return np.array([.9]) if size==1 else np.array([1.,0.,0.]) if size==3 else np.array([1.,1.,0.,0.])
    return call
for name,size in [('remnantmass',1),('remnantspin',1),('remnantspindirection',3),('remnantkick',4)]:
    setattr(mod,name,fit(name,size))
def inspiral(**kw):
    if mod.fail_inspiral and 'requested_outputs' not in kw:
        raise OverflowError('Range exceeds valid bounds')
    if 'requested_outputs' in kw:
        root_plus=0.200000000001 if mod.near_resonance else 0.6
        root_three=0.7 if mod.near_resonance else 1.e12
        result={k:np.array([[v,v]]) for k,v in [
            ('kappa',0.1),('a',10.),('e',0.),('uc',0.2),('u',0.2),
            ('deltachiminus',0.2),('deltachiplus',root_plus),('deltachi3',root_three)]}
        result.update({k:np.array([v]) for k,v in [
            ('chieff',0.05),('q',kw['q']),('chi1',kw['chi1']),('chi2',kw['chi2'])]})
        return result
    result={k:np.array([[v,v]]) for k,v in [('theta1',kw['theta1']),('theta2',kw['theta2']),('deltaphi',kw['deltaphi']),('deltachi',0.),('kappa',0.),('a',10.),('e',0.),('uc',0.),('u',0.)]}
    result.update({k:np.array([v]) for k,v in [('chieff',0.),('q',kw['q']),('chi1',kw['chi1']),('chi2',kw['chi2'])]})
    if mod.circular:
        assert kw['e']==0.
        result['e']=np.array([[0.]])
        result['deltaphi']=np.array([[kw['deltaphi'],0.375]])
    return result
mod.inspiral_precav=inspiral
def resonance(**kw):
    mod.resonance_calls += 1
    return np.array([0.2000000000005])
mod.deltachiresonance=resonance
def to_angles(**kw):
    return np.array([[0.8],[1.1],[0.25*np.sign(kw['cyclesign'])]])
mod.conserved_to_angles=to_angles
)PY");
    assert(calculate(0.5,false)); python("assert mod.last_phi>0");
    assert(calculate(-0.5,false)); python("assert mod.last_phi<0");
    assert(calculate());
    python("mod.circular=True");
    assert(BHMerger::calculateRemnant(2.,1.,3.,0.,0.,0.,-1.,0.,
                                     0.5,0.,0.,0.,0.5,0.,1.,10.));
    python("assert mod.last_phi==0.375; mod.circular=False");
    python("mod.fail_inspiral=True; mod.near_resonance=True; mod.resonance_calls=0");
    assert(calculate());
    python("assert mod.resonance_calls==1; assert abs(mod.last_phi-0.25)<1e-15");
    python("mod.near_resonance=False; mod.resonance_calls=0");
    assert(calculate());
    python("assert mod.resonance_calls==0; mod.fail_inspiral=False; mod.near_resonance=True");
    for(const char* name: {"remnantmass","remnantspin","remnantspindirection","remnantkick"}) {
        std::string code="mod.fail='"; code+=name; code+="'"; python(code.c_str());
        assert(!calculate()); python("mod.fail=''");
        code="mod.bad='"; code+=name; code+="'"; python(code.c_str());
        assert(!calculate()); python("mod.bad=''");
    }
    python("tracemalloc.start()");
    for(int i=0;i<100;++i) assert(calculate());
    python("gc.collect(); baseline=tracemalloc.get_traced_memory()[0]");
    for(int i=0;i<2000;++i) assert(calculate());
    python("gc.collect(); growth=tracemalloc.get_traced_memory()[0]-baseline; print('Python retained-memory growth:',growth); assert growth<65536; tracemalloc.stop()");
    assert(BHMerger::finalizePython()==0);
    if(host) {
        assert(Py_IsInitialized());
        python("assert 2+2==4");
        assert(Py_FinalizeEx()==0);
    }
    std::cout << "Python bridge regressions passed\n";
}
