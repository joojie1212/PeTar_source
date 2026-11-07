#!/usr/bin/env python3
import os
import sys

import numpy as np


repository = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(repository, "bhmerger", "precession-master"))

from precession.eccentricity import inspiral_precav


np.random.seed(7)
result = inspiral_precav(
    theta1=1.0,
    theta2=1.2,
    deltaphi=-0.5,
    a=np.array([1000.0, 900.0]),
    e=0.1,
    q=0.8,
    chi1=0.5,
    chi2=0.4,
    requested_outputs=["deltaphi"],
)
phase = np.atleast_1d(result["deltaphi"])
if not np.all(phase <= 0):
    raise AssertionError(f"negative input phase changed branch: {phase}")

print("precession cyclesign regression passed")

from precession.precession import _arccos_roundoff, eval_theta1
for sign in [-1.0, 1.0]:
    x = np.nextafter(sign, sign * np.inf)
    assert np.isfinite(_arccos_roundoff(x))
    # q=1, chi1=1 makes costheta1 = chieff + deltachi.
    angle = eval_theta1(deltachi=0., chieff=x, q=1., chi1=1.)
    assert np.all(np.isfinite(angle))
    assert np.allclose(angle, 0. if sign>0 else np.pi, atol=0., rtol=0.)
for value in [1.01, -1.01, np.nan, np.inf]:
    try:
        _arccos_roundoff(value)
    except ValueError:
        pass
    else:
        raise AssertionError(f"invalid arccos input accepted: {value}")
print("precession endpoint-roundoff regressions passed")

from precession.precession import (eval_S, eval_S1, eval_S2, eval_q,
    eval_kappa, deltachilimits_definition, deltachitildeav2)
assert eval_q(2., 2.) == 1.
q, chi1, chi2 = .8, .5, .4
s1, s2 = eval_S1(q, chi1), eval_S2(q, chi2)
assert np.allclose(eval_S(theta1=0., theta2=0., deltaphi=0.,
                        q=q, chi1=chi1, chi2=chi2), s1+s2)
assert np.allclose(eval_S(theta1=0., theta2=np.pi, deltaphi=0.,
                        q=q, chi1=chi1, chi2=chi2), abs(s1-s2))
for fn in [eval_S, eval_kappa]:
    try:
        fn()
    except TypeError:
        pass
    else:
        raise AssertionError('missing parameters accepted')
limits=deltachilimits_definition(q,chi1,chi2).ravel()
expected=(chi1+q*chi2)/(1+q)
assert np.allclose(limits,[-expected,expected])
angles=np.linspace(0,np.pi,101)
difference=(chi1*np.cos(angles[:,None])-q*chi2*np.cos(angles))/(1+q)
assert difference.min() >= limits[0] and difference.max() <= limits[1]

import mpmath as mp
mp.mp.dps=100
samples=np.r_[0.,np.logspace(-30,-2,40),.009999,.01,.010001,.1,.5,.99]
reference=[]
for value in samples:
    m=mp.mpf(float(value))
    reference.append(float(mp.mpf(3)/8 if m==0 else
        (2+m-2*(1+m)*mp.ellipe(m)/mp.ellipk(m))/(3*m*m)))
actual=deltachitildeav2(samples)
error=np.max(np.abs(actual-reference))
assert error < 5e-12,error
assert deltachitildeav2(0.) == .375
print('BH-21/22/23 regressions passed; max coefficient absolute error:',error)
