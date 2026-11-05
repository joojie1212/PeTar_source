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
