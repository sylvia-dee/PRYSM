import pytest

import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np

from psm.icecore import archive


def test_archive_diffusivity_one():
    
    rho = np.array([1])
    D = archive.diffusivity(rho)
    
    assert len(D) == 1
    assert np.round(D, 9) == np.round(1.2240425e-08, 9)
    assert D > 0


def test_archive_diffusivity_ones():

    rho = np.array([1, 1, 1, 1, 1])
    D = archive.diffusivity(rho)

    assert len(D) == 5
    assert np.round(D[1], 9) == np.round(np.mean(D), 9)


# def test_archive_diffusivity_natural():    
#     rho = np.array([1,2,3])
#     D = archive.diffusivity(rho)
