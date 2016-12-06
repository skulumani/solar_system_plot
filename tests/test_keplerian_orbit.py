"""Pytest for keplerian_orbit.py"""
from keplerian_orbit.keplerian_orbit import kepler_eq_E
import numpy as np

def test_kepler_eq_E():
    """
        A series of cases run in Matlab and copied here
    """
    M = np.deg2rad(110)
    ecc = 0.9
    E_matlab = 2.475786297687611 
    nu_matlab = 2.983273149717047

    E_python, nu_python, count_python = kepler_eq_E(M,ecc)

    np.testing.assert_allclose((E_python, nu_python),(E_matlab,nu_matlab))
    
