import numpy as np
from keplerian_orbit.coe import coe2rv

def test_coe2rv_equatorial_circular():
    """Test COE to RV for equatorial circular orbit around Earth"""

    p = 6378.137 # km
    ecc = 0.0
    inc = 0.0
    raan = 0.0
    arg_p = 0.0
    nu = 0.0
    mu = 398600.5 # km^3 /sec^2

    R_ijk_true = np.array([6378.137,0,0])
    V_ijk_true = np.sqrt(mu/p) * np.array([0,1,0])
    R_pqw_true = R_ijk_true
    V_pqw_true = V_ijk_true

    R_ijk, V_ijk, R_pqw, V_pqw = coe2rv(p,ecc,inc,raan,arg_p,nu, mu)

    np.testing.assert_array_almost_equal(R_ijk_true,R_ijk)
    np.testing.assert_array_almost_equal(V_ijk_true,V_ijk)
    np.testing.assert_array_almost_equal(V_pqw_true,V_pqw)
    np.testing.assert_array_almost_equal(V_pqw_true,V_pqw)

def test_coe2rv_polar_circular():
    """Test COE to RV for polar circular orbit around Earth"""

    p = 6378.137 # km
    ecc = 0.0
    inc = np.pi/2
    raan = 0.0
    arg_p = 0.0
    nu = 0.0
    mu = 398600.5 # km^3 /sec^2

    R_ijk_true = np.array([6378.137,0,0])
    V_ijk_true = np.sqrt(mu/p) * np.array([0.0,0.0,1])

    R_ijk, V_ijk, R_pqw, V_pqw = coe2rv(p,ecc,inc,raan,arg_p,nu, mu)

    np.testing.assert_array_almost_equal(R_ijk_true,R_ijk)
    np.testing.assert_array_almost_equal(V_ijk_true,V_ijk)

def test_coe2rv_equatorial_circular_quarter():
    """Test COE to RV for equatorial circular orbit around Earth"""

    p = 6378.137 # km
    ecc = 0.0
    inc = 0.0
    raan = 0.0
    arg_p = 0.0
    nu = np.pi/2
    mu = 398600.5 # km^3 /sec^2

    R_ijk_true = np.array([0,6378.137,0])
    V_ijk_true = np.sqrt(mu/p) * np.array([-1.0,0.0,0])

    R_ijk, V_ijk, R_pqw, V_pqw = coe2rv(p,ecc,inc,raan,arg_p,nu, mu)

    np.testing.assert_array_almost_equal(R_ijk_true,R_ijk)
    np.testing.assert_array_almost_equal(V_ijk_true,V_ijk)    

def test_coe2rv_equatorial_circular_half():
    """Test COE to RV for equatorial circular orbit around Earth"""

    p = 6378.137 # km
    ecc = 0.0
    inc = 0.0
    raan = 0.0
    arg_p = 0.0
    nu = np.pi
    mu = 398600.5 # km^3 /sec^2

    R_ijk_true = np.array([-p,0,0])
    V_ijk_true = np.sqrt(mu/p) * np.array([0.0,-1.0,0])

    R_ijk, V_ijk, R_pqw, V_pqw = coe2rv(p,ecc,inc,raan,arg_p,nu, mu)

    np.testing.assert_array_almost_equal(R_ijk_true,R_ijk)
    np.testing.assert_array_almost_equal(V_ijk_true,V_ijk) 