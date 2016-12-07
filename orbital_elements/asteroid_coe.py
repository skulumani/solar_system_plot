import numpy as np
from keplerian_orbit.keplerian_orbit import tof_delta_t, kepler_eq_E
from utilities.attitude import normalize

def asteroid_epoch(ast_flag):
    """
        This holds the orbital elements for the asteroids as taken from JPL

        Return the state at the JD epoch for use in a propogate function
    """

    if ast_flag == 0:
        # 2008 EV5
        a = 0.9582899238313918
        ecc = 0.08348599378460778    
        p =  a * (1-ecc**2)
        inc = np.deg2rad(7.436787362690259)
        raan = np.deg2rad(93.39122898787916)
        argp = np.deg2rad(234.8245876826614)
        M = np.deg2rad(3.409187469072454)
        JD_epoch = 2457800.5
    elif ast_flag == 1:
        # Itokawa
        a = 1.324163617639197
        ecc = 0.28011765678781    
        p =  a * (1-ecc**2)
        inc = np.deg2rad(1.62145641293925)
        raan = np.deg2rad(69.07992986350325)
        argp = np.deg2rad(162.8034822691509)
        M = np.deg2rad(131.4340297670125)
        JD_epoch = 2457800.5
    elif ast_flag == 2:
        # Bennu
        a = 1.126391026007489
        ecc = 0.2037451112033579    
        p =  a * (1-ecc**2)
        inc = np.deg2rad(6.034939195483961)
        raan = np.deg2rad(2.060867837066797)
        argp = np.deg2rad(66.22306857848962)
        M = np.deg2rad(101.7039476994243)
        JD_epoch = 2455562.5
    else:
        print("No such asteroid defined yet. Using Itokawa so nothing breaks instead")

        # Itokawa
        a = 1.324163617639197
        ecc = 0.28011765678781    
        p =  a * (1-ecc**2)
        inc = np.deg2rad(1.62145641293925)
        raan = np.deg2rad(69.07992986350325)
        argp = np.deg2rad(162.8034822691509)
        M = np.deg2rad(131.4340297670125)
        JD_epoch = 2457800.5

    # convert M to nu for later and output the coe
    E, nu, count = kepler_eq_E(M,ecc)

    epoch = (p,ecc,inc,raan,argp,normalize(nu,0,2*np.pi),JD_epoch)

    return epoch

def asteroid_coe(JD_curr,ast_flag):
    """
        Output the current COE for the chosen asteroid
    """

    # load the asteroid COE at the epoch
    (p,ecc,inc,raan,argp,nu_0, JD_epoch) = asteroid_epoch(ast_flag)

    # compute the delta t
    delta_t = (JD_curr - JD_epoch) * 86400
    mu = 1.32712440018e20 # m^3 / s^2
    mu = 1/149597870700**3 * mu # au^3 / sec^2

    # propogate to the current JD_curr
    (E_f, M_f, nu_f) = tof_delta_t(p,ecc,mu,nu_0,delta_t)

    # output the current COE
    coe = (p,ecc,inc,raan,argp,normalize(nu_f,0,2*np.pi))

    return coe