import numpy as np
from utilities.attitude import ROT1, ROT2, ROT3

def coe2rv(p,ecc,inc,raan,arg_p,nu,mu):
    """
       Purpose:
           - Convert the classical orbital elements (COEs) to position (R)
               and velocity (V) vectors.

       [R_ijk,V_ijk,R_pqw,V_pqw] = coe2rv(p,ecc,inc,raan,arg_p,nu,mu,arglat,truelon,lonper )

       Inputs:
           - p - semi-major axis (km)
           - ecc - eccentricity
           - raan - right acsension of the ascending node (rad) 0 < raan <2*pi
           - inc - inclination (rad) 0 < inc < pi
           - arg_p - argument of periapsis (rad) 0 < arg_p < 2*pi
           - nu - true anomaly (rad) 0 < nu < 2*pi
           - mu - gravitational parameter of central body (km^3/sec^2).
           - arglat - argument of latitude(CI) rad 0 <arglat < 2*pi
           - truelon - true longitude (CE) rad 0 < truelon < 2*pi
           - lonper - longitude of periapsis rad 0 < lonper < 2*pi

       Outpus:
           - R_ijk - position vector in inertial frame (km)
           - V_ijk - velocity vector in inertial frame (km/sec)
           - R_pqw - position vecotr in perifocal frame (km)
           - V_pqw - velocity vector in perifocal frame (km/sec)

       Dependencies:
           - ROT1 - elementary rotation about first axis
           - ROT2 - elementary rotation about second axis
           - ROT3 - elementary rotation about third axis

       Author:
           - Shankar Kulumani 18 Aug 2012
               - revised from code written at USAFA Fall 2007
           - Shankar Kulumani 18 Sept 2012
               - modified rotation matrix from PQW to ECI
           - Shankar Kulumani 29 Sept 2012
               - modified input to take semi-latus rectum to allow calculation
               for all orbit types
           - Shankar Kulumani 7 Dec 2014
               - loop for vector inputs
           - Shankar Kulumani 7 Dec 2016
               - Convert to python and remove vectorized inputs

    """

    tol = 1e-9

    # check eccentricity for special cases
    if ( ecc < tol ):
        # circular equatorial
        if (inc < tol) or ( abs(inc - np.pi)< tol ):
            arg_p = 0.0
            raan= 0.0
            # nu   = truelon
            nu = raan + arg_p + nu
        else:
            # circular inclined
            arg_p = 0.0
            nu = arg_p + nu
            # nu = arglat
        
    elif ( ( inc <tol) or (abs(inc-pi)<tol) ):    # elliptical equatorial
        arg_p = raan + arg_p
        raan = 0.0
        # arg_p=lonper
    else:
        print("Something crazy happened")
    
    
    cosnu = np.cos(nu)
    sinnu = np.sin(nu)
    
    radius = p/(1+ecc*cosnu)
    
    # semi latus rectum check
    if ( abs(p) < 0.0001):
        p = 0.0001
    
    
    # calculate postion and velocity in perifocal frame
    R_pqw = radius*np.array([cosnu,sinnu,0])
    V_pqw = np.sqrt(mu/p)*np.array([-sinnu,(ecc+cosnu),0])
    
    # define rotation matrix to rotate perifocal frame to inertial
    # PI =  [cos(raan) * cos(arg_p) - sin(raan) * cos(inc) * sin(arg_p) -cos(raan) * sin(arg_p) - sin(raan) * cos(inc) * cos(arg_p) sin(raan) * sin(inc);
    #         sin(raan) * cos(arg_p) + cos(raan) * cos(inc) * sin(arg_p) -sin(raan) * sin(arg_p) + cos(raan) * cos(inc) * cos(arg_p) -cos(raan) * sin(inc);
    #         sin(inc) * sin(arg_p) sin(inc) * cos(arg_p) cos(inc);]
    
    PI = ROT3(-raan).dot(ROT1(-inc)).dot(ROT3(-arg_p))
    
    # rotate postion and velocity vectors to inertial frame
    R_ijk = np.dot(PI,R_pqw)
    V_ijk = np.dot(PI,V_pqw)

    return (R_ijk,V_ijk,R_pqw,V_pqw)





