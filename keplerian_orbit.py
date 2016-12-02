import numpy as np 
from utilities import *
import pdb
import matplotlib.pyplot as plt 

from mpl_toolkits.mplot3d import Axes3D

def kepler_eq_E(M_in,ecc_in):
    """
    (E,nu,count) = kepler_eq_E(M,ecc)
    Purpose:
       - This function solves Kepler's equation for eccentric anomaly
       given a mean anomaly using a newton-rapson method.
           - Will work for elliptical/parabolic/hyperbolic orbits

    Inputs:
       - M - mean anomaly in rad -2*pi < M < 2*pi
       - ecc - eccentricity 0 < ecc < inf

    Outputs:
       - E - eccentric anomaly in rad 0 < E < 2*pi
       - nu - true anomaly in rad 0 < nu < 2*pi

    Dependencies:
       - none

    Author:
       - Shankar Kulumani 15 Sept 2012
           - rewritten from code from USAFA
           - solve for elliptical orbits add others later
       - Shankar Kulumani 29 Sept 2012
           - added parabolic/hyperbolic functionality
       - Shankar Kulumani 7 Dec 2014
          - added loop for vector inputs
       - Shankar Kulumani 2 Dec 2016
          - converted to python and removed the vector inputs

    References
       - USAFA Astro 321 LSN 24-25
       - Vallado 3rd Ed pg 72
    """
    tol = 1e-6
    max_iter = 50

    M = M_in
    ecc = ecc_in
    # eccentricity check
    """
        HYPERBOLIC ORBIT
    """
    if ecc-1.0  > tol: # eccentricity logic
        # initial guess
        if ecc < 1.6: # initial guess logic
            if M < 0.0 and (M > -np.pi or M > np.pi):
                E_0 = M - ecc
            else:
                E_0 = M + ecc
            
        else:
            if ecc < 3.6 and np.absolute(M) > np.pi:
                E_0 = M - np.sign(M)*ecc;
            else:
                E_0 = M/(ecc-1.0 );
        
        # netwon's method iteration to find hyperbolic anomaly
        count = 1
        E_1 = E_0 + ( (M - ecc*np.sinh(E_0)+ E_0) / (ecc*np.cosh(E_0) - 1.0 ) )
        while ((np.absolute(E_1-E_0) > tol ) and ( count <= max_iter )):
            E_0 = E_1
            E_1 = E_0 + ( (M - ecc*np.sinh(E_0)+ E_0) / (ecc*np.cosh(E_0) - 1.0 ) )
            count = count + 1
        
        E = E_0
        # find true anomaly
        sinv = -( np.sqrt( ecc*ecc-1.0  ) * np.sinh(E_1) ) / ( 1.0  - ecc*np.cosh(E_1) )
        cosv = ( np.cosh(E_1) - ecc ) / ( 1.0  - ecc*np.cosh(E_1) )
        nu  = np.arctan2( sinv,cosv );
    else:
        """
            PARABOLIC
        """
        if np.absolute(ecc-1.0) < tol: # parabolic logic
            count= 1
            
            S = 0.5  * (np.pi/2 - np.arctan( 1.5 * M ) )
            W = np.arctan( np.tan( S )**(1.0 /3.0 ) )
            
            E = 2.0 *1.0/np.tan(2.0 *W)
            
            nu = 2.0  * np.arctan(E)
        else:
            """
                ELLIPTICAl
            """
            if  ecc > tol:   # elliptical logic
                
                # determine intial guess for iteration
                if M > -np.pi and (M < 0 or M > np.pi):
                    E_0 = M - ecc
                else:
                    E_0 = M + ecc
                                    
                # newton's method iteration to find eccentric anomaly
                count= 1
                E_1 = E_0 + ( M - E_0 + ecc*np.sin(E_0) ) / ( 1.0  - ecc*np.cos(E_0) )
                while (( np.absolute(E_1-E_0) > tol ) and ( count <= max_iter )):
                    count = count + 1
                    E_0 = E_1
                    E_1 = E_0 + ( M - E_0 + ecc*np.sin(E_0) ) / ( 1.0  - ecc*np.cos(E_0) )
                
                E = E_0
                
                # find true anomaly
                sinv = ( np.sqrt( 1.0 -ecc*ecc ) * np.sin(E_1) ) / ( 1.0 -ecc*np.cos(E_1) )
                cosv = ( np.cos(E_1)-ecc ) / ( 1.0  - ecc*np.cos(E_1) )
                nu  = np.arctan2( sinv,cosv )
            else:
                """
                    CIRCULAR
                """
                # -------------------- circular -------------------
                count= 0
                nu = M
                E = M
                
            
        

    return (E,nu,count)

def conic_orbit(p,ecc, inc, raan, arg_p, nu_i, nu_f):
    """Plot conic orbit
        
        Purpose: 
           - Uses the polar conic equation to plot a conic orbit
        
        [x y z xs ys zs ] = conic_orbit(p,ecc,inc,raan,arg_p,nu_i,nu_f)
        
        Inputs: 
           - p - semi-major axis (km)
           - ecc - eccentricity
           - raan - right acsension of the ascending node (rad) 0 < raan <
           2*pi
           - inc - inclination (rad) 0 < inc < pi
           - arg_p - argument of periapsis (rad) 0 < arg_p < 2*pi
           - nu_i - initial true anomaly (rad) 0 < nu < 2*pi
           - nu_f - final true anomaly (rad) 0 < nu < 2*pi
        
        Outputs: 
           - none
        
        Dependencies: 
           - ROT1,ROT2,ROT3 - principle axis rotation matrices
        
        Author: 
           - Shankar Kulumani 1 Dec 2012
               - list revisions
           - Shankar Kulumani 5 Dec 2014
               - added outputs for orbit gui functions
        
        References
           - AAE532 
    """

    tol = 1e-9
    step = 1000
    
    # v = true anomaly
    if nu_f > nu_i:
        v = np.linspace(nu_i,nu_f,step)
    else:
        v = np.linspace(nu_i,nu_f+2*np.pi,step)
    
    if ecc - 1 > tol: # hyperbolic
        turn_angle = np.acos(-1.0/ecc)
        v = np.linespace(-turn_angle,turn_angle,step);

        if nu_i > pi:
            nu_i = nu_i-2*np.pi

        r = p/(1+ecc*np.cos(v))
        rs = p/(1+ecc*np.cos(nu_i))

    elif np.absolute(ecc-1) < tol: #parabolic
        v = np.linspace(-np.pi,np.pi,step);
        if nu_i > np.pi:
            nu_i = nu_i-2*np.pi
        
        r = p/2*(1+np.tan(v/2)**2);
        rs = p/2*(1+np.tan(nu_i/2)**2);
    else:
        # conic equation for elliptical orbit
        r = p/(1+ecc*cos(v));
        rs = p/(1+ecc*cos(nu_i));
    
    x = r*np.cos(v)
    y = r*np.sin(v)
    z = np.zeros_like(x)

    xs = rs*np.cos(nu_i)
    ys = rs*np.sin(nu_i)
    zs = 0
    # rotate orbit plane to correct orientation

     # M_rot = [cos(raan) * cos(arg_p) - sin(raan) * cos(inc) * sin(arg_p) -cos(raan) * sin(arg_p) - sin(raan) * cos(inc) * cos(arg_p) sin(raan) * sin(inc);
     #         sin(raan) * cos(arg_p) + cos(raan) * cos(inc) * sin(arg_p) -sin(raan) * sin(arg_p) + cos(raan) * cos(inc) * cos(arg_p) -cos(raan) * sin(inc);
     #         sin(inc) * sin(arg_p) sin(inc) * cos(arg_p) cos(inc);];
    dcm_pqw2eci = np.dot(np.dot(ROT3(-raan),ROT1(-inc)),ROT3(-arg_p));

    orbit_plane = np.dot(dcm_pqw2eci,np.array([x,y,z]));

    x = orbit_plane[1,:]
    y = orbit_plane[2,:]
    z = orbit_plane[3,:]

    sat_pos = np.dot(dcm_pqw2eci,np.array([xs,ys,zs]));

    xs = sat_pos(1)
    ys = sat_pos(2)
    zs = sat_pos(3)


if __name__ == "__main__":
    # test Kepler Equation solver
    M_in = np.array(0.5)
    ecc_in = np.array(0)
    E_out, nu_out, count_out = kepler_eq_E(M_in,ecc_in)

    print("E: %5.4f M: %5.4f count: %5.4f" %( E_out,nu_out,count_out))

    # define orbital elements
    p = 6578
    ecc = 0.2
    inc = np.deg2rad(10)
    raan = np.deg2rad(12)
    argp = np.deg2rad(15)
    nu = np.deg2rad(165)
    # call conic orbit plotter
    (x,y,z,xs,ys,zs) = conic_orbit(p,ecc, inc, raan, argp, nu, nu)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x,y,z,'b')
    ax.plot([xs],[ys],[zs],'ro')

    plt.show()
