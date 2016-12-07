import numpy as np 
import utilities.attitude as attitude

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
        r = p/(1+ecc*np.cos(v));
        rs = p/(1+ecc*np.cos(nu_i));
    
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
    dcm_pqw2eci = np.dot(np.dot(attitude.ROT3(-raan),attitude.ROT1(-inc)),attitude.ROT3(-arg_p));

    orbit_plane = np.dot(dcm_pqw2eci,np.array([x,y,z]));

    x = orbit_plane[0,:]
    y = orbit_plane[1,:]
    z = orbit_plane[2,:]

    sat_pos = np.dot(dcm_pqw2eci,np.array([xs,ys,zs]));

    xs = sat_pos[0]
    ys = sat_pos[1]
    zs = sat_pos[2]

    return (x,y,z,xs,ys,zs)

def nu2anom(nu,ecc):
    """
    [E M] = ecc_anomaly(nu,ecc)

       Purpose:
           - Calculates the eccentric and mean anomaly given eccentricity and
           true anomaly

       Inputs:
           - nu - true anomaly in rad -2*pi < nu < 2*pi
           - ecc - eccentricity of orbit 0 < ecc < inf

       Outputs:
           - E - (elliptical/parabolic/hyperbolic) eccentric anomaly in rad
               0 < E < 2*pi
           - M - mean anomaly in rad 0 < M < 2*pi

       Dependencies:
           - none

       Author:
           - Shankar Kulumani 5 Dec 2016
                - Convert to python
           - Shankar Kulumani 15 Sept 2012
               - modified from USAFA code and notes from AAE532
               - only elliptical case will add other later
           - Shankar Kulumani 17 Sept 2012
               - added rev check to reduce angle btwn 0 and 2*pi

       References
           - AAE532 notes
           - Vallado 3rd Ed
    """

    small = 1e-9

    if ecc <= small: # circular
         E = nu
         M = nu
    elif small < ecc and ecc <= 1-small: # elliptical
        sine = ( np.sqrt( 1.0 -ecc*ecc ) * np.sin(nu) ) / ( 1.0 +ecc*np.cos(nu) )
        cose = ( ecc + np.cos(nu) ) / ( 1.0  + ecc*np.cos(nu) )
        
        E   = np.arctan2( sine,cose )
        M   = E - ecc*np.sin(E)
        
        E = attitude.normalize(E,0,2*np.pi)
        M = attitude.normalize(M,0,2*np.pi)

    elif np.absolute(ecc-1) <= small: # parabolic
        B = np.tan(nu/2)

        E = B
        M = B + 1.0/3*B**3

        # E = revcheck(E);
        # M = revcheck(M);
    elif ecc > 1+small: # hyperbolic
            sine = ( np.sqrt(ecc**2-1) * np.sin(nu) ) / ( 1.0  + ecc*np.cos(nu) )
            H = np.arcsinh( sine )
            E = H
            M = ecc*np.sinh(H) - H
            
            # E = revcheck(E);
            # M = revcheck(M);
    else:
        print("Eccentricity is out of bounds 0 < ecc < inf")
    
    return (E, M)

def tof_delta_t(p,ecc,mu,nu_0,delta_t):
    """
        Propogate a COE into the future
    """

    tol = 1e-9
    # calculate initial eccentric anomaly and mean anomaly
    E_0, M_0 = nu2anom(nu_0,ecc)

    # check eccentricity
    if np.absolute(ecc-1) < tol: # parabolic
        n = 2*np.sqrt(mu/p**3)
    else:
        a = p/(1-ecc**2)
        # calculate mean motion
        n = np.sqrt(mu/a**3)
    

    # calculate mean anomaly after delta t

    M_f = M_0 + n * delta_t
    k = np.floor(M_f/(2*np.pi))
    M_f = M_f-2*np.pi*k
    # calculate eccentric anomaly from mean anomaly (newton iteration)

    E_f,nu_f,count = kepler_eq_E(M_f,ecc)

    return (E_f, M_f, nu_f)

