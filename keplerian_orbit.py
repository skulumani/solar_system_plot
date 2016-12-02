import numpy as np 

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

if __name__ == "__main__":
    # test Kepler Equation solver
    M_in = np.array(0.5)
    ecc_in = np.array(0)
    E_out, nu_out, count_out = kepler_eq_E(M_in,ecc_in)

    print("E: %5.4f M: %5.4f count: %5.4f" %( E_out,nu_out,count_out))