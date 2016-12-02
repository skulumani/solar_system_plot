import numpy as np 
from keplerian_orbit import *
from utilities import *

def element_approx(JD_curr, planet_flag):
    """
    Given current JD this function will output the current orbital elements for 
    all of the planets
    """

    # Find the JD centuries past J2000 epoch
    JD_J2000 = 2451545.0
    T = (JD_curr - JD_J2000)/36525

    # load the elements and element rates for the planets
    (a0,adot,e0,edot,inc0,incdot,meanL0,meanLdot,lonperi0,lonperidot,raan0,raandot,b,c,f,s) = planet_elements(JD_curr,planet_flag)

    # compute the current elements at this JD
    a = a0 + adot*T
    ecc = e0 + edot*T
    inc = inc0 + incdot*T
    L = meanL0 + meanLdot*T
    lonperi = lonperi0 + lonperidot*T
    raan = raan0 + raandot*T

    # compute argp and M/v to complete the element set
    argp = lonperi - raan
    M = L - lonperi + b*T**2 + c*np.cos(f*T) + s*np.sin(f*T)

    # solve kepler's equation to compute E and v
    E, nu, count = kepler_eq_E(M_in,ecc_in)
    # package into a vector and output


def planet_elements(JD,planet_flag):
    """
    This outputs the orbital element information data for the selected planet. 
    Also does some logic checking based on the desired time period
    """

    JD_1800AD = 2378497.000000 
    JD_2050AD = 2469808.000000

    JD_3000BC = 625674.000000
    JD_3000AD = 2816788.000000

    if JD < JD_1800AD or JD > JD_2050AD: # have to use the coarse model

        
        print("Date is outside 1800-2050. Need to implement the coarse model.")
    else: # date is between 1800-2050AD and use the fine model
        b = 0.0
        c = 0.0
        f = 0.0
        s = 0.0

        if planet_flag == 1: # mercury
            a0 = 0.38709927
            adot = 0.00000037

            e0 = 0.20563593
            edot = 0.00001906

            inc0 = 7.00497902
            incdot = -0.00594749

            meanL0 = 252.25032350
            meanLdot = 149472.67411175

            lonperi0 = 77.45779628
            lonperidot = 0.16047689

            raan0 = 48.33076593
            raandot = -0.12534081

        elif planet_flag == 2: # venus
            a0 = 0.72333566
            adot = 0.00000390

            e0 = 0.00677672 
            edot = -0.00004107

            inc0 = 3.39467605
            incdot = -0.00078890

            meanL0 = 181.97909950
            meanLdot = 58517.81538729

            lonperi0 = 131.60246718
            lonperidot = 0.00268329

            raan0 = 76.67984255
            raandot = -0.27769418
        elif planet_flag == 3: # earth moon barycenter
            a0 = 1.00000261
            adot = 0.00000562

            e0 = 0.01671123 
            edot = -0.00004392

            inc0 = -0.00001531
            incdot = -0.01294668

            meanL0 = 100.46457166
            meanLdot = 35999.37244981

            lonperi0 = 102.93768193
            lonperidot = 0.32327364

            raan0 = 0.0
            raandot = 0.0
        elif planet_flag == 4: # mars
            a0 = 1.52371034
            adot = 0.00001847

            e0 = 0.09339410 
            edot = 0.00007882

            inc0 = 1.84969142
            incdot = -0.00813131

            meanL0 = -4.55343205
            meanLdot = 19140.30268499

            lonperi0 = -23.94362959
            lonperidot = 0.44441088

            raan0 = 49.55953891
            raandot = -0.29257343
        elif planet_flag == 5: # jupiter
            a0 = 5.20288700
            adot = -0.00011607

            e0 = 0.04838624 
            edot = -0.00013253

            inc0 = 1.30439695 
            incdot = -0.00183714

            meanL0 = 34.39644051
            meanLdot = 3034.74612775

            lonperi0 = 14.72847983
            lonperidot = 0.21252668

            raan0 = 100.47390909
            raandot = 0.20469106
        elif planet_flag == 6: # saturn
            a0 = 9.53667594
            adot = -0.00125060

            e0 = 0.05386179 
            edot = -0.00050991

            inc0 = 2.48599187 
            incdot = 0.00193609

            meanL0 = 49.95424423
            meanLdot = 1222.49362201

            lonperi0 = 92.59887831
            lonperidot = -0.41897216

            raan0 = 113.66242448
            raandot = -0.28867794
        elif planet_flag == 7: # uranus
            a0 = 19.18916464
            adot = -0.00196176

            e0 = 0.04725744 
            edot = -0.00004397

            inc0 = 0.77263783 
            incdot = -0.00242939

            meanL0 = 313.23810451
            meanLdot = 428.48202785

            lonperi0 = 170.95427630
            lonperidot = 0.40805281

            raan0 = 74.01692503
            raandot = 0.04240589
        elif planet_flag == 8: # neptune
            a0 = 30.06992276
            adot = 0.00026291

            e0 = 0.00859048
            edot = 0.00005105

            inc0 = 1.77004347  
            incdot = 0.00035372

            meanL0 = -55.12002969
            meanLdot = 218.45945325

            lonperi0 = 44.96476227 
            lonperidot = -0.32241464

            raan0 = 131.78422574
            raandot = -0.00508664
        elif planet_flag == 9: # pluto
            a0 = 39.48211675
            adot = -0.00031596

            e0 = 0.24882730 
            edot = 0.00005170

            inc0 = 17.14001206  
            incdot = 0.00004818

            meanL0 = 238.92903833
            meanLdot = 145.20780515

            lonperi0 = 224.06891629 
            lonperidot = -0.04062942

            raan0 = 110.30393684
            raandot = -0.01183482
        else:
            print("Incorrect planet flag should be between 1 and 9")

        
        return (a0,adot,e0,edot,inc0,incdot,meanL0,meanLdot,lonperi0,lonperidot,raan0,raandot,b,c,f,s)

if __name__ == "__main__":
    JD_curr = 2457724.500000

    (a0,adot,e0,edot,inc0,incdot,meanL0,meanLdot,lonperi0,lonperidot,raan0,raandot,b,c,f,s) = planet_elements(JD_curr,3)
    # coe = element_approx(JD_curr,3)

# Table 2a.

# Keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000,
# valid for the time-interval 3000 BC -- 3000 AD.  NOTE: the computation of M for Jupiter through
# Pluto *must* be augmented by the additional terms given in Table 2b (below).

#                a              e               I                L            long.peri.      long.node.
#            AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
# ------------------------------------------------------------------------------------------------------
# Mercury   0.38709843      0.20563661      7.00559432      252.25166724     77.45771895     48.33961819
#           0.00000000      0.00002123     -0.00590158   149472.67486623      0.15940013     -0.12214182
# Venus     0.72332102      0.00676399      3.39777545      181.97970850    131.76755713     76.67261496
#          -0.00000026     -0.00005107      0.00043494    58517.81560260      0.05679648     -0.27274174
# EM Bary   1.00000018      0.01673163     -0.00054346      100.46691572    102.93005885     -5.11260389
#          -0.00000003     -0.00003661     -0.01337178    35999.37306329      0.31795260     -0.24123856
# Mars      1.52371243      0.09336511      1.85181869       -4.56813164    -23.91744784     49.71320984
#           0.00000097      0.00009149     -0.00724757    19140.29934243      0.45223625     -0.26852431
# Jupiter   5.20248019      0.04853590      1.29861416       34.33479152     14.27495244    100.29282654
#          -0.00002864      0.00018026     -0.00322699     3034.90371757      0.18199196      0.13024619
# Saturn    9.54149883      0.05550825      2.49424102       50.07571329     92.86136063    113.63998702
#          -0.00003065     -0.00032044      0.00451969     1222.11494724      0.54179478     -0.25015002
# Uranus   19.18797948      0.04685740      0.77298127      314.20276625    172.43404441     73.96250215
#          -0.00020455     -0.00001550     -0.00180155      428.49512595      0.09266985      0.05739699
# Neptune  30.06952752      0.00895439      1.77005520      304.22289287     46.68158724    131.78635853
#           0.00006447      0.00000818      0.00022400      218.46515314      0.01009938     -0.00606302
# Pluto    39.48686035      0.24885238     17.14104260      238.96535011    224.09702598    110.30167986
#           0.00449751      0.00006016      0.00000501      145.18042903     -0.00968827     -0.00809981
# ------------------------------------------------------------------------------------------------------



# Table 2b.

# Additional terms which must be added to the computation of M
# for Jupiter through Pluto, 3000 BC to 3000 AD, as described
# in the related document.

#                 b             c             s            f 
# ---------------------------------------------------------------
# Jupiter   -0.00012452    0.06064060   -0.35635438   38.35125000
# Saturn     0.00025899   -0.13434469    0.87320147   38.35125000
# Uranus     0.00058331   -0.97731848    0.17689245    7.67025000
# Neptune   -0.00041348    0.68346318   -0.10162547    7.67025000
# Pluto     -0.01262724
# ---------------------------------------------------------------