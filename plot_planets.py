import numpy as np 
from datetime import datetime
from utilities.time import date2jd
from utilities.attitude import normalize
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from keplerian_orbit.keplerian_orbit import conic_orbit, kepler_eq_E  
from orbital_elements.planet_coe import planet_coe

def plot_planets(JD):
    # function to draw all of the planets

    planet_names = ('Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # loop over all the planets and compute the COE
    for planet_flag in range(4):
        # calculate the conic orbit for each planet 
        p,ecc,inc,raan,argp,nu = planet_coe(JD,planet_flag)

        (x,y,z,xs,ys,zs) = conic_orbit(p,ecc, inc, raan, argp, nu, nu)
        # plot the planet to figure
        ax.plot(x,y,z,'b')
        ax.plot([xs],[ys],[zs],'ro')
        ax.text(xs,ys,zs,planet_names[planet_flag])

    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    ax.set_zlim([-5,5])
    plt.axis('equal')    
    plt.show()


if __name__ == "__main__":

    today_date = datetime.today()
    yr = today_date.year
    mon = today_date.month
    day = today_date.day
    hr = today_date.hour
    minute = today_date.minute
    sec = today_date.second

    JD_curr,MJD = date2jd(yr,mon,day,hr,minute,sec)
    planet_flag = 3
    p,ecc,inc,raan,argp,nu = planet_coe(JD_curr,planet_flag)

    print("p: %16.16f au" % p)
    print("a: %16.16f au" % (p/(1-ecc**2)))
    print("ecc: %16.16f " % ecc)
    print("inc: %16.16f deg" % np.rad2deg(inc))
    print("raan: %16.16f deg" % np.rad2deg(raan))
    print("argp: %16.16f deg" % np.rad2deg(argp))
    print("nu: %16.16f deg" % np.rad2deg(nu))

    plot_planets(JD_curr)


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