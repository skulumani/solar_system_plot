# script to plot the orbits of several asteroids
import numpy as np 
from datetime import datetime
from utilities.time import date2jd
from utilities.attitude import normalize
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from keplerian_orbit.keplerian_orbit import conic_orbit, kepler_eq_E  
from orbital_elements.planet_coe import planet_coe

from orbital_elements.asteroid_coe import asteroid_coe

from keplerian_orbit.coe import coe2rv

km2au = 1/149597870.700

# plot all the planets
today_date = datetime.today()
yr = today_date.year
mon = today_date.month
day = today_date.day
hr = today_date.hour
minute = today_date.minute
sec = today_date.second

JD_curr,MJD = date2jd(yr,mon,day,hr,minute,sec)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot the planets
planet_names = ('Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto')
# loop over all the planets and compute the COE
for planet_flag in range(4):
    # calculate the conic orbit for each planet 
    p,ecc,inc,raan,argp,nu = planet_coe(JD_curr,planet_flag)

    (x,y,z,xs,ys,zs) = conic_orbit(p,ecc, inc, raan, argp, nu, nu)
    # plot the planet to figure
    ax.plot(x,y,z,'b')
    ax.plot([xs],[ys],[zs],'ro')
    ax.text(xs,ys,zs,planet_names[planet_flag])
   

asteroid_names = ('EV5','Itokawa','Bennu')
# plot each asteroid
for ast_flag in range(3):
    (p,ecc,inc,raan,argp,nu) = asteroid_coe(JD_curr,ast_flag)

    (x,y,z,xs,ys,zs) = conic_orbit(p,ecc, inc, raan, argp, nu, nu)
    ax.plot(x,y,z,'g')
    ax.plot([xs],[ys],[zs],'ro')
    ax.text(xs,ys,zs,asteroid_names[ast_flag])

    # convert COE to RV and print to screen
    # make sure the units of the COEs are consistent
    p = 1/km2au * p # convert to km
    mu = 1.32712440018e11 # km^3/sec^2

    r_ijk, v_ijk, r_pqw, v_pqw = coe2rv(p,ecc,inc,raan,argp,nu,mu)

    # print to screen/text file
    print("%s state wrt Sol barycenter" % asteroid_names[ast_flag])
    print("    Pos (km): %16.16f %16.16f %16.16f" % (r_ijk[0], r_ijk[1], r_ijk[2]))
    print("    Vel (km/sec): %16.16f %16.16f %16.16f" % (v_ijk[0], v_ijk[1], v_ijk[2]))
# print the state of the Earth
p,ecc,inc,raan,argp,nu = planet_coe(JD_curr,3)
r_ijk, v_ijk, r_pqw, v_pqw = coe2rv(p*1/km2au,ecc,inc,raan,argp,nu,mu)

print("")
print("Earth State wrt sol barycenter")
print("    Pos (km): %16.16f %16.16f %16.16f" % (r_ijk[0], r_ijk[1], r_ijk[2]))
print("    Vel (km/sec): %16.16f %16.16f %16.16f" % (v_ijk[0], v_ijk[1], v_ijk[2]))

ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])
ax.set_title('Solar System JD %10.3f' % JD_curr)
plt.axis('off')
ax.view_init(azim=0, elev=90)
plt.show()
# transform RV to Sun Earth 3BP and nondimensionalize
