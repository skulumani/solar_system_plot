# script to plot the orbits of several asteroids
import numpy as np 
from datetime import datetime
from utilities.time import date2jd
from utilities.attitude import normalize
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from keplerian_orbit.keplerian_orbit import conic_orbit, kepler_eq_E, tof_delta_t
from orbital_elements.planet_coe import planet_coe

from orbital_elements.asteroid_coe import asteroid_coe

from keplerian_orbit.coe import coe2rv


km2au = 1/149597870.700
au2km = 1/km2au
mu = 1.32712440018e11 # km^3/sec^2

def write_to_file():
    """
        Write 10 periods of the J2000 state to a text file
    """
    # loop over the asteriods and compute the state
    asteroid_names = ('EV5','Itokawa','Bennu')
    # plot each asteroid
    for ast_flag in range(3):
        (p,ecc,inc,raan,argp,nu) = asteroid_coe(JD_curr,ast_flag)
        p = p*au2km

        # compute period of orbit
        period = 10*2*np.pi*np.sqrt((p/(1-ecc**2))**3/mu) # period in seconds

        with open(asteroid_names[ast_flag] + ".txt", "w") as text_file:

            print("Asteroid: {} state wrt Sol barycenter ( t(sec) x(km) y(km) z(km) vx(km/sec) vy(km/sec) vz(km/sec)".format(asteroid_names[ast_flag]), file=text_file)
            # loop over nu and compute COE2RV and print to text file
            time_span = np.arange(0,period,86400)
            for idx,t_curr in enumerate(time_span):

                # propogate epoch to current t
                nu_curr = tof_delta_t(p,ecc,mu,nu,t_curr)[2]
                # convert COE to RV
                r_ijk, v_ijk, r_pqw, v_pqw = coe2rv(p,ecc,inc,raan,argp,nu_curr,mu)
                # print to text file
                print("%16.16f %16.16f %16.16f %16.16f %16.16f %16.16f %16.16f" % (t_curr, r_ijk[0], r_ijk[1], r_ijk[2], v_ijk[0], v_ijk[1], v_ijk[2]), file=text_file)

    return 0

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
for planet_flag in range(9):
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

ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])
ax.set_title('Solar System JD %10.3f' % JD_curr)
plt.axis('off')
ax.view_init(azim=0, elev=90)
plt.show()


