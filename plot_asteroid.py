# script to plot the orbits of several asteroids
import numpy as np

# load the COEs for the asteroids
# 2008 EV5
coe_ev5 = (0.9582899238313918*(1 - 0.08348599378460778**2),
    0.08348599378460778,np.deg2rad(7.436787362690259),
    np.deg2rad(93.39122898787916),np.deg2rad(234.8245876826614),
    np.deg2rad(90))
# Itokawa

# Bennu

# plot all the planets

# plot each asteroid

# convert COE to RV

# transform RV to Sun Earth 3BP and nondimensionalize


import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import keplerian_orbit.keplerian_orbit as kpl

# test Kepler Equation solver
M_in = np.array(0.5)
ecc_in = np.array(0)
E_out, nu_out, count_out = kpl.kepler_eq_E(M_in,ecc_in)

print("E: %5.4f M: %5.4f count: %5.4f" %( E_out,nu_out,count_out))

# define orbital elements
p = 6578
ecc = 0.2
inc = np.deg2rad(10)
raan = np.deg2rad(12)
argp = np.deg2rad(15)
nu = np.deg2rad(165)
# call conic orbit plotter
(x,y,z,xs,ys,zs) = kpl.conic_orbit(p,ecc, inc, raan, argp, nu, nu)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z,'b')
ax.plot([xs],[ys],[zs],'ro')
plt.axis('equal')
plt.show()

# test True anomaly converter
print("Convert \\nu to E and M")