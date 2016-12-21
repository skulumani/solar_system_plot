## Solar System Plot

This package will plot the positions of all planets (and Pluto) in the solar system.
In addition, this will also read in the orbital elements of any other objects, 
such as asteroids, or spacecraft, and plot those relative to the planets.

Finally, the state history of these asteroids is output to a text file for use
in other simulations.

To use just simply run `plot_planets.py` to generate a plot of the solar system
planets. 
You can use `plot_asteroid.py` to plot the planets with an asteroid. 
Currently, it is setup for three asteroids

1. Bennu
2. Itokawa
3. 2008 EV5

## Testing

You need to install `pytest` to run the self tests. 

Then navigate to the root directory and run `pytest` 