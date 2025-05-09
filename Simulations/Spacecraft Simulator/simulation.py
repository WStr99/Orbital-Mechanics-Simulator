# Perform basic orbital manuvers to change COEs
# Add attitude adjustment w Euler here or in another program
# Simulate HCW here or in another program

# Import libraries
import numpy as np
from scipy.integrate import solve_ivp

# Import classes
from visualization import visualization
from spacecraft import spacecraft

# Main simulation class
class simulation:

    def __init__(self):
        # Constants
        self.mu = 398600.0  # Earth's gravitational parameter (km^3/s^2)
        self.Re = 6378.0    # Earth radius (km)

        # Time variables
        self.t = 0 # Starting time of simulation
        self.dt = 100 # Simulates every 100 seconds

        # For drag simulation 
        self.Cd = 2.2
        self.A = 3.0        # m^2
        self.m = 1000.0     # kg
        self.rho0 = 3.614e-13
        self.h0 = 700.0     # km
        self.H = 88.667     # km

        # J2 perturbation
        self.J2 = 1.08263e-3

        # Thrust magnitude in Newtons
        self.thrust_magnitude = 90


    def simulate(self):

        vis = visualization()

        # Rocket setup
        rocket = spacecraft('Rocket', 'red')
        rocket.setInitialDistance(7000.0, 0.0, 0.0)
        rocket.setInitialVelocity(0, 7, -1)
        rocket.setMass(12545) # Approx max of SpaceX Dragon

        # ISS setup
        iss = spacecraft('ISS', 'green')
        iss.simulateISS(self.mu)

        while True:

            rocketState = rocket.orbit(self.t, self.dt)
            issState = iss.orbit(self.t, self.dt)

            vis.visualize(rocket.name, rocket.color, rocketState[0], rocketState[1], rocketState[2], rocket.getCOEs())
            vis.visualize(iss.name, iss.color,  issState[0], issState[1], issState[2])

# Main method
if __name__ == "__main__":
    sim = simulation()
    sim.simulate()

# Handle while loop in spacecraft class
# Handle visualization in spacecraft class
# Apply thrust from 
# Create PID controller

