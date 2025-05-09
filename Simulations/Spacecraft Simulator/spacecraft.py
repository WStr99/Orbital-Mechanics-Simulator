# Import libraries
import numpy as np
from scipy.integrate import solve_ivp

# Earth simulation class
class Earth:

    def __init__(self):
        # Constants
        self.mu = 398600.0  # Earth's gravitational parameter (km^3/s^2)
        self.Re = 6378.0    # Earth radius (km)

        # For drag simulation 
        self.Cd = 2.2
        self.A = 3.0        # m^2
        self.rho0 = 3.614e-13
        self.h0 = 700.0     # km
        self.H = 88.667     # km

        # J2 perturbation
        self.J2 = 1.08263e-3

    # Return gravitational constant
    def getMu(self):
        return self.mu

    # Return radius of the Earth
    def getR(self):
        return self.Re
    
class spacecraft(Earth):

    def __init__(self, name, color):

        # Inherit variables from Earth class
        super().__init__()
        
        # Name of spacecraft
        self.name = name
        self.color = color

        # Distance and Velocity class variables
        self.r = []
        self.v = []

        # Mass (kg)
        self.mass = 0

        # Max thrust magnitude in Newtons
        self.max_thrust = 90

        # COEs
        self.a = 0 # semi-major axis
        self.e = 0 # eccentricity
        self.i_deg = 0 # inclination
        self.Omega_deg = 0 # RAAN
        self.omega_deg = 0 # argument of periapsis
        self.nu_deg = 0 # true anomaly

    # Getter function
    def getCOEs(self):
        
        # Returns list of COEs
        COE = [self.a, 
               self.e,  
               self.i_deg, 
               self.Omega_deg,  
               self.omega_deg,
               self.nu_deg]
        return COE

    # Setter for distance
    def setInitialDistance(self, x, y, z):
        self.r = [x, y, z]

    # Setter for velocity
    def setInitialVelocity(self, x, y, z):
        self.v = [x, y, z]

    # Setter for mass (kg)
    def setMass(self, mass):
        self.mass = mass

    # Getter for distance
    def getInitialDistance(self):
        return self.r

    # Getter for velocity
    def getInitialVelocity(self):
        return self.v

    # Getter for mass
    def getMass(self):
        return self.mass
    
    # Returns state vector
    def getStateVector(self):
        return np.concatenate((self.r, self.v))
    
    def twoBodyEquations(self, y):

        # Distance and velocity vectors
        r = y[:3]
        v = y[3:]

        norm_r = np.linalg.norm(r)
        x, y_pos, z = r

        # Gravity
        a_gravity = -self.mu * r / norm_r**3

        # J2 perturbation
        z2 = z * z
        r2 = norm_r * norm_r
        factor = 1.5 * self.J2 * self.mu * (self.Re**2) / norm_r**5
        ax_J2 = factor * x * (5 * z2 / r2 - 1)
        ay_J2 = factor * y_pos * (5 * z2 / r2 - 1)
        az_J2 = factor * z * (5 * z2 / r2 - 3)
        a_J2 = np.array([ax_J2, ay_J2, az_J2])

        # Atmospheric drag
        alt = norm_r - self.Re
        rho = self.rho0 * np.exp(-(alt - self.h0) / self.H)
        v_mag = np.linalg.norm(v)
        a_drag = -0.5 * self.Cd * self.A / self.mass * rho * v_mag * v

        # Total acceleration
        a_total = a_gravity + a_J2 + a_drag

        return np.concatenate((v, a_total))

    def stepSimulation(self, t0, y0, dt):
        sol = solve_ivp(lambda t, y: self.twoBodyEquations(y), 
                        [t0, t0 + dt], y0, method='RK45', 
                        rtol=1e-9, atol=1e-9)
        return sol.y[:, -1]

    # Approximate ISS initial position and velocity
    def simulateISS(self, mu):

        # Set mass
        self.mass = 400000
        
        # Set distance
        self.r = np.array([6798.0, 0.0, 0.0])  # km #6798.0, 0.0, 0.0 (8798.0, 0.0, 0.0)
        v_mag = np.sqrt(mu / np.linalg.norm(self.r))  # ≈ 7.67 km/s
        inclination = np.radians(51.64)

        # Set velocity
        v_y = v_mag * np.cos(inclination)
        v_z = v_mag * np.sin(inclination)
        self.v = np.array([0.0, v_y, v_z]) #[0.0, v_y, v_z] ([5, 7, 1])

    # Computing orbital elements
    def ComputeCOE(self):
        
        # Orbital Elements
        r = np.linalg.norm(self.r)
        v = np.linalg.norm(self.v)
        h_vec = np.cross(self.r, self.v)
        h = np.linalg.norm(h_vec)
        e_vec = (np.cross(self.v, h_vec) / self.mu) - (self.r / r)
        energy = v**2 / 2 - self.mu / r
        K = np.array([0, 0, 1])
        n_vec = np.cross(K, h_vec)
        n = np.linalg.norm(n_vec)

        # Eccentricity
        e = np.linalg.norm(e_vec)
        
        # Semi-major axis
        a = -(self.mu) / (2 * energy)
        
        # Inclination
        i = np.arccos(h_vec[2] / h)
        
        # Right Ascension of the Ascending Node (RAAN)
        Omega = np.arccos(n_vec[0] / n) if n != 0 else 0
        if n_vec[1] < 0:
            Omega = 2*np.pi - Omega
        
        # Argument of periapsis
        omega = np.arccos(np.dot(n_vec, e_vec) / (n * e)) if n != 0 and e != 0 else 0
        if e_vec[2] < 0:
            omega = 2*np.pi - omega
        
        # True anomaly
        nu = np.arccos(np.dot(e_vec, self.r) / (e * r)) if e != 0 else 0
        if np.dot(self.r, self.v) < 0:
            nu = 2 * np.pi - nu

        # Update class variables
        self.a = a
        self.e = e

        # Convert radians to degrees
        self.i_deg = np.degrees(i)
        self.Omega_deg = np.degrees(Omega)
        self.omega_deg = np.degrees(omega)
        self.nu_deg = np.degrees(nu)

    def addThrust(self):

        # Convert inclination difference to radians
        #delta_i_rad = np.radians(iss_coe[2] - rocket_coe[2])

        # Compute scalar delta-v for inclination change
        #v_mag = np.linalg.norm(rocket_v)
        #delta_v_mag = 2 * v_mag * np.sin(delta_i_rad / 2)

        # Choose thrust direction orthogonal to current velocity to change inclination (simplified)
        thrust_dir = np.cross(self.r, self.v)
        thrust_dir = thrust_dir / np.linalg.norm(thrust_dir)

        # Compute acceleration from thrust
        thrust = (self.max_thrust / self.mass) * thrust_dir

        # Apply acceleration over a time step dt to get delta-v
        delta_t = 10.0  # seconds (or whatever your sim time step is)
        delta_v = thrust * delta_t

        # Update velocity
        v_new = self.v + delta_v

        return v_new
    
    def orbit(self, t, dt):

        # Get current state vector of spacecraft
        stateV = self.getStateVector()

        # Update velocity to include added thrust
        stateV[3:] = self.addThrust()

        # Simulate new position with two body equation
        stateV = self.stepSimulation(t, stateV, dt)

        # Update new position
        self.r = stateV[:3]
        self.v = stateV[3:]

        # Compute new COEs
        self.ComputeCOE()

        # Return new state vector (position and velocity)
        return stateV
    
