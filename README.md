# Orbital-Mechanics-Simulator
Orbital Mechanics Simulator in Python simulating a two body orbit with Newton's laws of motion.

This project simulates and visualizes a 3D satellite orbit around Earth using Python. It demonstrates classical orbital mechanics concepts including elliptical trajectories, orbital inclination, and the computation of key orbital elements.

## Features

3D simulation of Earth-centered elliptical orbits using `solve_ivp` from `scipy`
Real-time animation of satellite motion along a numerically integrated orbit
Visualization of:
  - Perigee and apogee (fastest and slowest orbital points)
  - Earth as one orbital focus, with computed second focus and ellipse center
  - Equatorial plane and orbit inclination
Calculation of classical orbital elements (COEs):
  - Semi-major axis (a)
  - Eccentricity (e)
  - Inclination (i)
  - RAAN (Ω)
  - Argument of periapsis (ω)
  - True anomaly (ν)

# Customization

You can tune initial position (r0) and velocity (v0) vectors inside orbit_3d_animated.py to control:
    - Orbit inclination (e.g. set to ~150°)
    - Orbit shape (eccentricity)
    - Initial speed and orientation

## Screenshots / Media

Orbital-diagram.png
    - detailed diagram of COEs
Orbital-Simulation-1/2.MOV
    - Video of orbital simulation

## ⚙️ How to Run

1. Clone the repo:
   ```bash
   git clone https://github.com/yourusername/orbital-simulation.git
   cd orbital-simulation


## Install Requirments

pip install -r requirements.txt

## Run the simulation

python two_body_simulator.py
