# Orbital-Mechanics-Simulator
Orbital Mechanics Simulator in Python simulating both a two-body orbit around Earth and an interplanetary transfer from Earth to Mars using Newton’s laws and classical orbital mechanics.
This project includes:
  - A 3D satellite orbit simulator demonstrating elliptical motion and orbital elements around Earth.
  - An orbital transfer simulation illustrating an interplanetary maneuver using Hohmann transfer principles and an approximate Lambert solution.


## Features

Earth-Centered Two-Body Orbit (two_body_simulator.py)

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

Interplanetary Transfer: Earth to Mars (orbit_transfer.py)

Models Earth, Mars, Sun, and rocket interactions using Newtonian gravity
  - Implements an interplanetary transfer based on:
  - Hohmann Transfer concept
  - Approximate Lambert’s problem solution for the rocket’s injection vector
Real-time 3D animation of:
  - Earth and Mars orbiting the Sun
  - Rocket transferring from Earth to Mars
  - Ideal for understanding interplanetary trajectories and mission planning concepts

# Customization

In two_body_simulator.py, you can tune initial position (r0) and velocity (v0) vectors inside orbit_3d_animated.py to control:
  - Orbit inclination (e.g. set to ~150°)
  - Orbit shape (eccentricity)
  - Initial speed and orientation

In orbit_transfer.py, you can experiment with:
  - Launch vector v_rocket
  - Initial spacing between Earth and Mars
  - Integration span (transfer time window)

## Screenshots / Media

Orbital-diagram.png
    - detailed diagram of COEs

Orbital-Simulation-1/2.MOV
  - Video of orbital simulation
    
Orbital-Tranfer.mov
  - Video of orbital tarnsfer from Earth to Mars

<img width="582" alt="two_body" src="https://github.com/user-attachments/assets/b2d19a8d-d2b7-432c-acb4-b2f1bcd5c9b7" />

<img width="582" alt="orbital_transfer" src="https://github.com/user-attachments/assets/da625694-8af5-4544-a01c-d61e0cbf0a0c" />

## Known Bugs
  - in orbital_transfer.py, program crashes when saatelite position == earth's position (divide by 0 error), as a result the satelite must start outside of earth's orbit and will not be fully captured by Mars' orbit.

## Future Improvements 
- in two_body_simulator.py the ability to input thrust in real-time will be added in order to vizualize how this affects orbit
- in orbital_transfer.py, known bug will be fixed and Lambert's problem will be able to be solved within the program rather than hard coded

## ⚙️ How to Run

1. Clone the repo:
   ```bash
   git clone https://github.com/yourusername/orbital-simulation.git
   cd Orbital-Simulatior

## Install Requirments

pip install -r requirements.txt

## Run the simulation

python two_body_simulator.py
python orbit_transfer.py
