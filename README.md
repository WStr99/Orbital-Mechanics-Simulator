# Orbital-Mechanics-Simulator
This project is an Object-Oriented modular Python suite for simulating and visualizing spacecraft orbits, orbital maneuvers, and interplanetary transfers. It includes support for gravity modeling, atmospheric drag, J2 perturbations, real-time COE tracking, and multiple simulation modes like Earth orbit and Earth-to-Mars transfer.

## Features

Earth-Centric Orbit Simulation (Spacecraft simulation)
- Models gravitational effects, atmospheric drag, and J2 perturbation.
- Simulates n number of spacecraft and ISS orbits with thrust input.
- Real-time 3D visualization using matplotlib.
- Thrust modeled in orthagonal direction to change inclination
- Computes real-time COEs
- Object Oriented Design

<img width="698" alt="Screenshot 2025-05-08 at 11 19 14 PM" src="https://github.com/user-attachments/assets/5928168a-b4c1-4675-889d-b326b46ca0be" />


Hohmann Transfer Simulator (test_scripts/orbital_transfer.py)
- Models an interplanetary transfer from Earth to Mars.
- Multi-body integration of Sun, Earth, Mars, and Rocket.
- Animates orbital paths using approximate Lambert solution.

<img width="700" alt="Screenshot 2025-05-08 at 11 34 49 PM" src="https://github.com/user-attachments/assets/a5d4223e-df14-4ae6-8d59-df4e6a8d5579" />


Basic Two-Body Simulator (test_scripts/two_body_simulatior.py)
- Accurate RK45 integration of elliptical orbits.
- COE computation and animation of orbital motion.

<img width="667" alt="Screenshot 2025-05-08 at 11 37 48 PM" src="https://github.com/user-attachments/assets/880d8096-2388-4982-b231-9dc2773a7531" />

## Known Bugs
  - in orbital_transfer.py, program crashes when saatelite position == earth's position (divide by 0 error), as a result the satelite must start outside of earth's orbit and will not be fully captured by Mars' orbit.

## Future Improvements 
- in two_body_simulator.py the ability to input thrust in real-time will be added in order to vizualize how this affects orbit
- in orbital_transfer.py, known bug will be fixed and as a result Lambert's problem & Hohmann Transfer will be able to be solved fully within the program rather than hard coded
- in orbital_transfer.py, other planets in the solar system will be added
- in in orbital_transfer.py, the ability to add new bodies to the solar system will be added

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
