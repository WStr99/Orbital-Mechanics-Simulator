import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

# Constants
G = 6.67430e-11  # m^3/kg/s^2
m_sun = 1.989e30
m_earth = 5.972e24
m_mars = 6.4171e23
m_rocket = 5.972e24/6
masses = [m_sun, m_earth, m_mars, m_rocket]  # rocket has no mass

# Orbital parameters
r_earth_mag = 1.496e11  # m
r_mars_mag = 2.279e11   # m
v_earth_mag = 29780     # m/s
v_mars_mag = 24077      # m/s

# Initial positions
r_sun = np.array([0.0, 0.0, 0.0])
r_earth = np.array([r_earth_mag, 0.0, 0.0])

# Rotate Mars position 44.4 degrees ahead of Earth
angle_deg = 44.4
angle_rad = np.radians(angle_deg)
r_mars = r_mars_mag * np.array([np.cos(angle_rad), np.sin(angle_rad), 0.0])

# Initial velocities
v_sun = np.array([0.0, 0.0, 0.0])
v_earth = np.array([0.0, v_earth_mag, 0.0])
v_mars = v_mars_mag * np.array([-np.sin(angle_rad), np.cos(angle_rad), 0.0])

# Approximate Lambert solution for rocket (hardcoded)
v_rocket = np.array([2200.0, 35580.0, 0.0])
r_rocket = np.array([r_earth_mag + 30000000.0, 0.0, 0.0])

# State vector
initial_state = np.concatenate([
    r_sun, v_sun,
    r_earth, v_earth,
    r_mars, v_mars,
    r_rocket, v_rocket
])

def compute_accelerations(positions):
    n = len(positions)
    accelerations = np.zeros_like(positions)
    for i in range(n):
        for j in range(len(masses)):
            if i != j:
                r_ij = positions[j] - positions[i]
                dist = np.linalg.norm(r_ij)
                accelerations[i] += G * masses[j] * r_ij / dist**3
    return accelerations

def deriv(t, state):
    positions = np.array([
        state[0:3], state[6:9], state[12:15], state[18:21]
    ])
    velocities = np.array([
        state[3:6], state[9:12], state[15:18], state[21:24]
    ])
    accelerations = compute_accelerations(positions)
    return np.concatenate([np.concatenate([velocities[i], accelerations[i]]) for i in range(4)])

# Integrate
t_span = (0, 3.154e7)  # 1 year
t_eval = np.linspace(*t_span, 2000)
sol = solve_ivp(deriv, t_span, initial_state, t_eval=t_eval, rtol=1e-9)

# Extract positions
sun_pos = sol.y[0:3].T
earth_pos = sol.y[6:9].T
mars_pos = sol.y[12:15].T
rocket_pos = sol.y[18:21].T

# =================================================== Plotting ===================================================

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-3e11, 3e11)
ax.set_ylim(-3e11, 3e11)
ax.set_zlim(-1e10, 1e10)

sun_plot = ax.scatter([], [], [], color='yellow', s=300, label='Sun')
earth_plot = ax.scatter([], [], [], color='blue', s=100, label='Earth')
mars_plot = ax.scatter([], [], [], color='red', s=100, label='Mars')
rocket_plot = ax.scatter([], [], [], color='green', s=20, label='Rocket')

earth_trail, = ax.plot([], [], [], 'b--', linewidth=0.5)
mars_trail, = ax.plot([], [], [], 'r--', linewidth=0.5)
rocket_trail, = ax.plot([], [], [], 'g--', linewidth=0.5)

def update(frame):
    sun_plot._offsets3d = ([sun_pos[frame][0]], [sun_pos[frame][1]], [sun_pos[frame][2]])
    earth_plot._offsets3d = ([earth_pos[frame][0]], [earth_pos[frame][1]], [earth_pos[frame][2]])
    mars_plot._offsets3d = ([mars_pos[frame][0]], [mars_pos[frame][1]], [mars_pos[frame][2]])
    rocket_plot._offsets3d = ([rocket_pos[frame][0]], [rocket_pos[frame][1]], [rocket_pos[frame][2]])

    earth_trail.set_data(earth_pos[:frame+1, 0], earth_pos[:frame+1, 1])
    earth_trail.set_3d_properties(earth_pos[:frame+1, 2])
    mars_trail.set_data(mars_pos[:frame+1, 0], mars_pos[:frame+1, 1])
    mars_trail.set_3d_properties(mars_pos[:frame+1, 2])
    rocket_trail.set_data(rocket_pos[:frame+1, 0], rocket_pos[:frame+1, 1])
    rocket_trail.set_3d_properties(rocket_pos[:frame+1, 2])

    return sun_plot, earth_plot, mars_plot, rocket_plot, earth_trail, mars_trail, rocket_trail

ani = animation.FuncAnimation(
    fig,
    update,
    frames=range(0, len(t_eval), 10),
    interval=10,
    blit=False
)

ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Rocket Transfer: Earth to Mars (Lambert Approx)')
ax.legend()
plt.show()