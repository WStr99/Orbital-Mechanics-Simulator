# Import libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

class visualization:
    
    def __init__(self):

        # Create figure and axis
        self.fig = plt.figure(figsize=(10, 8), facecolor='black')
        self.ax = self.fig.add_subplot(111, projection='3d', facecolor='black')

        # Earth at origin
        self.ax.scatter(0, 0, 0, color='blue', s=1000)
        self.ax.scatter(0, 0, 0, color='blue', s=100, label='Earth') # For plotting

        # Labels and formatting
        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlabel("X (km)", color='white')
        self.ax.set_ylabel("Y (km)", color='white')
        self.ax.set_zlabel("Z (km)", color='white')
        self.ax.set_title("Real-Time Orbit Simulation", color='white')
        self.ax.tick_params(colors='white')
        self.ax.legend()

        # Fixed Earth-centered view
        view_limit = 10000  # km
        self.ax.set_xlim(-view_limit, view_limit)
        self.ax.set_ylim(-view_limit, view_limit)
        self.ax.set_zlim(-view_limit, view_limit)

        # Storage for multiple spacecraft: {id: (marker, trail, [xs, ys, zs])}
        self.spacecraft_data = {}

        # COE text placeholder (top-left corner in axes coordinates)
        self.coe_text = self.ax.text2D(0.02, 0.98, "", transform=self.ax.transAxes, color='white', fontsize=10, verticalalignment='top', family='monospace')

        # Enable interactive mode
        plt.ion()
        plt.show()

    # Register a new spacecraft with a unique ID and color 
    def add_spacecraft(self, name, color='red'):
        # Add trail and satellite marker
        marker, = self.ax.plot([], [], [], '^', color=color, markersize=7)
        trail, = self.ax.plot([], [], [], linestyle=':', color='white', linewidth=1)

        # Placeholder label (will be updated later)
        label = self.ax.text(0, 0, 0, name, color='white', fontsize=7)

        # Store data for this spacecraft
        self.spacecraft_data[name] = {
            'marker': marker,
            'trail': trail,
            'label': label,
            'xs': [],
            'ys': [],
            'zs': []
    }

    # Visualize spacecraft orbits
    def visualize(self, spacecraft_id, color, x, y, z, coe=None):
        
        # Adds spacecraft if previously unknown
        if spacecraft_id not in self.spacecraft_data:
            self.add_spacecraft(spacecraft_id, color)

        # Stores orbit data (Can limit amount of data stored)
        data = self.spacecraft_data[spacecraft_id]
        data['xs'].append(x)
        data['ys'].append(y)
        data['zs'].append(z)

        # Update marker and trail
        data['marker'].set_data([x], [y])
        data['marker'].set_3d_properties([z])
        data['trail'].set_data(data['xs'][-20:], data['ys'][-20:])
        data['trail'].set_3d_properties(data['zs'][-20:])

        # Update label position
        data['label'].set_position((x, y))
        data['label'].set_3d_properties(z + 100)

        # Update COE display (top-left text box)
        if coe:
            elements = ['a', 'e', 'i', 'Ω', 'ω', 'ν']
            coe_lines = [f"{name} = {val:.3f}" for name, val in zip(elements, coe)]
            title = f"{spacecraft_id.upper()} COEs:"
            self.coe_text.set_text(title + "\n" + "\n".join(coe_lines))

        # Show animation
        self.fig.canvas.draw()
        plt.pause(0.001)