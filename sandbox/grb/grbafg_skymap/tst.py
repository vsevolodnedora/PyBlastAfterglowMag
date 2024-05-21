import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_3d_histogram_with_marginals(hist, x_array, y_array):
    """
    Plots a 3D histogram and its marginal histograms.

    Parameters:
    - hist: 2D histogram data.
    - x_array: Values along the X axis.
    - y_array: Values along the Y axis.
    """
    fig = plt.figure(figsize=(10, 8))

    # Main 3D histogram
    ax = fig.add_axes([0.0, 0.0, 0.7, 0.7], projection='3d')
    X, Y = np.meshgrid(x_array, y_array)
    ax.plot_surface(X, Y, hist.T, cmap='viridis')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Count')
    ax.view_init(elev=10, azim=-70, roll=0)

    # Top subplot for the X-axis marginal
    ax2 = fig.add_axes([0.1, 0.7, 0.6, 0.2])
    ax2.bar(x_array, np.sum(hist, axis=1), width=x_array[1]-x_array[0])
    ax2.set(xlim=ax.get_xlim(), ylabel='Count')

    # Right subplot for the Y-axis marginal
    ax3 = fig.add_axes([0.7, 0.1, 0.2, 0.6])
    ax3.barh(y_array, np.sum(hist, axis=0), height=y_array[1]-y_array[0])
    ax3.set(ylim=ax.get_ylim(), xlabel='Count')

    plt.tight_layout()
    plt.show()

# Example usage:
hist_data = np.random.rand(10, 10)  # Replace with your 2D histogram data
x = np.linspace(0, 9, 10)  # Replace with your x-axis data
y = np.linspace(0, 9, 10)  # Replace with your y-axis data
# plot_3d_histogram_with_marginals(hist_data, x, y)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create random heatmap data
def random_heatmap_data():
    return np.random.rand(10, 10)

fig = plt.figure()
G = gridspec.GridSpec(2, 2)  # 2x2 grid

# Function to plot 3D line and heatmaps
def plot_3D_and_heatmaps(ax):
    ax.plot3D([1, 2, 3], [1, 2, 3], [1, 2, 3])

    x = np.linspace(0, 3, 11)  # Create an array of 11 points from 0 to 3
    z = np.linspace(0, 3, 11)
    X, Z = np.meshgrid(x, z)

    y_positions = [0.5, 1.5, 2.5]  # Y-positions to stack the heatmaps
    for y in y_positions:
        Y = np.ones_like(X) * y
        ax.pcolormesh(X, Y, Z, shading='auto', cmap='hot', alpha=0.7, antialiased=True)

# Top-left subplot
ax1 = fig.add_subplot(G[0, 0], projection='3d')
plot_3D_and_heatmaps(ax1)

# Top-right subplot
ax2 = fig.add_subplot(G[0, 1], projection='3d')
plot_3D_and_heatmaps(ax2)

# Bottom-left subplot
ax3 = fig.add_subplot(G[1, 0], projection='3d')
plot_3D_and_heatmaps(ax3)

plt.tight_layout()
plt.show()
