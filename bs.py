import numpy as np
import matplotlib.pyplot as plt

def burning_ship_fractal(xmin, xmax, ymin, ymax, width, height, max_iterations):
    x = np.linspace(xmin, xmax, width).reshape((1, width))
    y = np.linspace(ymin, ymax, height).reshape((height, 1))
    c = x + y * 1j
    z = np.zeros_like(c)

    for i in range(max_iterations):
        z = (np.abs(z.real) + 1j * np.abs(z.imag))**2 + c

    burning_ship = (np.abs(z) < 2).astype(int)

    return burning_ship

# Set the parameters for the fractal
xmin, xmax = -2, 2
ymin, ymax = -2, 2
width, height = 1800, 1800
max_iterations = 100

# Generate the Burning Ship fractal
burning_ship = burning_ship_fractal(xmin, xmax, ymin, ymax, width, height, max_iterations)

# Plot the fractal
plt.imshow(burning_ship, extent=(xmin, xmax, ymin, ymax), cmap='hot', interpolation='bilinear')
plt.colorbar()
plt.title('Burning Ship Fractal')
plt.show()
