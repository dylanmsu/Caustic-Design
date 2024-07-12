import matplotlib.pyplot as plt
from xml.etree import ElementTree as ET

# Read the pointset data from the file
filename = '/home/dylan/lena.dat'
with open(filename, 'r') as file:
    lines = file.readlines()

# Extract the x and y coordinates from each line
x_coords = []
y_coords = []
for line in lines:
    x, y = map(float, line.strip().split())
    x_coords.append(x)
    y_coords.append(y)

# Set the width and height of the SVG image (in inches)
width, height = 8, 6

# Create a new figure with the specified dimensions and white background
plt.figure(figsize=(width, height), facecolor='white')

# Plot the points
plt.scatter(x_coords, y_coords, s=5, color='black')
plt.axis('equal')
plt.axis('off')

plt.show()