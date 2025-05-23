import os
import numpy as np
import matplotlib.pyplot as plt

# Define the first filename function. Extracts data from C++ output CSV files.
def filename(j, h):
    base_path = r"C:\Users\Noah Donald\Documents\Junior Year\Research\Summer Research\Data\site_num_24"
    return os.path.join(base_path, f"J_{j}", f"H_{h}", "TopoEE.dat")

# Define the second filename function. Extracts data from larger C++ output CSV files.
def filename2(j, h):
    base_path = r"C:\Users\Noah Donald\Documents\Junior Year\Research\Summer Research\Data\Total_Kitaev"
    return os.path.join(base_path, f"J_{j}", f"H_{h}", "TopoEE.dat")

# Define the topological entanglement entropy function
def Mag(J, H):
    file_path = filename(J, H)
    try:
        # Read the file
        with open(file_path, 'r') as f:
            lines = f.readlines()
        # Extract the first number from the third line
        mag = float(lines[2].strip().split()[0])
        return [H, J, mag]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return [H, J, 0]  # Return 0 for mag if file reading fails

# Define the Mag2 function
def Mag2(J, H):
    file_path = filename2(J, H)  # Note: Using filename2 for Mag2
    try:
        # Read the file (assuming text file with at least 3 lines)
        with open(file_path, 'r') as f:
            lines = f.readlines()
        # Extract the first number from the third line
        mag = float(lines[2].strip().split()[0])
        return [H, mag]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return [H, 0]  # Return 0 for mag if file reading fails

# Generate List2 for contour plot
H_range = np.arange(0, 0.25, 0.01)  # 0 to 0.24 inclusive with step 0.01
J_range = np.arange(0, 4.01, 0.01)  # 0 to 4 inclusive with step 0.01
List = [[Mag(J, H) for J in J_range] for H in H_range]

# Convert List2 to arrays for plotting
H_vals, J_vals, mag_vals = [], [], []
for row in List:
    for entry in row:
        H_vals.append(entry[0])
        J_vals.append(entry[1])
        mag_vals.append(entry[2])

# Create topological entanglement entropy contour plot
plt.figure(figsize=(8, 6))
contour = plt.tricontourf(H_vals, J_vals, mag_vals, levels=50, cmap='inferno')
plt.colorbar(contour)
plt.tick_params(labelsize=16)
plt.xlabel('H', fontsize=16)
plt.ylabel('J', fontsize=16)
plt.title('Topological Entanglement Entropy', fontsize=16)
plt.show()

# Generate Lists for line plots that cut across the countour plot with different K_z values
H_range_line = np.arange(0, 0.41, 0.01)  # 0 to 0.4 inclusive with step 0.01
List11 = [Mag2(0.5, H) for H in H_range_line]
List21 = [Mag2(1, H) for H in H_range_line]
List31 = [Mag2(3, H) for H in H_range_line]

# Extract H and topological entanglement entropy values for plotting
H11, mag11 = [x[0] for x in List11], [x[1] for x in List11]
H21, mag21 = [x[0] for x in List21], [x[1] for x in List21]
H31, mag31 = [x[0] for x in List31], [x[1] for x in List31]

# Create line plot with markers
plt.figure(figsize=(8, 6))
plt.plot(H11, mag11, label='Kz=0.5', linewidth=2)
plt.plot(H21, mag21, label='Kz=1', linewidth=2)
plt.plot(H31, mag31, label='Kz=3', linewidth=2)
plt.scatter(H11, mag11, s=20)
plt.scatter(H21, mag21, s=20)
plt.scatter(H31, mag31, s=20)
plt.legend(fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel('H', fontsize=16)
plt.ylabel('Topological Entanglement Entropy', fontsize=16)
plt.title('Cuts across topological entanglement entropy countour plot', fontsize=16)
plt.show()