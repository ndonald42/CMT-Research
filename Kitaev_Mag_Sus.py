import os
import numpy as np
import matplotlib.pyplot as plt

# Define filename function
def get_filename(j, h):
    base_path = r"C:\Users\Noah Donald\Documents\Junior Year\Research\Summer Research\Data\Total_Kitaev"
    return os.path.join(base_path, f"J_{j}", f"H_{h}", "OutC")

# Define magnetization function with output format option
def Mag(J, H, output_format='H'):
    file_path = get_filename(J, H)
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        # Find index of " Ritz values: "
        ritz_index = next(i for i, line in enumerate(lines) if " Ritz values: " in line)
        # Extract 24 lines after ritz_index + 2
        Sx, Sy, Sz = 0, 0, 0
        for i in range(24):
            line = lines[ritz_index + 2 + i].strip().split()
            # Extract numbers after dropping first 2 characters of each component
            Sx += float(line[0][2:])
            Sy += float(line[1][2:])
            Sz += float(line[2][2:])
        mag = np.sqrt(Sx**2 + Sy**2 + Sz**2)
        return [H, mag] if output_format == 'H' else [J, mag]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return [H, 0] if output_format == 'H' else [J, 0]

# Section 1: Magnetic susceptibility for fixed J
H_range = np.arange(0, 0.40, 0.01)  # 0 to 0.39 inclusive with step 0.01
deriv_kz05 = [
    [(Mag(0.5, H + 0.01, 'H')[1] - Mag(0.5, H, 'H')[1]) / 0.01 + (H - 1), 0]
    for H in H_range
]
deriv_kz1 = [
    [(Mag(1, H + 0.01, 'H')[1] - Mag(1, H, 'H')[1]) / 0.01 + (H - 1), 0]
    for H in H_range
]
deriv_kz3 = [
    [(Mag(3, H + 0.01, 'H')[1] - Mag(3, H, 'H')[1]) / 0.01 + (H - 1), 0]
    for H in H_range
]

# Plot Section 1
plt.figure(figsize=(8, 6))
H_deriv05 = [x[0] for x in deriv_kz05]
mag_deriv05 = [x[1] for x in deriv_kz05]
H_deriv1 = [x[0] for x in deriv_kz1]
mag_deriv1 = [x[1] for x in deriv_kz1]
H_deriv3 = [x[0] for x in deriv_kz3]
mag_deriv3 = [x[1] for x in deriv_kz3]

plt.plot(H_deriv05, mag_deriv05, label='Kz=0.5', linewidth=2)
plt.plot(H_deriv1, mag_deriv1, label='Kz=1', linewidth=2)
plt.plot(H_deriv3, mag_deriv3, label='Kz=3', linewidth=2)
plt.scatter(H_deriv05, mag_deriv05, s=20)
plt.scatter(H_deriv1, mag_deriv1, s=20)
plt.scatter(H_deriv3, mag_deriv3, s=20)
plt.legend(fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel('H - 1', fontsize=16)
plt.ylabel('\chi', fontsize=16)
plt.title('Magnetic Susceptibility vs H', fontsize=16)
plt.show()

# Section 2: Magnetization for fixed J
H_range = np.arange(0, 0.41, 0.01)  # 0 to 0.4 inclusive with step 0.01
mag_j0 = [Mag(0, H, 'H') for H in H_range]
mag_j05 = [Mag(0.5, H, 'H') for H in H_range]
mag_j1 = [Mag(1, H, 'H') for H in H_range]

# Plot Section 2
plt.figure(figsize=(8, 6))
H_j0 = [x[0] for x in mag_j0]
mag_vals_j0 = [x[1] for x in mag_j0]
H_j05 = [x[0] for x in mag_j05]
mag_vals_j05 = [x[1] for x in mag_j05]
H_j1 = [x[0] for x in mag_j1]
mag_vals_j1 = [x[1] for x in mag_j1]

plt.plot(H_j0, mag_vals_j0, label='J=0', linewidth=2)
plt.plot(H_j05, mag_vals_j05, label='J=0.5', linewidth=2)
plt.plot(H_j1, mag_vals_j1, label='J=1', linewidth=2)
plt.scatter(H_j0, mag_vals_j0, s=20)
plt.scatter(H_j05, mag_vals_j05, s=20)
plt.scatter(H_j1, mag_vals_j1, s=20)
plt.legend(fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel('H', fontsize=16)
plt.ylabel('Magnetization', fontsize=16)
plt.title('Magnetization vs H', fontsize=16)
plt.show()

# Section 3: Magnetization and magnetic susceptibility for fixed H
J_range = np.arange(0, 1.01, 0.05)  # 0 to 1 inclusive with step 0.05
mag_h0 = [Mag(J, 0, 'J') for J in J_range]
mag_h005 = [Mag(J, 0.05, 'J') for J in J_range]
mag_h01 = [Mag(J, 0.1, 'J') for J in J_range]
deriv_h0 = [
    [(Mag(J + 0.05, 0, 'J')[1] - Mag(J, 0, 'J')[1]) / 0.05 + (J - 1), 0]
    for J in J_range
]
deriv_h005 = [
    [(Mag(J + 0.05, 0.05, 'J')[1] - Mag(J, 0.05, 'J')[1]) / 0.05 + (J - 1), 0]
    for J in J_range
]
deriv_h01 = [
    [(Mag(J + 0.05, 0.1, 'J')[1] - Mag(J, 0.1, 'J')[1]) / 0.05 + (J - 1), 0]
    for J in J_range
]

# Plot Section 3: Magnetization
plt.figure(figsize=(8, 6))
J_h0 = [x[0] for x in mag_h0]
mag_vals_h0 = [x[1] for x in mag_h0]
J_h005 = [x[0] for x in mag_h005]
mag_vals_h005 = [x[1] for x in mag_h005]
J_h01 = [x[0] for x in mag_h01]
mag_vals_h01 = [x[1] for x in mag_h01]

plt.plot(J_h0, mag_vals_h0, label='H=0', linewidth=2)
plt.plot(J_h005, mag_vals_h005, label='H=0.05', linewidth=2)
plt.plot(J_h01, mag_vals_h01, label='H=0.1', linewidth=2)
plt.scatter(J_h0, mag_vals_h0, s=20)
plt.scatter(J_h005, mag_vals_h005, s=20)
plt.scatter(J_h01, mag_vals_h01, s=20)
plt.legend(fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel('J', fontsize=16)
plt.ylabel('Magnetization', fontsize=16)
plt.title('Magnetization vs J', fontsize=16)
plt.show()

# Plot Section 3: Magnetic Susceptibility
plt.figure(figsize=(8, 6))
J_deriv0 = [x[0] for x in deriv_h0]
mag_deriv0 = [x[1] for x in deriv_h0]
J_deriv005 = [x[0] for x in deriv_h005]
mag_deriv005 = [x[1] for x in deriv_h005]
J_deriv01 = [x[0] for x in deriv_h01]
mag_deriv01 = [x[1] for x in deriv_h01]

plt.plot(J_deriv0, mag_deriv0, label='H=0', linewidth=2)
plt.plot(J_deriv005, mag_deriv005, label='H=0.05', linewidth=2)
plt.plot(J_deriv01, mag_deriv01, label='H=0.1', linewidth=2)
plt.scatter(J_deriv0, mag_deriv0, s=20)
plt.scatter(J_deriv005, mag_deriv005, s=20)
plt.scatter(J_deriv01, mag_deriv01, s=20)
plt.legend(fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel('J - 1', fontsize=16)
plt.ylabel('\chi', fontsize=16)
plt.title('Magnetic Susceptibility vs J', fontsize=16)
plt.show()