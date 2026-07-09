import numpy as np
import pyPLUTO as pluto
import matplotlib.pyplot as plt
import matplotlib.animation as animation

w_dir = "/home/domin/PLUTO/Test_Problems/HD/Wind_Tunnel/"

# Physical units
UNIT_DENSITY  = 3.5e-12
UNIT_LENGTH   = 1.496e13
UNIT_VELOCITY = 1.0e5        # update this once you check Mstar in pluto.ini
UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY**2

# Read valid snapshots from dbl.out
with open(w_dir + "dbl.out") as f:
    lines = [l for l in f.readlines() if l.strip() and not l.startswith('#')]

valid_ns = [int(l.split()[0]) for l in lines]
n_snaps  = len(valid_ns)
print(f"Valid snapshots: {valid_ns}")

# Load all snapshots
rho_cube  = []
vmag_cube = []
prs_cube  = []

for ns in valid_ns:
    D = pluto.Load(ns, path=w_dir)
    rho_cube.append(D.d_vars['rho'].T  * UNIT_DENSITY)
    vmag_cube.append(np.sqrt(D.d_vars['vx1'].T**2 + D.d_vars['vx2'].T**2) * UNIT_VELOCITY)
    prs_cube.append(D.d_vars['prs'].T  * UNIT_PRESSURE)


x1 = D.x1   # r in AU
x2 = D.x2   # θ in radians
print(f"Grid: x1 shape={x1.shape}, x2 shape={x2.shape}")

fig, axes = plt.subplots(1, 3, figsize=(16, 4))

# Initial plots
im0 = axes[0].pcolormesh(x1, x2, rho_cube[0],  cmap='inferno', shading='auto')
im1 = axes[1].pcolormesh(x1, x2, vmag_cube[0], cmap='plasma')
im2 = axes[2].pcolormesh(x1, x2, prs_cube[0],  cmap='viridis')

plt.colorbar(im0, ax=axes[0], label=f'ρ × {UNIT_DENSITY:.2e} g/cm³')
plt.colorbar(im1, ax=axes[1], label=f'|v| × {UNIT_VELOCITY:.2e} cm/s')
plt.colorbar(im2, ax=axes[2], label=f'p × {UNIT_PRESSURE:.2e} dyne/cm²')

axes[0].set_title('Density')
axes[1].set_title('Velocity')
axes[2].set_title('Pressure')

for ax in axes:
    ax.set_xlabel('r (AU)')
    ax.set_ylabel('θ (rad)')

suptitle = fig.suptitle('')

def update(i):
    im0.set_array(rho_cube[i].ravel())
    im1.set_array(vmag_cube[i].ravel())
    im2.set_array(prs_cube[i].ravel())
    suptitle.set_text(f'Snapshot {i:04d}')
    return im0, im1, im2, suptitle

ani = animation.FuncAnimation(fig, update, frames=n_snaps, interval=500, blit=False)
ani.save('wind_simulation_MS_stage.mp4', writer='ffmpeg', fps=2, dpi=150)
plt.close()
print("Saved wind_simulation_MS_stage.mp4!")
