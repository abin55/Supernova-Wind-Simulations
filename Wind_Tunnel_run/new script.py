import numpy as np
import pyPLUTO as pluto
import matplotlib.pyplot as plt
import matplotlib.animation as animation

w_dir = "/home/domin/PLUTO/Test_Problems/HD/Wind_Tunnel/"

# Physical units
UNIT_DENSITY  = 3.5e-12
UNIT_LENGTH   = 1.496e13
UNIT_VELOCITY = 1.0e5        
UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY**2

with open(w_dir + "dbl.out") as f:
    lines = [l for l in f.readlines() if l.strip() and not l.startswith('#')]
valid_ns = [int(l.split()[0]) for l in lines]
n_snaps  = len(valid_ns)
print(f"Valid snapshots: {valid_ns}")


rho_cube  = []
vmag_cube = []
prs_cube  = []
for ns in valid_ns:
    D = pluto.Load(ns, path=w_dir)
    rho_cube.append(D.d_vars['rho'].T  * UNIT_DENSITY)
    vmag_cube.append(np.sqrt(D.d_vars['vx1'].T**2 + D.d_vars['vx2'].T**2) * UNIT_VELOCITY)
    prs_cube.append(D.d_vars['prs'].T  * UNIT_PRESSURE)


ns0 = valid_ns[0]
D0  = pluto.Load(ns0, path=w_dir)
r     = D0.x1  
theta = D0.x2  
print(f"Grid: r shape={r.shape}, theta shape={theta.shape}")

def centers_to_edges(x):
    edges = np.zeros(len(x) + 1)
    edges[1:-1] = 0.5 * (x[1:] + x[:-1])
    edges[0]  = x[0]  - (edges[1] - x[0])
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return edges

r_edges     = centers_to_edges(r)
theta_edges = centers_to_edges(theta)

R_edges, THETA_edges = np.meshgrid(r_edges, theta_edges, indexing='ij')
X_edges = R_edges * np.cos(THETA_edges)   # AU
Y_edges = R_edges * np.sin(THETA_edges)   # AU


def to_plot_shape(field):
    return field.

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

im0 = axes[0].pcolormesh(X_edges, Y_edges, to_plot_shape(rho_cube[0]),
                          cmap='inferno', shading='flat')
im1 = axes[1].pcolormesh(X_edges, Y_edges, to_plot_shape(vmag_cube[0]),
                          cmap='plasma', shading='flat')
im2 = axes[2].pcolormesh(X_edges, Y_edges, to_plot_shape(prs_cube[0]),
                          cmap='viridis', shading='flat')

plt.colorbar(im0, ax=axes[0], label=f'ρ × {UNIT_DENSITY:.2e} g/cm³')
plt.colorbar(im1, ax=axes[1], label=f'|v| × {UNIT_VELOCITY:.2e} cm/s')
plt.colorbar(im2, ax=axes[2], label=f'p × {UNIT_PRESSURE:.2e} dyne/cm²')

axes[0].set_title('Density')
axes[1].set_title('Velocity')
axes[2].set_title('Pressure')
for ax in axes:
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_aspect('equal') 

suptitle = fig.suptitle('')

def update(i):

    im0.set_array(to_plot_shape(rho_cube[i]).ravel())
    im1.set_array(to_plot_shape(vmag_cube[i]).ravel())
    im2.set_array(to_plot_shape(prs_cube[i]).ravel())
    suptitle.set_text(f'Snapshot {i:04d}')
    return im0, im1, im2, suptitle

ani = animation.FuncAnimation(fig, update, frames=n_snaps, interval=500, blit=False)
ani.save('wind_simulation_cartesian.mp4', writer='ffmpeg', fps=2, dpi=150)
plt.close()
print("Saved wind_simulation_cartesian.mp4!")
