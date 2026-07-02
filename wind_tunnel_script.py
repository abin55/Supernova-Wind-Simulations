import numpy as np
import pyPLUTO as pluto
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re

UNIT_DENSITY  = 1.67e-24
UNIT_LENGTH   = 1.496e13
UNIT_VELOCITY = 1.0e5
UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY**2

base    = "/home/domin/PLUTO/Test_Problems/HD/"
n_steps = 51

# Field functions
def rho_field(D):
    return D.d_vars['rho'].T * UNIT_DENSITY

def vmag_field(D):
    vx1 = D.d_vars['vx1'].T
    vx2 = D.d_vars['vx2'].T
    return np.sqrt(vx1**2 + vx2**2) * UNIT_VELOCITY

# Load
def load_field_cube(folder, field_func):
    w_dir = f"{base}{folder}/"
    D0 = pluto.Load(0, path=w_dir)
    NY, NX = D0.d_vars['rho'].T.shape
    cube = np.zeros((n_steps, NY, NX))
    for i in range(n_steps):
        D = pluto.Load(i, path=w_dir)
        cube[i] = field_func(D)
    return cube, D0.x1 * UNIT_LENGTH, D0.x2 * UNIT_LENGTH

# Comparison animation
def make_comparison(folder_a, label_a, folder_b, label_b,
                    field_func, field_key, cmap, outname):

    cube_a, x1, x2 = load_field_cube(folder_a, field_func)
    cube_b, _,  _  = load_field_cube(folder_b, field_func)

    # Each panel uses its OWN scale
    vmin_a, vmax_a = cube_a.min(), cube_a.max()
    vmin_b, vmax_b = cube_b.min(), cube_b.max()

    if field_key == 'vmag':
        cb_label = f"|v| × v₀ = {UNIT_VELOCITY:.2e} cm/s"
    elif field_key == 'rho':
        cb_label = f"ρ × ρ₀ = {UNIT_DENSITY:.2e} g/cm³"
    elif field_key == 'prs':
        cb_label = f"p × p₀ = {UNIT_PRESSURE:.2e} dyne/cm²"

    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    im0 = axes[0].pcolormesh(x1, x2, cube_a[0], cmap=cmap, vmin=vmin_a, vmax=vmax_a)
    axes[0].set_title(f"{label_a}\nmax = {vmax_a/UNIT_VELOCITY:.1f} × v₀  (v₀ = {UNIT_VELOCITY:.2e} cm/s)")
    plt.colorbar(im0, ax=axes[0], label=cb_label)

    im1 = axes[1].pcolormesh(x1, x2, cube_b[0], cmap=cmap, vmin=vmin_b, vmax=vmax_b)
    axes[1].set_title(f"{label_b}\nmax = {vmax_b/UNIT_VELOCITY:.1f} × v₀  (v₀ = {UNIT_VELOCITY:.2e} cm/s)")
    plt.colorbar(im1, ax=axes[1], label=cb_label)

    for ax in axes:
        ax.set_xlabel(f'X × L₀ = {UNIT_LENGTH:.3e} cm')
        ax.set_ylabel(f'Y × L₀ = {UNIT_LENGTH:.3e} cm')
        ax.set_aspect('equal')

    suptitle = fig.suptitle('')

    def update(i):
        im0.set_array(cube_a[i].ravel())
        im1.set_array(cube_b[i].ravel())
        suptitle.set_text(f'step {i:04d}')
        return im0, im1, suptitle

    ani = animation.FuncAnimation(fig, update, frames=n_steps,
                                  interval=100, blit=False)
    ani.save(outname, writer='ffmpeg', fps=3, dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")

# Run 
make_comparison("Wind_Tunnel_VX1_4.0", "VX1=4.0",
                "Wind_Tunnel_VX1_50", "VX1=50.0",
                vmag_field, 'vmag', "plasma", "compare_vx1_4_vs_50_vmag.mp4")

make_comparison("Wind_Tunnel_VX1_4.0", "VX1=4.0",
                "Wind_Tunnel_VX1_50", "VX1=50.0",
                rho_field, 'rho', "inferno", "compare_vx1_4_vs_50_rho.mp4")

print("Done!")
