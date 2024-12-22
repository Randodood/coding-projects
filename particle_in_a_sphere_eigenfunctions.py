import matplotlib.pyplot as plt
import numpy as np
from scipy import special, optimize

# User inputs
try:
    n = int(input('Enter n (0 < n): '))
    l = int(input('Enter l (0 <= l < n): '))
    m = int(input('Enter m (-l <= m <= l): '))
    resolution = int(input('Enter plot resolution: '))
except:
    print('Non-integer value')

# Error conditions
if (n <= 0 or l < 0 or l >= n or m < -l or m > l or resolution <= 1):
    raise ValueError('Invalid quantum numbers')

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(8,7.5))
ax = fig.add_subplot(projection='3d')
ax.set_axis_off()

# Meshes for equatorial and azimuthal angles
t = np.linspace(0, np.pi, resolution)
p = np.linspace(0, 2*np.pi, resolution)
P, T = np.meshgrid(p, t)

# Quantum physics convention for spherical harmonics
if m >= 0:
    sph_lm = special.sph_harm(m,l,P,T).real
    modifier = r"$\Re$"
else:
    sph_lm = special.sph_harm(m,l,P,T).imag
    modifier = r"$\Im$"

# Negative radius correction
sph_lm_list = list(sph_lm)
for i in range (len(sph_lm_list)):
    for j in range (len(sph_lm_list[0])):
        if sph_lm_list[i][j] < 0:
                T[i][j] = T[i][j] + np.pi
            
# Code found online for solving for roots of spherical Bessel function
def spherical_bessel_zeros():
    if l != 0:
        x_grid = np.linspace(l, l+ 2*n*(np.pi*(np.log(l)+1)), resolution)
    else:
        x_grid = np.linspace(l, l+ 2*n*np.pi, resolution)
    y_grid = special.spherical_jn(l, x_grid)
    
    diffs = np.sign(y_grid)[1:] - np.sign(y_grid)[:-1]
    index_0s = np.where(diffs)[0][:n]
    x_0s = x_grid[index_0s]
    
    def fn(x):
        return special.spherical_jn(l, x)
    return [optimize.root(fn, x0).x[0] for x0 in x_0s]

zeros = spherical_bessel_zeros()
opacity = [(1/(i+1)) for i in range(n)]
scale = np.array(zeros)/zeros[n-1]

# Convsersion from spherical coordinates to cartesian coordinates
X = sph_lm*np.sin(T)*np.cos(P)
Y = sph_lm*np.sin(T)*np.sin(P)
Z = sph_lm*np.cos(T)

# Iteration over nodal surfaces of n values
for j in range(n):
    ax.plot_surface(scale[j]*X,scale[j]*Y,scale[j]*Z, alpha = opacity[j], cmap = plt.cm.YlGnBu_r)

# Plot finalization
axis_limit = 0.7*sph_lm.max()
plt.title(r"$\{(r,\theta,\phi)\in[0,a]\times[0,2\pi]\times[0,\pi]|$"+modifier+r"$[\Psi_{nlm}(r, \theta, \phi)] = 0\}$"+
          '\n'+r"$n="+str(n)+", l="+str(l)+", m="+str(m)+"$", loc="center", fontsize=18)
ax.set_zlim(-0.75*axis_limit,0.75*axis_limit)
ax.set_xlim(-axis_limit,axis_limit)
ax.set_ylim(-axis_limit,axis_limit)

plt.show()