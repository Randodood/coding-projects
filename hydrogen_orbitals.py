import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
from scipy import special

# User inputs
try:
    n = int(input('Enter n (0 < n): '))
    l = int(input('Enter l (0 <= l < n): '))
    m = int(input('Enter m (-l <= m <= l): '))
    resolution = int(input('Enter plot resolution: '))
    nb_points = int(input('Enter number of points: '))
except:
    print('Non-integer value')

# Appropriate error conditions
if (n <= 0 or l < 0 or l >= n or m < -l or m > l or resolution <= 1 or nb_points <= 0):
    raise ValueError('Invalid numbers')

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(8,7.5))
ax = fig.add_subplot(projection='3d')
ax.set_axis_off()

# Physical variables and constants
a_0 = 1
r_limit = 30*a_0

dr = r_limit/resolution
dphi = 2*np.pi/resolution

r = np.linspace(0, r_limit, resolution)
theta = np.linspace(0, 2*np.pi, resolution)
phi  = np.linspace(0, np.pi, resolution)

# Discretization and normalization of r and phi probablitity density functions
def prob_r(n_, l_, r_):
    vec_1 = r_**(2*l_+2)*np.exp(-2*r_/(n_*a_0))
    vec_2 = np.polyval(special.genlaguerre(n_-l_-1,2*l_+1),2*r_/(n_*a_0))**2
    P = vec_1*vec_2
    return P/np.sum(P)

def prob_phi(m_, l_, phi_):
    P=(special.lpmv(m_, l_, np.cos(phi_))**2)*np.sin(phi_)
    return P/np.sum(P)

# Generation of points based on probability functions
r_gen = np.random.choice(r, size=nb_points, p=prob_r(n,l,r))
theta_gen = np.random.choice(theta, size=nb_points)
phi_gen = np.random.choice(phi, size=nb_points, p=prob_phi(m, l, phi))

# Conversion from spherical coordinates to cartesian coordinates
x = r_gen*np.sin(phi_gen)*np.cos(theta_gen)
y = r_gen*np.sin(phi_gen)*np.sin(theta_gen)
z = r_gen*np.cos(phi_gen)

# View distance
axis_limit = np.mean(r_gen) + (max(r_gen) - np.mean(r_gen))/2

# Radial color gradient
red_col = (1 - r_gen/max(r_gen))
green_col = np.full_like(r_gen, 1)
blue_col = r_gen/max(r_gen)

color = np.zeros((nb_points,3))
color[:, 0] = red_col
color[:, 1] = green_col
color[:, 2] = blue_col

color_gradient = np.array([(1-i/255, 1, i/255) for i in range(256)])
radius_colormap = ListedColormap(color_gradient)

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0,axis_limit/a_0), cmap=radius_colormap),
                    ax=ax, orientation='vertical', label=r"$r(a_0)$")

cbar.set_label(label = r"$r$ $(a_0)$", fontsize = 15)

# Plot finalization
ax.scatter(x, y, z, s=0.025, c=color)

plt.title(r"$P(X\in V)=\iiint_V |\Psi_{nlm}(r,\theta,\phi)|^2 r^2\sin (\phi) dr d\theta d\phi$"+
          '\n'+r"$n="+str(n)+", l="+str(l)+", m="+str(m)+"$", loc="center", fontsize=18)
ax.set_zlim(-0.75*axis_limit,0.75*axis_limit)
ax.set_xlim(-axis_limit,axis_limit)
ax.set_ylim(-axis_limit,axis_limit)

plt.show()