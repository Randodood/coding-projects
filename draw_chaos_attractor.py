import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random

try:
    nb_of_lines = int(input('Enter number of points: '))
    plot_speed = int(input('Enter integer point speed (recommended: 20): '))
    resolution = int(input('Enter time resolution (recommended: 1000): '))
    tail_length = int(input('Enter integer trail length: '))
except:
    raise ValueError('Invalid integers')
if nb_of_lines <= 0 or plot_speed <= 0 or resolution <= 0 or tail_length <= 0:
    raise ValueError('Invalid values')
try:
    disp = float(input('Enter dispersion of generated points: '))
except:
    raise ValueError('Not a number')
attractors = ['lorenz', 'dadras', 'finance', 'nose-hoover', 'tsucs', 'wang-sun', 'lorenz83', 'chen', 'thomas']
attractor = input('Enter name of attractor\nBuilt in attractors: '+str(attractors)+": ")
if attractor not in attractors:
    raise ValueError('Unknown attractor')

# Time-step and center
delta_t = 1/resolution
center = [0, 0, 0]

# Attractor differential equations
if attractor == 'lorenz':
    title = 'Lorenz Attractor'
    params = [10, 28, 8/3]
    subtitle = r'$\sigma=$'+str(params[0])+r', $\rho=$'+str(params[1])+r', $\beta=$'+str(round(params[2],3))
    def dxdt(xs_n, t_n, fparams):
        return fparams[0]*(xs_n[1] - xs_n[0])

    def dydt(xs_n, t_n, fparams):
        return xs_n[0]*(fparams[1] - xs_n[2]) - xs_n[1]

    def dzdt(xs_n, t_n, fparams):
        return xs_n[0]*xs_n[1] - fparams[2]*xs_n[2]

if attractor == 'dadras':
    title = 'Dadras Attractor'
    params = [3, 27, 1.7, 2, 9]
    subtitle ='a='+str(params[0])+', b='+str(params[1])+', c='+str(params[2])+', d='+str(params[3])+', e='+str(params[4])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e = fparams
        return y - a*x + b*y*z

    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e = fparams
        return c*y - x*z + z

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e = fparams
        return d*x*y - e*z
    
if attractor == 'finance':
    title = 'Finance Attractor'
    params = [0.0001, 0.2, 1.1]
    subtitle ='a='+str(params[0])+', b='+str(params[1])+', c='+str(params[2])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c = fparams
        return (-a + 1/b)*x + z + x*y
    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c = fparams
        return -b*y - x**2

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c = fparams
        return -x-c*z
    
if attractor == 'nose-hoover':
    title = 'Nosé-Hoover Attractor'
    params = [1.5]
    subtitle ='a='+str(params[0])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a = fparams[0]
        return y

    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a = fparams[0]
        return -x + y*z

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a = fparams[0]
        return a - y**2
    
if attractor == 'tsucs':
    title = 'Three-Scroll Unified Chaotic System'
    params = [32.48, 45.84, 1.18, 0.13, 0.57, 14.7]
    subtitle ='a='+str(params[0])+', b='+str(params[1])+', c='+str(params[2])+', d='+str(params[3])+', e='+str(params[4])+', f='+str(params[5])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return a*(y - x) + d*x*z

    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return b*x - x*z + f*y

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return c*z + x*y - e*x**2
    
if attractor == 'wang-sun':
    title = 'Wang-Sun Attractor'
    params = [0.2, 0.01, 1, -0.4, -1, -1]
    subtitle ='a='+str(params[0])+', b='+str(params[1])+', c='+str(params[2])+', d='+str(params[3])+', e='+str(params[4])+', f='+str(params[5])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return a*x + c*y*z
    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return b*x + d*y - x*z
    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, c, d, e, f = fparams
        return e*z + f*x*y
    
if attractor == 'lorenz83':
    title = 'Lorenz83 Attractor'
    params = [0.95, 7.91, 4.83, 4.66]
    subtitle ='a='+str(params[0])+', b='+str(params[1])+', f='+str(params[2])+', g='+str(params[3])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, f, g = fparams
        return -a*x - y**2 - z**2 + a*f
    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, f, g = fparams
        return -y + x*y - b*x*z + g

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, f, g = fparams
        return -z + b*x*y + x*z
    
if attractor == 'chen':
    title = 'Chen Attractor'
    params = [5, -10, -0.38]
    subtitle = r'$\alpha$='+str(params[0])+r', $\beta=$'+str(params[1])+', $\delta=$'+str(params[2])
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, d = fparams
        return a*x - y*z
    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, d = fparams
        return b*y + x*z

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        a, b, d = fparams
        return d*z + x*y/3

if attractor == 'thomas':
    title = 'Thomas Attractor'
    params = [0.208186]
    subtitle = 'b='+str(params[0])
    from math import sin
    def dxdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        b, = fparams
        return sin(y) - b*x

    def dydt(xs_n, t_n, fparams):
        x, y, z = xs_n
        b, = fparams
        return sin(z) - b*y

    def dzdt(xs_n, t_n, fparams):
        x, y, z = xs_n
        b, = fparams
        return sin(x) - b*z

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(8,7.5))
ax = fig.add_subplot(projection='3d')
fig.suptitle(title, color = 'w', fontsize = 25)
low_index = 0

# Colors
colors = {0:'b', 1:'g', 2:'r', 3:'c', 4:'m', 5:'y', 6:'w'}

# Generation of random points
lines = [[([center[0] + (random.random()-0.5)*2*disp], [center[1] + (random.random()-0.5)*2*disp],
           [center[2] + (random.random()-0.5)*2*disp]),[0], random.randint(0, 6)] for i in range(nb_of_lines)]

# Runge-Kutta method
def rk4_nd(xs: tuple, t: list, fs: tuple, delta_t: float, fparams: list = []) -> None:
    xs_n = np.array([xs[i][-1] for i in range(len(xs))])
    t_n = t[-1]
    k1 = np.array([f(xs_n, t_n, fparams) for f in fs])
    k2 = np.array([f(xs_n + delta_t*k1/2, t_n + delta_t/2, fparams) for f in fs])
    k3 = np.array([f(xs_n + delta_t*k2/2, t_n + delta_t/2, fparams) for f in fs])
    k4 = np.array([f(xs_n + delta_t*k3, t_n + delta_t, fparams) for f in fs])
    
    xs_new = xs_n + (k1 + 2*k2 + 2*k3 + k4)*delta_t/6
    t_new = t_n + delta_t
    for i in range(len(xs)):
        xs[i].append(xs_new[i])
    t.append(t_new)

# Plot animation
def init():
    ax.set_axis_off()
    ax.set_title(subtitle, color = 'w', fontsize = 20, y = 1.1)
    if len(lines[0][0][0]) - tail_length < 0:
        low_index = 0
    else:
        low_index = len(lines[0][0][0]) - tail_length
    for line in lines:
        ax.plot(line[0][0][low_index:-1], line[0][1][low_index:-1], line[0][2][low_index:-1], color = colors[line[2]])
    
def update_frames(frame):
    for line in lines:
        for i in range(plot_speed):
            rk4_nd(line[0], line[1], (dxdt, dydt, dzdt), delta_t, params)
    ax.clear()
    init()

ani = FuncAnimation(fig, update_frames, frames=range(1), interval=1, init_func = init, blit=False)

plt.show()