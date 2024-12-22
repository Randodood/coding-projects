import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import sparse

# User inputs
try: 
    time_interval = float(input('Enter time interval (recommended: less than 5): '))
except:
    raise ValueError('Must be a number')
if (time_interval <= 0):
    raise ValueError('Negative time interval')
try:
    nb_frames = int(input('Enter number of time frames : '))
    eigenstates_resolution = int(input('Enter number of eigenstates (recommended: 150): '))
    resolution = int(input('Enter discretization resolution (recommended: 100): '))
    contourf_resolution = int(input('Enter graph resolution (recommended: 100): '))
except:
    raise ValueError('Non-integer value')
try:
    x_position = float(input('Enter initial x position (-30 < x < 30): '))
    y_position = float(input('Enter initial y position (-30 < y < 30): '))
    stdev = float(input('Enter initial standard deviation: '))
    kx0 = float(input('Enter initial x momentum: '))
    ky0 = float(input('Enter initial y momentum: '))
except:
    raise ValueError('Must be a number')
if (x_position <= -30 or x_position >= 30 or y_position <= -30 or y_position >= 30):
    raise ValueError('Outside of domain')
if (stdev <= 0):
    raise ValueError('Standard deviation must be positive')

# Physical constants
time = np.linspace(0, time_interval, nb_frames)
xa = -30
xb = 30
ya = -30
yb = 30

h_bar = 1
m = 1

x, y = np.meshgrid(np.linspace(xa, xb, resolution), np.linspace(ya, yb, resolution))
x_for_calc, y_for_calc = np.meshgrid(np.linspace(xa, xb, resolution)[1:-1], np.linspace(ya, yb, resolution)[1:-1])
mat_size = resolution -2

delta_x = (xb - xa)/resolution
delta_y = (yb - ya)/resolution

# Potential functions
def harmonic(X, Y):
    return 0.1*np.sqrt(X**2 + Y**2)

def free(X, Y):
    return 0*(X + Y)

def double_slit(X, Y):
    temp = np.full_like(X, 0)
    for p in range(len(X)):
        for q in range((len(Y))):
            if np.abs(X[p][q]) < 1 and np.abs(Y[p][q]) < 2.5:
                temp[p][q] = 1000
            elif np.abs(X[p][q]) < 1 and np.abs(Y[p][q]) > 7.5:
                temp[p][q] = 1000
    return temp

def barrier(X, Y):
    temp = np.full_like(X, 0)
    for i in range(len(X)):
        for j in range(len(X)):
            if np.abs(X[i][j]) < 1:
                temp[i][j] = 25
    return temp

def circular(X, Y):
    temp = np.full_like(X, 0)
    a = xa
    b = ya
    for i in range(len(X)):
        for j in range((len(Y))):
            if (X[i][j]**2)/(a**2) + (Y[i][j]**2)/(b**2) > 1:
                temp[i][j] = 1000
    return temp

potentials = {"harmonic": harmonic, "free": free, "double_slit": double_slit, "barrier": barrier, "circular": circular}
potential = input('Enter potential function\nBuilt in functions: '+str(list(potentials.keys()))+": ")
if (potential not in potentials.keys()):
    raise ValueError('Invalid potential')
current_potential = potentials[potential]

# Initial gaussian wavefunction
p0 = np.exp(-0.5*((x-x_position)**2+(y-y_position)**2)/(stdev**2)) 
p0 = p0/(np.sum(p0)*delta_x*delta_y)
psi0 = np.sqrt(p0)*np.exp(-1j*(kx0*x + ky0*y))

# Discretization of the Hamiltonian operator in 2D
def get_V(func, X, Y):
    return func(X, Y)

diag = np.ones(mat_size)
diags = np.array([diag, -2*diag, diag])
Dx = sparse.spdiags(diags/(delta_x**2), np.array([-1, 0, 1]), mat_size, mat_size)
Dy = sparse.spdiags(diags/(delta_y**2), np.array([-1, 0, 1]), mat_size, mat_size)
T = -((h_bar**2)/(2*m))*sparse.kronsum(Dy, Dx)

V = get_V(current_potential, x_for_calc, y_for_calc)
U = sparse.diags(V.reshape(mat_size**2), (0))

H = T + U

# Solving for eigenvalues and for eigenvectors
eigenvalues, eigenvectors = sparse.linalg.eigsh(H, k=eigenstates_resolution, which = 'SM')
energies_normalized = eigenvalues/eigenvalues[0]

# Eigenvector coefficients for approxiating initial wavefunction
coefs = []
for i in range(eigenstates_resolution):
    coef = np.sum(eigenvectors.T[i].reshape(mat_size, mat_size)*psi0[1:-1,1:-1])*delta_x*delta_y
    coefs.append(coef)

# Plot Initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(8,7.5))
potential_3D = fig.add_subplot(2, 2, 1)
psi_squared = fig.add_subplot(2,2,2)
psi_real = fig.add_subplot(2, 2, 3)
psi_imag = fig.add_subplot(2, 2, 4)
plt.subplots_adjust(wspace = 0.25, hspace = 0.5)

#Plot animation
def init():
    potential_3D.set_title(r'$V(x,y)$')
    potential_3D.set_xlim(xa, xb)
    potential_3D.set_xlabel(r'$x$')
    potential_3D.set_ylabel(r'$y$')
    potential_3D.set_ylim(ya, yb)
    potential_3D.contourf(x, y, get_V(current_potential, x, y), 100,  cmap = 'plasma')

    psi = np.full_like(x, 0, dtype = np.complex128)
    for j in range(eigenstates_resolution):
        psi[1:-1,1:-1] += coefs[j]*eigenvectors.T[j].reshape(mat_size, mat_size)
    psi2 = psi_squared.contourf(x, y, (psi*np.conj(psi)).real, 100, cmap='inferno')
    psi_real.contourf(x, y, psi.real, contourf_resolution, cmap='viridis')
    psi_imag.contourf(x, y, psi.imag, contourf_resolution, cmap='magma')
    return psi2

def update_frames(frame):
    psi_squared.clear()
    psi_squared.set_title(r'$|\Psi(x,y,t)|^2$')
    psi_squared.set_xlim(xa, xb)
    psi_squared.set_xlabel(r'$x$')
    psi_squared.set_ylabel(r'$y$')
    psi_squared.set_ylim(ya, yb)
    
    psi_real.clear()
    psi_real.set_title(r'$\Re(\Psi(x,y,t))$')
    psi_real.set_xlim(xa, xb)
    psi_real.set_xlabel(r'$x$')
    psi_real.set_ylabel(r'$y$')
    psi_real.set_ylim(ya, yb)
    
    psi_imag.clear()
    psi_imag.set_title(r'$\Im(\Psi(x,y,t))$')
    psi_imag.set_xlim(xa, xb)
    psi_imag.set_xlabel(r'$x$')
    psi_imag.set_ylabel(r'$y$')
    psi_imag.set_ylim(ya, yb)
    
    psi = np.full_like(x, 0, dtype = np.complex128)
    for k in range(eigenstates_resolution):
        psi[1:-1,1:-1] += coefs[k]*eigenvectors.T[k].reshape(mat_size, mat_size)*np.exp(-1j*energies_normalized[k]*frame/h_bar)
    psi2 = psi_squared.contourf(x, y, (psi*np.conj(psi)).real, 100, cmap='inferno')
    psi_real.contourf(x, y, psi.real, contourf_resolution, cmap='viridis')
    psi_imag.contourf(x, y, psi.imag, contourf_resolution, cmap='magma')
    return psi2

ani = FuncAnimation(fig, update_frames, frames = time, interval=1, init_func = init, blit=False)

plt.show()