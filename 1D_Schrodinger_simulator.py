import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches

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
    resolution = int(input('Enter resolution (recommended: 1000): '))
except:
    raise ValueError('Non-integer value')
try:
    position = float(input('Enter initial position (-30 < x < 30): '))
    stdev = float(input('Enter initial standard deviation: '))
    k0 = float(input('Enter initial momentum: '))
except:
    raise ValueError('Must be a number')
if (position <= -30 or position >= 30):
    raise ValueError('Outside of domain')
if (stdev <= 0):
    raise ValueError('Standard deviation must be positive')

# Physical constants
h_bar = 1
m = 1
left_lim = -30
right_lim = 30
y_a, y_b = -0.75, 1.25 

x = np.linspace(left_lim, right_lim, resolution)
X = x[1:-1]
delta_x = X[1] - X[0]

# Potential functions
def harmonic(x):
    return 0.05*(x**2)

def free(x):
    return np.full_like(x, 0)

def barrier(x):
    temp = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] > -2.5 and x[i] < 2.5:
            temp[i] = 15
    return temp

potentials = {"harmonic": harmonic, "free": free, "barrier": barrier}
potential = input('Enter potential function\nBuilt in functions: '+str(list(potentials.keys()))+": ")
if (potential not in potentials.keys()):
    raise ValueError('Invalid potential')
current_function = potentials[potential]

# Discretization of Hamiltonian operator 
def get_V_vector(func):
    return np.array(func(x))

def get_V_matrix(func):
    return np.diag(func(X))

def get_T_matrix():
    const = -(h_bar**2)/(2*m*delta_x**2)
    one_temp = list(np.zeros((len(X), len(X))))
    for i in range(len(X)):
        if i+1 < len(X):
            one_temp[i][i+1] = 1
    T = np.array(one_temp) + np.array(one_temp).T -2*np.eye(len(X))
    return const*T

def get_Hamiltonian(func):
    return get_T_matrix() + get_V_matrix(func)

# Finding eigenvalues and eigenvectors of the Hamiltonian
val, vec = np.linalg.eig(get_Hamiltonian(current_function))
z = np.argsort(val)
z = z[0:eigenstates_resolution]
energies_unit_E0 = (val[z]/val[z][0])

# Gaussian initial probability function
p0 = np.exp(-0.5*((x-position)/stdev)**2)
p0 = p0/(np.dot(p0,np.ones(len(p0)))*delta_x)

# Initial wavefunction
psi0 = np.sqrt(p0)*np.exp(1j*k0*x)

# Normalization of eigenvectors
eigenstates = []
for i in range(len(z)):
    y = []
    y = np.append(y, vec[:,z[i]])
    y = np.append(y, 0)
    y = np.insert(y, 0, 0)
    y = y/np.sqrt(np.dot(y,y)*delta_x)
    eigenstates.append(y)

# Eigenvector coefficients for approxiating initial wavefunction
coef = []
for j in range(len(eigenstates)):
    coef.append(np.dot(eigenstates[j],psi0)*delta_x)

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(8,7.5))
plt.xlim(left_lim, right_lim)
plt.ylim(y_a, y_b)
plt.xlabel(r'$x$')
plt.ylabel(r'$Amplitude$')
yellow_patch = mpatches.Patch(color='yellow', label=r'$V(x)$')
black_patch = mpatches.Patch(color='white', label=r'$|\Psi(x,t)|^2$')
blue_patch = mpatches.Patch(color='blue', label=r'$\Re[\Psi(x,t)]$')
red_patch = mpatches.Patch(color='red', label=r'$\Im[\Psi(x,t)]$')
patches = [yellow_patch, black_patch, blue_patch, red_patch]

# Plot animation
def init():
    plt.plot(x, get_V_vector(current_function), color = 'y')
    plt.legend(handles = patches)
    plt.xlim(left_lim, right_lim)
    psi = np.zeros(len(eigenstates[0]), dtype = np.complex128)
    for k in range(len(eigenstates)):
        psi += coef[k]*eigenstates[k]
    plt.plot(x, psi.imag, color = 'r')
    plt.plot(x, psi.real, color = 'b')
    psi2 = plt.plot(x, (psi**2).real, color = 'k')
    return psi2
    
def update_psi(frame):
    plt.clf()
    plt.plot(x, get_V_vector(current_function), color = 'y')
    plt.xlim(left_lim, right_lim)
    plt.ylim(y_a, y_b)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$Amplitude$')
    plt.legend(handles = patches)
    psi = np.zeros(len(eigenstates[0]), dtype = np.complex128)
    for l in range(len(eigenstates)):
        psi += coef[l]*eigenstates[l]*np.exp(-1j*energies_unit_E0[l]*frame/h_bar)
    plt.plot(x, psi.real, color = 'b')
    plt.plot(x, psi.imag, color = 'r')
    plt.plot(x, (psi*np.conj(psi)).real, color = 'w')
    
ani = FuncAnimation(fig, update_psi, frames=np.linspace(0, time_interval, nb_frames), interval=1, init_func = init, blit=False)

plt.show()