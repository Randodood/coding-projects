import matplotlib.pyplot as plt
import numpy as np
import tqdm as tqdm
from matplotlib.widgets import Slider, RadioButtons

# Spatial and temporal discrete steps
h = 0.001
dt = 0.001

# Initial path positions and velocities, and acceleration lists
X_list = [[-1], [-1], [-1], [-1], [-0.6], [-0.2], [0.2], [0.6]]
Y_list = [[0.6], [0.2], [-0.2], [-0.6], [-1], [-1], [-1], [-1]]
X_dot_list = [[1], [1], [1], [1], [0], [0], [0], [0]]
Y_dot_list = [[0], [0], [0], [0], [1], [1], [1], [1]]
X_ddot_list, Y_ddot_list = [[] for i in range(len(X_list))], [[] for i in range(len(X_list))]

# UI elements
colors_list = ['deepskyblue', 'deepskyblue', 'deepskyblue', 'deepskyblue', 'lime', 'lime', 'lime', 'lime']
geodesics_extrinsic_list, geodesics_intrinsic_list = [None for i in range(len(X_list))], [None for i in range(len(X_list))]
slider_states = [True, True, True, False, False]
equations = {r'$Plane$': r'$f(x, y) = \alpha(\beta x + \gamma y)$', 
             r'$Quadratic$': r'$f(x, y) = \alpha(\beta(x-\delta)^2 + \gamma(y-\epsilon)^2)$',
             r'$Gaussian$': r'$f(x, y) = \alpha e^{-(\beta(x-\delta)^2 + \gamma(y-\epsilon)^2)}$', 
             r'$Sinusoidal$': r'$f(x, y) = \alpha\sin(\beta\pi xy)$',
             r'$Radial$': r'$f(x, y) = \alpha\sqrt{2-\beta x^2-\gamma y^2}$'}
equations_slider_states = {r'$Plane$': [True, True, True, False, False], 
                           r'$Quadratic$': [True, True, True, True, True],
                           r'$Gaussian$': [True, True, True, True, True], 
                           r'$Sinusoidal$': [True, True, False, False, False],
                           r'$Radial$': [True, True, True, False, False]}

# Manifold x and y meshes
X_mesh, Y_mesh = np.meshgrid(np.linspace(-1, 1, 25), np.linspace(-1, 1, 25))

# Initial plane equation
def f(x, y):
    return x + y

# Assign equation with parameters
def get_f(label, a, b, c, d, e):
    if label == r'$Plane$':
        def f(x, y):
            return a*(b*x + c*y) + 0*d + 0*e
    elif label == r'$Quadratic$':
        def f(x, y):
            global slider_states
            return a*(b*(x-d)**2 + c*(y-e)**2)
    elif label == r'$Gaussian$':
        def f(x, y):
            return a*np.exp(-(b*(x-d)**2 + c*(y-e)**2))
    elif label == r'$Sinusoidal$':
        def f(x, y):
            return a*np.sin(b*np.pi*x*y) + 0*c + 0*d + 0*e
    elif label == r'$Radial$':
        def f(x, y):
            return a*np.sqrt(2-b*x**2-c*y**2) + 0*d + 0*e 
    return f

# Numerical approach to calculate geodesics
def get_metric_properties(F, x, y):
    
    del_x_F = (F(x+h,y) - F(x-h,y)) / (2*h)
    del_y_F = (F(x,y+h) - F(x,y-h)) / (2*h)
    del_xx_F = (F(x+2*h,y) - 2*F(x,y) + F(x-2*h,y)) / (2*h)**2
    del_yy_F = (F(x,y+2*h) - 2*F(x,y) + F(x,y-2*h)) / (2*h)**2
    del_xy_F = (F(x+h, y+h) + F(x-h, y-h) - F(x-h, y+h) - F(x+h, y-h)) / (2*h)**2
    
    g = np.array([[1 + del_x_F**2, del_x_F*del_y_F],
                  [del_x_F*del_y_F, 1 + del_y_F**2]])
    g_inv = np.linalg.inv(g)
    
    del_g = [[[2*del_x_F*del_xx_F, del_xx_F*del_y_F + del_x_F*del_xy_F],
              [del_xx_F*del_y_F + del_x_F*del_xy_F, 2*del_y_F*del_xy_F]],
             
             [[2*del_x_F*del_xy_F, del_xy_F*del_y_F + del_x_F*del_yy_F],
              [del_xy_F*del_y_F + del_x_F*del_yy_F, 2*del_y_F*del_yy_F]]]
    return g_inv.tolist(), del_g

def get_christoffel_symbol(k, i, j, g_inv, del_g):
    # Form: \Gamma^k_{ij}
    gamma_kij = 0
    for ksi in range(2):
        gamma_kij += 0.5*g_inv[ksi][k]*(del_g[i][j][ksi] + del_g[j][ksi][i] - del_g[ksi][i][j])
    return gamma_kij

def get_next_ddot(F, x, y, x_dot_now, y_dot_now):
    
    g_inv, del_g = get_metric_properties(F, x, y)
    gamma_xxx = get_christoffel_symbol(0, 0, 0, g_inv, del_g)
    gamma_yyy = get_christoffel_symbol(1, 1, 1, g_inv, del_g)
    gamma_xxy = get_christoffel_symbol(0, 0, 1, g_inv, del_g)
    gamma_yxy = get_christoffel_symbol(1, 1, 0, g_inv, del_g)
    gamma_xyy = get_christoffel_symbol(0, 1, 1, g_inv, del_g)
    gamma_yxx = get_christoffel_symbol(1, 0, 0, g_inv, del_g)
    
    x_ddot_next = -gamma_xxx*x_dot_now**2 - 2*gamma_xxy*x_dot_now*y_dot_now - gamma_xyy*y_dot_now**2
    y_ddot_next = -gamma_yxx*x_dot_now**2 - 2*gamma_yxy*x_dot_now*y_dot_now - gamma_yyy*y_dot_now**2
    return x_ddot_next, y_ddot_next
    
def step_integrate(func_now, func_dot_now):
    func_next = func_now + func_dot_now*dt
    return func_next

def get_geodesic_curve(F, X, Y, X_dot, Y_dot, X_ddot, Y_ddot):
    while np.abs(X[-1]) <= 1 and np.abs(Y[-1]) <= 1 and len(X) < 10000:
        x_ddot_next, y_ddot_next = get_next_ddot(F, X[-1], Y[-1], X_dot[-1], Y_dot[-1])
        X_ddot.append(x_ddot_next), Y_ddot.append(y_ddot_next)
        X_next, Y_next = step_integrate(X[-1], X_dot[-1]), step_integrate(Y[-1], Y_dot[-1])
        X.append(X_next), Y.append(Y_next)
        X_dot_next, Y_dot_next = step_integrate(X_dot[-1], X_ddot[-1]), step_integrate(Y_dot[-1], Y_ddot[-1])
        X_dot.append(X_dot_next), Y_dot.append(Y_dot_next)

# Z coordinate of geodesics
def get_Z_list(F, X, Y):
    Z_list = []
    for i in range(len(X)):
        Z_list.append(F(X[i], Y[i]))
    return Z_list

# Reset positions and velocities for different equation
def reset_lists(index):
    X_list[index], Y_list[index] = [X_list[index][0]], [Y_list[index][0]]
    X_dot_list[index], Y_dot_list[index] = [X_dot_list[index][0]], [Y_dot_list[index][0]] 
    X_ddot_list[index], Y_ddot_list[index] = [], []

# Plot initiation
plt.style.use('dark_background')
fig = plt.figure(figsize=(10, 7.5))

# Extrinsic view
ax_extrinsic = plt.subplot2grid((24,16), (3,0), colspan = 7, rowspan = 14, projection = '3d')
ax_extrinsic.set_title('Extrinsic View', color = 'w')
ax_extrinsic.set_xlim(-1, 1), ax_extrinsic.set_ylim(-1, 1)
ax_extrinsic.set_zlim(np.min(f(X_mesh, Y_mesh))-np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))), 
                          np.max(f(X_mesh, Y_mesh))+np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))))
ax_extrinsic.axis('off')

# Intrinsic view
ax_intrinsic = plt.subplot2grid((24,16), (3,9), colspan = 7, rowspan = 14)
ax_intrinsic.set_title('Intrinsic View', color = 'w')
ax_intrinsic.set_xlim(-1, 1), ax_intrinsic.set_ylim(-1, 1)
ax_intrinsic.set_xlabel(r'$x$'), ax_intrinsic.set_ylabel(r'$y$')
ax_intrinsic.set_xticks([-1, -0.5, 0, 0.5, 1]), ax_intrinsic.set_yticks([-1, -0.5, 0, 0.5, 1])

# Initial plotted geodesics
for i in tqdm.tqdm(range(len(X_list)), desc = "Initializing geodesic curves..."):
    get_geodesic_curve(f, X_list[i], Y_list[i], X_dot_list[i], Y_dot_list[i], X_ddot_list[i], Y_ddot_list[i])
    geodesics_extrinsic_list[i], = ax_extrinsic.plot(X_list[i], Y_list[i], get_Z_list(f, X_list[i], Y_list[i]), 
                                                     color = colors_list[i], linewidth = 0.75)
    geodesics_intrinsic_list[i], = ax_intrinsic.plot(X_list[i], Y_list[i], color = colors_list[i], linewidth = 1)
    reset_lists(i)

# Initial manifold
manifold = ax_extrinsic.plot_surface(X_mesh, Y_mesh, f(X_mesh, Y_mesh), cmap = "winter", alpha = 0.4)

# Radio buttons to select equation
ax_radio = plt.subplot2grid((24, 18), (19, 0), colspan = 3, rowspan = 5)
ax_radio.set_facecolor('k')
radio = RadioButtons(ax_radio, 
                     (r'$Plane$', r'$Quadratic$', r'$Gaussian$', r'$Sinusoidal$', r'$Radial$'),
                     radio_props = {'s': [64, 64, 64, 64, 64], 
                    'color': ['lightgoldenrodyellow', 'lightgoldenrodyellow', 'lightgoldenrodyellow', 'lightgoldenrodyellow', 'lightgoldenrodyellow']})

# Parameter sliders
ax_a = plt.subplot2grid((24,18), (19,4), colspan = 14)
a_slider = Slider(
    ax = ax_a,
    label = r'$\alpha$',
    valmin = -5,
    valmax = 5,
    valstep = 0.1,
    valinit = 1,
    facecolor = 'lightgoldenrodyellow'
)
a_slider.label.set_color('w')
a_slider.label.set_size(12.5)
a_slider.valtext.set_size(12.5)
a_slider.valtext.set_color('w')

ax_b = plt.subplot2grid((24,18), (20,4), colspan = 14)
b_slider = Slider(
    ax = ax_b,
    label = r'$\beta$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 1,
    facecolor = 'lightgoldenrodyellow'
)
b_slider.label.set_color('w')
b_slider.label.set_size(12.5)
b_slider.valtext.set_size(12.5)
b_slider.valtext.set_color('w')

ax_c = plt.subplot2grid((24,18), (21,4), colspan = 14)
c_slider = Slider(
    ax = ax_c,
    label = r'$\gamma$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 1,
    facecolor = 'lightgoldenrodyellow'
)
c_slider.label.set_color('w')
c_slider.label.set_size(12.5)
c_slider.valtext.set_size(12.5)
c_slider.valtext.set_color('w')

ax_d = plt.subplot2grid((24,18), (22,4), colspan = 14)
d_slider = Slider(
    ax = ax_d,
    label = r'$\delta$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'gray',
)
d_slider.set_active(False)
d_slider.label.set_color('gray')
d_slider.label.set_size(12.5)
d_slider.valtext.set_size(12.5)
d_slider.valtext.set_color('gray')

ax_e = plt.subplot2grid((24,18), (23,4), colspan = 14)
e_slider = Slider(
    ax = ax_e,
    label = r'$\epsilon$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'gray',
)
e_slider.set_active(False)
e_slider.label.set_color('gray')
e_slider.label.set_size(12.5)
e_slider.valtext.set_size(12.5)
e_slider.valtext.set_color('gray')

sliders = [a_slider, b_slider, c_slider, d_slider, e_slider]

# Display current equations
ax_equation = plt.subplot2grid((24,18), (0,0), colspan = 18, rowspan = 2)
ax_equation.set_xlim(-1, 1), ax_equation.set_ylim(-1, 1)
ax_equation.axis('off')
equation = ax_equation.text(0, 1, r'$f(x, y) = \alpha(\beta x + \gamma y)$', color = 'lightgoldenrodyellow', 
                            ha = 'center', va = 'center', fontsize = 20)

# Update plots and UI based on input
def update_slider(val):
    for i in range(len(X_list)):
        reset_lists(i)
        f = get_f(radio.value_selected, a_slider.val, b_slider.val, c_slider.val, d_slider.val, e_slider.val)
        get_geodesic_curve(f, X_list[i], Y_list[i], X_dot_list[i], Y_dot_list[i], X_ddot_list[i], Y_ddot_list[i])
        geodesics_extrinsic_list[i].set_data(X_list[i], Y_list[i])
        geodesics_extrinsic_list[i].set_3d_properties(get_Z_list(f, X_list[i], Y_list[i]))
        geodesics_intrinsic_list[i].set_data(X_list[i], Y_list[i])
    global manifold
    manifold.remove()
    manifold = ax_extrinsic.plot_surface(X_mesh, Y_mesh, f(X_mesh, Y_mesh), cmap = "winter", alpha = 0.4)
    ax_extrinsic.set_zlim(np.min(f(X_mesh, Y_mesh))-np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))), 
                          np.max(f(X_mesh, Y_mesh))+np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))))
        
def change_slider_state(index):
    if not slider_states[index]:
        sliders[index].label.set_color('gray')
        sliders[index].valtext.set_color('gray')
        sliders[index].poly.set_fc('gray')
        sliders[index].set_active(False)
    else:
        sliders[index].label.set_color('w')
        sliders[index].valtext.set_color('w')
        sliders[index].poly.set_fc('lightgoldenrodyellow')
        sliders[index].set_active(True)
        
def update_radio(label):
    for i in range(len(X_list)):
        equation.set_text(equations[label])
        reset_lists(i)
        f = get_f(label, a_slider.val, b_slider.val, c_slider.val, d_slider.val, e_slider.val)
        get_geodesic_curve(f, X_list[i], Y_list[i], X_dot_list[i], Y_dot_list[i], X_ddot_list[i], Y_ddot_list[i])
        geodesics_extrinsic_list[i].set_data(X_list[i], Y_list[i])
        geodesics_extrinsic_list[i].set_3d_properties(get_Z_list(f, X_list[i], Y_list[i]))
        geodesics_intrinsic_list[i].set_data(X_list[i], Y_list[i])
    global slider_states
    slider_states = equations_slider_states[label]
    for i in range(len(sliders)):
        change_slider_state(i)
    global manifold
    manifold.remove()
    manifold = ax_extrinsic.plot_surface(X_mesh, Y_mesh, f(X_mesh, Y_mesh), cmap = "winter", alpha = 0.4)
    ax_extrinsic.set_zlim(np.min(f(X_mesh, Y_mesh))-np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))), 
                          np.max(f(X_mesh, Y_mesh))+np.exp(-(np.max(f(X_mesh, Y_mesh)-np.min(f(X_mesh, Y_mesh))))))
    fig.canvas.draw_idle()
        
a_slider.on_changed(update_slider)
b_slider.on_changed(update_slider)
c_slider.on_changed(update_slider)
d_slider.on_changed(update_slider)
e_slider.on_changed(update_slider)
radio.on_clicked(update_radio)

plt.show()