import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, TextBox
import numpy as np
from scipy.sparse import csr_array, lil_array, csr_matrix
from scipy.sparse.linalg import spsolve, inv
from matplotlib.patches import FancyArrow

# Discrete functions and operators
def get_quantities(L, n):
    dx = L/n
    vec_one, vec_x = np.ones(n), np.linspace(0, L, n)
    mat_derivative = lil_array((n, n))
    mat_derivative.setdiag(-1, -1), mat_derivative.setdiag(1, 0)
    mat_derivative = mat_derivative.tocsc()/dx
    mat_integral = inv(mat_derivative)
    return dx, vec_one, vec_x, mat_integral

def get_heaviside(vec_x, c):   
    return np.where(vec_x >= c, 1, 0)

def get_w(x, expression, x1, x2):
    limits = np.array([x1, x2])
    try:
        new_w_ydata = eval(expression, {'np': np}, {'x': x})
    except:
        new_w_ydata = 0*x
    if np.max(np.abs(new_w_ydata)) > 10:
        new_w_ydata = new_w_ydata/np.max(np.abs(new_w_ydata))*10
    truncated_y_data = new_w_ydata*((get_heaviside(x, np.min(limits)) - get_heaviside(x, np.max(limits))))
    w_disp_x = x[truncated_y_data != 0]
    w_disp_y = truncated_y_data[truncated_y_data != 0]
    w_disp.set_data(w_disp_x, w_disp_y)
    if len(w_disp_x) == 0 or len(w_disp_y) == 0 or w_disp_x[0] == w_disp_x[-1] or w_disp_y[0] == 0 or w_disp_y[-1] == 0:
        w_arrow_1.set_visible(False), w_arrow_2.set_visible(False)
    else:
        w_arrow_1.set_visible(True), w_arrow_2.set_visible(True)
        w_arrow_1.set_data(x = w_disp_x[0], y = w_disp_y[0], dx = 0, dy = -w_disp_y[0])
        w_arrow_2.set_data(x = w_disp_x[-1], y = w_disp_y[-1], dx = 0, dy = -w_disp_y[-1])
    return truncated_y_data

def draw_moment(M, cx, n, right = True, color = 'orangered', M_obj = None):
    angle = np.pi*M/10
    angles = np.linspace(-angle, angle, n)
    y_circum = 1.5*np.sin(angles)
    alpha = 1
    if M == 0:
        alpha = 0
        x_circum, y_circum, dy, dx = [0], [0], 0, 0
    else:
        dy = 0.01*np.sin(-(angle+ M*np.pi/(2*np.abs(M))))
        if right:
            x_circum = 0.75*np.cos(angles) + cx
            dx = 0.005*np.cos(-(angle + M*np.pi/(2*np.abs(M))))
        else:
            x_circum = -0.75*np.cos(angles) + cx
            dx =  -0.005*np.cos(-(angle + M*np.pi/(2*np.abs(M))))
    if M_obj is None:
        body, = ax_beam.plot(x_circum, y_circum, color = color, linewidth = 1, alpha = alpha)
        head = FancyArrow(x_circum[0], y_circum[0], dx, dy, color = color, head_width = 0.2, head_length = 0.2, 
                          width = 0.00025, alpha = alpha)
        ax_beam.add_patch(head)
        position = ax_beam.scatter(cx, 0, color = color, alpha = alpha, s = 10)
        return body, head, position
    else:
        M_obj[0].set_data(x_circum, y_circum), M_obj[0].set_color(color), M_obj[0].set_alpha(alpha)
        M_obj[1].set_data(x = x_circum[0], y = y_circum[0], dx = dx, dy = dy), M_obj[1].set_color(color)
        M_obj[1].set_alpha(alpha)
        M_obj[2].set_offsets([cx, 0]), M_obj[2].set_color(color), M_obj[2].set_alpha(alpha)

def draw_force(F, cx, color = 'deepskyblue', F_obj = None):
    alpha = 1
    if F == 0:
        alpha = 0
    if F_obj is None:
        F_arrow = FancyArrow(cx, F, 0, -F, color = color, head_width = 0.2, head_length = 0.2, width = 0.00025, 
                             alpha = alpha)
        ax_beam.add_patch(F_arrow)
        return F_arrow
    else:
        F_obj.set_data(x = cx, y = F, dx = 0, dy = -F), F_obj.set_color(color), F_obj.set_alpha(alpha)

def general_update():
    EI = EI_slider.val
    n = resolution_slider.val
    dx, vec_one, x, mat_integral = get_quantities(L_init, n)
    a_w, b_w = dist_force_limit_1_slider.val, dist_force_limit_2_slider.val
    w_expression = dist_force_magnitude_textbox.text_disp.get_text()
    M, x_M = moment_magnitude_slider.val, moment_position_slider.val
    F, x_F = force_magnitude_slider.val, force_position_slider.val
    vec_w = get_w(x, w_expression, a_w, b_w)
    draw_moment(M, x_M, n, M_obj = M_disp)
    draw_force(F, x_F, F_obj = F_disp)
    u_M, u_F = get_heaviside(x, x_M), get_heaviside(x, x_F)
    find_values = select_condition(conditions_radio.value_selected)
    Ay, M_A, By, M_B, C1, C2 = find_values(dx, EI, x, vec_one, mat_integral, vec_w, u_M, u_F, F, M, x_F, x_M, L_init)
    theta, nu = get_deflection_slope(EI, x, vec_one, mat_integral, vec_w, u_M, u_F, Ay, M_A, F, M, C1, C2)
    slope.set_data(x, theta), deflection.set_data(x, nu)
    ax_deflection.set_ylim(min(nu)-np.exp(-((max(nu)-min(nu))/2)**2), max(nu)+np.exp(-((max(nu)-min(nu))/2)**2))
    ax_slope.set_ylim(min(theta)-np.exp(-((max(theta)-min(theta))/2)**2), max(theta)+np.exp(-((max(theta)-min(theta))/2)**2))
    
def get_M_ext(vec_x, vec_one, mat_int, vec_w, u_M, u_F, Ay, M_A, F, M):
    return Ay*vec_x - M_A*vec_one - mat_int@mat_int@vec_w - F*mat_int@u_F + M*u_M   

def get_deflection_slope(EI, vec_x, vec_one, mat_int, vec_w, u_M, u_F, Ay, M_A, F, M, C1, C2):
    theta = (1/EI)*mat_int@get_M_ext(vec_x, vec_one, mat_int, vec_w, u_M, u_F, Ay, M_A, F, M) + C1*vec_one
    nu = mat_int@theta + C2*vec_one
    return theta, nu

def select_condition(condition):
    if condition == '$Fixed-Free$':
        def find_values(dx, EI, vec_x, vec_one, mat_int, vec_w, u_M, u_F, F, M, x_F, x_M, L):
            Ay = F +dx*vec_w@vec_one
            M_A = M + dx*vec_w@vec_x + x_F*F
            C1, C2, By, M_B = 0, 0, 0, 0
            return Ay, M_A, By, M_B, C1, C2
    elif condition == '$Fixed-Pinned$':
        pass
    return find_values
    
def textbox_update(expression):
    general_update()
    
def slider_update(val):
    general_update()

# UI
plt.style.use('dark_background')
fig = plt.figure(figsize=(10,7.5))

ax_beam = plt.subplot2grid((12,16), (0,0), colspan = 7, rowspan = 7)
ax_beam.set_xlim(-1, 11), ax_beam.set_ylim(-11, 11), ax_beam.set_xlabel(r'$x$'), ax_beam.set_xticks([0, 2, 4, 6, 8, 10])
ax_beam.set_title(r'$Beam$')

ax_deflection = plt.subplot2grid((12,16), (0,9), colspan = 7, rowspan = 3)
ax_deflection.set_xlim(0, 10), ax_deflection.set_ylabel(r'$\nu(x)$')
ax_deflection.set_title(r'$Deflection$')

ax_slope = plt.subplot2grid((12,16), (4,9), colspan = 7, rowspan = 3)
ax_slope.set_xlim(0, 10), ax_slope.set_xlabel(r'$x$'), ax_slope.set_ylabel(r'$\theta (x)$')
ax_slope.set_title(r'$Slope$')

ax_conditions = plt.subplot2grid((24, 18), (20, 2), colspan = 4, rowspan = 4)
ax_conditions.set_facecolor('k')
conditions_radio = RadioButtons(ax_conditions, 
                     (r'$Fixed-Free$', r'$Fixed-Pinned$', r'$Pinned-Pinned$', r'$Fixed-Fixed$'),
                     radio_props = {'s': [64, 64, 64, 64], 
                    'color': ['lightgoldenrodyellow', 'lightgoldenrodyellow', 'lightgoldenrodyellow', 'lightgoldenrodyellow']})

ax_resolution = plt.subplot2grid((24, 18), (20, 0), rowspan = 4)
resolution_slider = Slider(
    ax = ax_resolution,
    label = r'$n$',
    valmin = 10,
    valmax = 100,
    valstep = 1,
    valinit = 100,
    facecolor = 'lightgoldenrodyellow',
    orientation = 'vertical'
)

ax_EI = plt.subplot2grid((24, 18), (20, 1), rowspan = 4)
EI_slider = Slider(
    ax = ax_EI,
    label = r'$EI$',
    valmin = 100,
    valmax = 1000,
    valstep = 1,
    valinit = 100,
    facecolor = 'lightgoldenrodyellow',
    orientation = 'vertical'
)

ax_moment_magnitude = plt.subplot2grid((24, 18), (23, 7), colspan = 4)
moment_magnitude_slider = Slider(
    ax = ax_moment_magnitude,
    label = r'$M$',
    valmin = -10,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'orangered',
)
moment_magnitude_slider.label.set_color('orangered')
moment_magnitude_slider.valtext.set_color('orangered')

ax_moment_position = plt.subplot2grid((24, 18), (23, 14), colspan = 4)
moment_position_slider = Slider(
    ax = ax_moment_position,
    label = r'$x_M$',
    valmin = 0,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'red',
)
moment_position_slider.label.set_color('orangered')
moment_position_slider.valtext.set_color('orangered')

ax_force_magnitude = plt.subplot2grid((24, 18), (22, 7), colspan = 4)
force_magnitude_slider = Slider(
    ax = ax_force_magnitude,
    label = r'$F$',
    valmin = -10,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'deepskyblue',
)
force_magnitude_slider.label.set_color('deepskyblue')
force_magnitude_slider.valtext.set_color('deepskyblue')

ax_force_position = plt.subplot2grid((24, 18), (22, 14), colspan = 4)
force_position_slider = Slider(
    ax = ax_force_position,
    label = r'$x_F$',
    valmin = 0,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'deepskyblue',
)
force_position_slider.label.set_color('deepskyblue')
force_position_slider.valtext.set_color('deepskyblue')

ax_dist_force_limit_1 = plt.subplot2grid((24, 18), (21, 7), colspan = 4)
dist_force_limit_1_slider = Slider(
    ax = ax_dist_force_limit_1,
    label = r'$a_w$',
    valmin = 0,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'lime',
)
dist_force_limit_1_slider.label.set_color('lime')
dist_force_limit_1_slider.valtext.set_color('lime')

ax_dist_force_limit_2 = plt.subplot2grid((24, 18), (21, 14), colspan = 4)
dist_force_limit_2_slider = Slider(
    ax = ax_dist_force_limit_2,
    label = r'$b_w$',
    valmin = 0,
    valmax = 10,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'lime',
)
dist_force_limit_2_slider.label.set_color('lime')
dist_force_limit_2_slider.valtext.set_color('lime')

ax_dist_force_magnitude = plt.subplot2grid((24, 18), (20, 7), colspan = 11)
dist_force_magnitude_textbox = TextBox(
    ax = ax_dist_force_magnitude,
    label = r'$w$',
    initial = '10*np.sin(np.pi*x/10)',
    color = 'lime',
    textalignment='center'
)
dist_force_magnitude_textbox.label.set_color('lime')
dist_force_magnitude_textbox.text_disp.set_color('k')

EI_init = 100
n_init = 100
L_init = 10
condition = '$Fixed-Free$'

dx, vec_one, x, mat_integral = get_quantities(L_init, n_init)
w_init = np.full_like(x, 0)
w_disp, = ax_beam.plot([], [], color = 'lime', linewidth = 1)
w_arrow_1 = FancyArrow(0, 0, 0, 0, color = 'lime', head_width = 0.2, head_length = 0.2, width = 0.00025)
w_arrow_2 = FancyArrow(0, 0, 0, 0, color = 'lime', head_width = 0.2, head_length = 0.2, width = 0.00025)
ax_beam.add_patch(w_arrow_1), ax_beam.add_patch(w_arrow_2)
vec_w = get_w(x, '10*np.sin(np.pi*x/10)', 0, 0)
M_disp, F_disp = draw_moment(0, 0, n_init), draw_force(0, 0)
ax_beam.add_patch(F_disp)
M_A_disp = draw_moment(3.5, 0, n_init, right = False, color = 'lightgoldenrodyellow')
M_B_disp = draw_moment(3.5, 10, n_init, color = 'lightgoldenrodyellow')
M_A_disp[2].set_visible(False), M_B_disp[2].set_visible(False)
beam = ax_beam.plot([0, 10], [0, 0], color = 'w', linewidth = 2, alpha = 0.35, linestyle = 'dashed')
F_A_disp = FancyArrow(0, -5, 0, 3.75, color = 'lightgoldenrodyellow', head_width = 0.2, head_length = 0.2, width = 0.00025)
F_B_disp = FancyArrow(10, -5, 0, 3.75, color = 'lightgoldenrodyellow', head_width = 0.2, head_length = 0.2, width = 0.00025)
ax_beam.add_patch(F_A_disp), ax_beam.add_patch(F_B_disp)
ax_beam.text(0, 0, r'$A$', color = 'lightgoldenrodyellow', fontsize = 12.5, ha = 'center', va = 'center')
ax_beam.text(10, 0, r'$B$', color = 'lightgoldenrodyellow', fontsize = 12.5, ha = 'center', va = 'center')
u_M, u_F = get_heaviside(x, 8), get_heaviside(x, 2)
find_values = select_condition(condition)
Ay, M_A, By, M_B, C1, C2 = find_values(dx, EI_init, x, vec_one, mat_integral, vec_w, u_M, u_F, 0, 0, 0, 0, L_init)
theta, nu = get_deflection_slope(EI_init, x, vec_one, mat_integral, vec_w, u_M, u_F, Ay, M_A, 0, 0, C1, C2)
slope, = ax_slope.plot(x, theta, color = 'lightgoldenrodyellow')
deflection, = ax_deflection.plot(x, nu, color = 'lightgoldenrodyellow')
ax_deflection.set_ylim(min(nu)-np.exp(-((max(nu)-min(nu))/2)**2), max(nu)+np.exp(-((max(nu)-min(nu))/2)**2))
ax_slope.set_ylim(min(theta)-np.exp(-((max(theta)-min(theta))/2)**2), max(theta)+np.exp(-((max(theta)-min(theta))/2)**2))

resolution_slider.on_changed(slider_update)
EI_slider.on_changed(slider_update)
dist_force_limit_1_slider.on_changed(slider_update)
dist_force_limit_2_slider.on_changed(slider_update)
dist_force_magnitude_textbox.on_submit(textbox_update)
moment_magnitude_slider.on_changed(slider_update)
moment_position_slider.on_changed(slider_update)
force_magnitude_slider.on_changed(slider_update)
force_position_slider.on_changed(slider_update)
dist_force_magnitude_textbox.set_val('10*np.sin(np.pi*x/10)')

plt.show()