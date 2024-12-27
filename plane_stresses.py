import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import matplotlib.patches as mpatches

# Rotation matrix function
def rot_matrix(theta_rad):
    return np.array([[np.cos(theta_rad), np.sin(theta_rad)],[-np.sin(theta_rad), np.cos(theta_rad)]])

init_theta_deg = 0
wedge_p1 = np.array([1, 1])
wedge_p2 = np.array([1, -1])
wedge_p3 = np.array([-1, -1])
wedge_points = np.array([wedge_p1, wedge_p2, wedge_p3, wedge_p1])

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(10,7.5))

# Wedge shaped stress element
ax_elmt = plt.subplot2grid((12,16), (0,0), colspan = 7, rowspan = 7)
wedge, = ax_elmt.plot(wedge_points[:,0], wedge_points[:,1], color = 'w', alpha = 0.35, linestyle = 'dashed')
ax_elmt.set_xlim(-3, 3)
ax_elmt.set_ylim(-3, 3)
ax_elmt.axis(False)
ax_elmt.set_title('Stress Element')
ax_elmt.set_xticks([]), ax_elmt.set_yticks([])
x_dir = mpatches.FancyArrow(-3, -3, 1, 0, color = 'w', head_width = 0.15, head_length = 0.15)
ax_elmt.add_patch(x_dir)
ax_elmt.text(-1.5, -3, r'$x$', color = 'w', ha = 'center', va = 'center')
y_dir = mpatches.FancyArrow(-3, -3, 0, 1, color = 'w', head_width = 0.15, head_length = 0.15)
ax_elmt.add_patch(y_dir)
ax_elmt.text(-3, -1.5, r'$y$', color = 'w', ha = 'center', va = 'center')

# Initial stress values
norm_x = 0.0
norm_y = 0.0
shear = 0.0

# Stresses initially in x
pos_x = np.array([1.1,0])
u_sigma_x = np.array([1,0])
u_tau_xy = np.array([0,-1])
sigma_x = mpatches.FancyArrow(pos_x[0], pos_x[1], norm_x*u_sigma_x[0], norm_x*u_sigma_x[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(sigma_x)
tau_xy = mpatches.FancyArrow(pos_x[0], pos_x[1], shear*u_tau_xy[0], shear*u_tau_xy[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(tau_xy)
T_x = norm_x*(u_sigma_x + u_tau_xy)
traction_x = mpatches.FancyArrow(pos_x[0], pos_x[1], T_x[0], T_x[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(traction_x)

# Stresses ininitially in y
pos_y = np.array([0,-1.1])
u_sigma_y = np.array([0,-1])
u_tau_yx = np.array([1,0])
sigma_y = mpatches.FancyArrow(pos_y[0], pos_y[1], norm_y*u_sigma_y[0], norm_y*u_sigma_y[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(sigma_y)
tau_yx = mpatches.FancyArrow(pos_y[0], pos_y[1], shear*u_tau_yx[0], shear*u_tau_yx[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(tau_yx)
T_y = norm_y*(u_sigma_y + u_tau_yx)
traction_y = mpatches.FancyArrow(pos_y[0], pos_y[1], T_y[0], T_y[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(traction_y)

# Stress transformations
# uses the counter-clockwise positive convention for shear stress
def get_sigma_n(s_x, s_y, t_xy, theta_rad):
    return 0.5*(s_x + s_y) + 0.5*(s_x-s_y)*np.cos(2*theta_rad) - t_xy*np.sin(2*theta_rad)

def get_tau_nt(s_x, s_y, t_xy, theta_rad):
    return -0.5*(s_x - s_y)*np.sin(2*theta_rad) - t_xy*np.cos(2*theta_rad)

# Stresses on the diagonal surface
pos_d = 0.1*np.array([-1,1])/np.sqrt(2)
u_sigma_d = np.array([-1,1])
u_tau_d = np.array([1,1])
norm_d = get_sigma_n(norm_x, norm_y, shear, np.radians(init_theta_deg)+5*np.pi/4)
shear_d = get_tau_nt(norm_x, norm_y, shear, np.radians(init_theta_deg)+5*np.pi/4)
sigma_d = mpatches.FancyArrow(pos_d[0], pos_d[1], norm_d*u_sigma_d[0], norm_d*u_sigma_d[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(sigma_d)
tau_d = mpatches.FancyArrow(pos_d[0], pos_d[1], shear_d*u_tau_d[0], shear_d*u_tau_d[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(tau_d)
T_d = norm_d*u_sigma_d + shear_d*u_tau_d
traction_d = mpatches.FancyArrow(pos_d[0], pos_d[1], T_d[0], T_d[1], color = 'k', head_width = 0.1, head_length = 0.1)
ax_elmt.add_patch(traction_d)

# Rotation slider
ax_rot = plt.subplot2grid((24,9), (22,0), colspan = 9, rowspan = 2)
rot_slider = Slider(
    ax = ax_rot,
    label = r'$\theta$ [deg]',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = np.radians(init_theta_deg),
    facecolor = 'lightgoldenrodyellow')
rot_slider.label.set_color('lightgoldenrodyellow')
rot_slider.label.set_size(20)
rot_slider.valtext.set_size(20)
rot_slider.valtext.set_color('lightgoldenrodyellow')

# Stress sliders
ax_sigma_x = plt.subplot2grid((24,18), (19,0), colspan = 4, rowspan = 2)
sigma_x_slider = Slider(
    ax = ax_sigma_x,
    label = r'$\sigma_x$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'lime'
)
sigma_x_slider.label.set_color('lime')
sigma_x_slider.label.set_size(20)
sigma_x_slider.valtext.set_size(20)
sigma_x_slider.valtext.set_color('lime')

ax_sigma_y = plt.subplot2grid((24,18), (19,7), colspan = 4, rowspan = 2)
sigma_y_slider = Slider(
    ax = ax_sigma_y,
    label = r'$\sigma_y$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'lime'
)
sigma_y_slider.label.set_color('lime')
sigma_y_slider.label.set_size(20)
sigma_y_slider.valtext.set_size(20)
sigma_y_slider.valtext.set_color('lime')

ax_tau_xy = plt.subplot2grid((24,18), (19,14), colspan = 4, rowspan = 2)
tau_xy_slider = Slider(
    ax = ax_tau_xy,
    label = r'$\tau_{xy}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'deepskyblue'
)
tau_xy_slider.label.set_color('deepskyblue')
tau_xy_slider.label.set_size(20)
tau_xy_slider.valtext.set_size(20)
tau_xy_slider.valtext.set_color('deepskyblue')

# Mohr circle
def mohr_circle(s_x, s_y, t_xy):
    C = np.array([0.5*(s_x + s_y), 0])
    R = np.array([np.sqrt(0.25*(s_x - s_y)**2 + t_xy**2),0])
    points = np.array([C + R@rot_matrix(theta) for theta in np.linspace(0, 2*np.pi, 100)]+[C+R])
    return points

ax_mohr = plt.subplot2grid((12,16), (0,9), colspan = 7, rowspan = 7)
ax_mohr.set_xlim(-2, 2), ax_mohr.set_ylim(-2, 2)
ax_mohr.set_xlabel(r'$\sigma$'), ax_mohr.set_ylabel(r'$\tau$')
ax_mohr.set_xticks([-2, -1, 0, 1, 2]), ax_mohr.set_yticks([-2, -1, 0, 1, 2])
ax_mohr.set_title('Mohr Circle')
init_mohr_circle = mohr_circle(norm_x, norm_y, shear)
circle, = ax_mohr.plot(init_mohr_circle[:,0], init_mohr_circle[:,1], color = 'w', alpha = 0.35, linestyle = 'dashed')
nt_line, = ax_mohr.plot([norm_x,norm_y], [shear, -shear], color = 'w', alpha = 0.35)

# Face identifiers
alpha_elmt = ax_elmt.text(0.75*pos_x[0], 0.75*pos_x[1], r'$\alpha$', color = 'y', ha = 'center', va = 'center')
beta_elmt = ax_elmt.text(0.75*pos_y[0], 0.75*pos_y[1], r'$\beta$', color = 'y', ha = 'center', va = 'center')
alpha_mohr = ax_mohr.text(norm_x, -shear, '', color = 'y', ha = 'center', va = 'center')
beta_mohr = ax_mohr.text(norm_y, shear, '', color = 'y', ha = 'center', va = 'center')
center = ax_mohr.scatter((norm_x + norm_y)/2, 0, color = 'y', s= 5)

# Stress display
ax_sigma_alpha = plt.subplot2grid((24, 18), (16, 0), colspan = 4, rowspan = 3)
ax_sigma_alpha.set_xlim(-1, 1), ax_sigma_alpha.set_ylim(-1, 1)
sigma_alpha_text = ax_sigma_alpha.text(0, 0, r'$\sigma_{\alpha}=$'+str(norm_x), color = 'lime', ha = 'center', va = 'center', fontsize = 20)
ax_sigma_alpha.set_xticks([]), ax_sigma_alpha.set_yticks([])
ax_sigma_alpha.axis(False)
ax_sigma_beta = plt.subplot2grid((24, 18), (16, 7), colspan = 4, rowspan = 3)
ax_sigma_beta.set_xlim(-1, 1), ax_sigma_beta.set_ylim(-1, 1)
sigma_beta_text = ax_sigma_beta.text(0, 0, r'$\sigma_{\beta}=$'+str(norm_y), color = 'lime', ha = 'center', va = 'center', fontsize = 20)
ax_sigma_beta.set_xticks([]), ax_sigma_beta.set_yticks([])
ax_sigma_beta.axis(False)
ax_tau_alpha_beta = plt.subplot2grid((24, 18), (16, 14), colspan = 4, rowspan = 3)
ax_tau_alpha_beta.set_xlim(-1, 1), ax_tau_alpha_beta.set_ylim(-1, 1)
tau_alpha_beta_text = ax_tau_alpha_beta.text(0, 0, r'$\tau_{\alpha\beta}=$'+str(shear), color = 'deepskyblue', ha = 'center', va = 'center', fontsize = 20)
ax_tau_alpha_beta.set_xticks([]), ax_tau_alpha_beta.set_yticks([])
ax_tau_alpha_beta.axis(False)

# Update function
def update(val):
    angle_rad = np.radians(rot_slider.val)
    Q = rot_matrix(angle_rad)
    sigma_n = get_sigma_n(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val, angle_rad)
    tau_nt = get_tau_nt(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val, angle_rad)
    sigma_t = get_sigma_n(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val, angle_rad+np.pi/2)
    sigma_dn = get_sigma_n(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val, angle_rad+5*np.pi/4)
    tau_dnt = get_tau_nt(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val, angle_rad+5*np.pi/4)
    
    sigma_x.set_data(x=(pos_x@Q)[0], y=(pos_x@Q)[1], dx=(sigma_n*u_sigma_x@Q)[0], dy=(sigma_n*u_sigma_x@Q)[1])
    if np.abs(sigma_n) < 0.0001:
        sigma_x.set_color('k')
    else:
        sigma_x.set_color('limegreen')
    tau_xy.set_data(x=(pos_x@Q)[0], y=(pos_x@Q)[1], dx=(tau_nt*u_tau_xy@Q)[0], dy=(tau_nt*u_tau_xy@Q)[1])
    if np.abs(tau_nt) < 0.0001:
        tau_xy.set_color('k')
    else:
        tau_xy.set_color('deepskyblue')
    traction_x.set_data(x=(pos_x@Q)[0], y=(pos_x@Q)[1], dx=((sigma_n*u_sigma_x+tau_nt*u_tau_xy)@Q)[0], dy=((sigma_n*u_sigma_x+tau_nt*u_tau_xy)@Q)[1])
    if np.abs(sigma_n) < 0.0001 and np.abs(tau_nt) < 0.0001:
        traction_x.set_color('k')
    else:
        traction_x.set_color('r')
    
    sigma_y.set_data(x=(pos_y@Q)[0], y=(pos_y@Q)[1], dx=(sigma_t*u_sigma_y@Q)[0], dy=(sigma_t*u_sigma_y@Q)[1])
    if np.abs(sigma_t) < 0.0001:
        sigma_y.set_color('k')
    else:
        sigma_y.set_color('limegreen')
    tau_yx.set_data(x=(pos_y@Q)[0], y=(pos_y@Q)[1], dx=(tau_nt*u_tau_yx@Q)[0], dy=(tau_nt*u_tau_yx@Q)[1])
    if np.abs(tau_nt) < 0.0001:
        tau_yx.set_color('k')
    else:
        tau_yx.set_color('deepskyblue')
    traction_y.set_data(x=(pos_y@Q)[0], y=(pos_y@Q)[1], dx=((sigma_t*u_sigma_y+tau_nt*u_tau_yx)@Q)[0], dy=((sigma_t*u_sigma_y+tau_nt*u_tau_yx)@Q)[1])
    if np.abs(sigma_t) < 0.0001 and np.abs(tau_nt) < 0.0001:
        traction_y.set_color('k')
    else:
        traction_y.set_color('r')
    
    sigma_d.set_data(x=(pos_d@Q)[0], y=(pos_d@Q)[1], dx=(sigma_dn*u_sigma_d@Q)[0], dy=(sigma_dn*u_sigma_d@Q)[1])
    if np.abs(sigma_dn) < 0.0001:
        sigma_d.set_color('k')
    else:
        sigma_d.set_color('limegreen')
    tau_d.set_data(x=(pos_d@Q)[0], y=(pos_d@Q)[1], dx=(tau_dnt*u_tau_d@Q)[0], dy=(tau_dnt*u_tau_d@Q)[1])
    if np.abs(tau_dnt) < 0.0001:
        tau_d.set_color('k')
    else:
        tau_d.set_color('deepskyblue')
    traction_d.set_data(x=(pos_d@Q)[0], y=(pos_d@Q)[1], dx=((sigma_dn*u_sigma_d+tau_dnt*u_tau_d)@Q)[0], dy=((sigma_dn*u_sigma_d+tau_dnt*u_tau_d)@Q)[1])
    if np.abs(sigma_dn) < 0.0001 and np.abs(tau_dnt) < 0.0001:
        traction_d.set_color('k')
    else:
        traction_d.set_color('r')                           
    wedge.set_data((wedge_points@Q)[:,0], (wedge_points@Q)[:,1])
    
    circle.set_data(mohr_circle(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val)[:,0],
                    mohr_circle(sigma_x_slider.val, sigma_y_slider.val, tau_xy_slider.val)[:,1])
    nt_line.set_data([sigma_n, sigma_t], [-tau_nt, tau_nt])
    
    alpha_mohr_pos = np.array([sigma_n, -tau_nt])
    beta_mohr_pos = np.array([sigma_t, tau_nt])
    
    if np.all(alpha_mohr_pos == beta_mohr_pos):
        alpha_mohr.set_text('')
        beta_mohr.set_text('')
        center.set_color('y')
    else:
        alpha_mohr.set_text(r'$\alpha$')
        beta_mohr.set_text(r'$\beta$')
        alpha_mohr.set_position(alpha_mohr_pos)
        beta_mohr.set_position(beta_mohr_pos)
        center.set_offsets([(sigma_n + sigma_t)/2, 0]) 
        center.set_color('k')  
    
    alpha_elmt.set_position((0.75*pos_x@Q))
    beta_elmt.set_position((0.75*pos_y@Q))
    
    sigma_alpha_text.set_text(r'$\sigma_{\alpha}=$'+str(round(sigma_n, 3)))
    sigma_beta_text.set_text(r'$\sigma_{\beta}=$'+str(round(sigma_t, 3)))
    tau_alpha_beta_text.set_text(r'$\tau_{\alpha\beta}=$'+str(round(-tau_nt, 3)))
    
    fig.canvas.draw_idle()

# Making the sliders update
rot_slider.on_changed(update)
sigma_x_slider.on_changed(update)
sigma_y_slider.on_changed(update)
tau_xy_slider.on_changed(update)
plt.show()