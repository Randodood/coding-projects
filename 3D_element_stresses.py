import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

I_HAT, J_HAT, K_HAT = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])
SLIDER_FONTSIZE = 12.5
STRESS_VECTOR_OFFSET = 1.1
ANGLE = np.linspace(0, 2*np.pi, 100)
MINIMUM_LETTER_DISTANCE = 0.05

def get_rotation_matrix(yaw, pitch, roll):
    psi = np.radians(yaw)
    theta = np.radians(pitch)
    phi = np.radians(roll)

    Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

    return np.matmul(Rz, np.matmul(Ry, Rx))

def get_stress_tensor(sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_zx, tau_yz):
    return np.array([[sigma_xx, tau_xy, tau_zx], [tau_xy, sigma_yy, tau_yz], [tau_zx, tau_yz, sigma_zz]])

def rotate_objects(stress_tensor, rotation_matrix):
    stress_tensor_prime = np.matmul(rotation_matrix, np.matmul(stress_tensor, rotation_matrix.T))
    basis = np.array([I_HAT, J_HAT, K_HAT])
    basis_prime = np.matmul(basis, rotation_matrix)
    i_prime, j_prime, k_prime = basis_prime[0], basis_prime[1], basis_prime[2]

    return stress_tensor_prime, i_prime, j_prime, k_prime

def draw_half_cube(i_hat, j_hat, k_hat):
    half_cube_p1 = i_hat + j_hat + k_hat
    half_cube_p2 = i_hat + j_hat - k_hat
    half_cube_p3 = i_hat - j_hat + k_hat
    half_cube_p4 = -i_hat + j_hat + k_hat
    half_cube_p5 = i_hat - j_hat - k_hat
    half_cube_p6 = -i_hat -j_hat + k_hat
    half_cube_p7 = -i_hat + j_hat - k_hat
    half_cube_points1 = np.array([half_cube_p1, half_cube_p2, half_cube_p5, half_cube_p3, half_cube_p1, half_cube_p4, 
                                  half_cube_p6, half_cube_p3])
    half_cube_points2 = np.array([half_cube_p2, half_cube_p7, half_cube_p4])
    half_cube1.set_data(half_cube_points1[:,0], half_cube_points1[:,1]), half_cube1.set_3d_properties(half_cube_points1[:,2])
    half_cube2.set_data(half_cube_points2[:,0], half_cube_points2[:,1]), half_cube2.set_3d_properties(half_cube_points2[:,2])
    global alpha_elmt_text, beta_elmt_text, gamma_elmt_text
    alpha_elmt_text.remove(), beta_elmt_text.remove(), gamma_elmt_text.remove()
    alpha_elmt_text = ax_elmt.text(i_hat[0], i_hat[1], i_hat[2], r'$\alpha$', ha = 'center', va = 'center', color = 'y')
    beta_elmt_text = ax_elmt.text(j_hat[0], j_hat[1], j_hat[2], r'$\beta$', ha = 'center', va = 'center', color = 'y')
    gamma_elmt_text = ax_elmt.text(k_hat[0], k_hat[1], k_hat[2], r'$\gamma$', ha = 'center', va = 'center', color = 'y')

def draw_stress(u_hat, v_hat, stress_uv, vec_stress_uv):
    endpoint = STRESS_VECTOR_OFFSET*u_hat + stress_uv*v_hat
    vec_stress_uv.set_data([STRESS_VECTOR_OFFSET*u_hat[0], endpoint[0]], [STRESS_VECTOR_OFFSET*u_hat[1], endpoint[1]])
    vec_stress_uv.set_3d_properties([STRESS_VECTOR_OFFSET*u_hat[2], endpoint[2]])
    
def draw_traction(u_hat, v_hat, w_hat, stress_u, stress_uv, stress_wu, vec_traction_u):
    endpoint = STRESS_VECTOR_OFFSET*u_hat + stress_u*u_hat + stress_uv*v_hat + stress_wu*w_hat
    vec_traction_u.set_data([STRESS_VECTOR_OFFSET*u_hat[0], endpoint[0]], [STRESS_VECTOR_OFFSET*u_hat[1], endpoint[1]])
    vec_traction_u.set_3d_properties([STRESS_VECTOR_OFFSET*u_hat[2], endpoint[2]])
    
def draw_stresses(i_hat, j_hat, k_hat, stress_tensor):
    sigma_alpha, sigma_beta, sigma_gamma = stress_tensor[0,0], stress_tensor[1,1], stress_tensor[2,2]
    tau_alpha_beta, tau_beta_gamma, tau_gamma_alpha = stress_tensor[0,1], stress_tensor[1,2], stress_tensor[2,0]
    
    draw_stress(i_hat, i_hat, sigma_alpha, vec_sigma_alpha)
    draw_stress(j_hat, j_hat, sigma_beta, vec_sigma_beta)
    draw_stress(k_hat, k_hat, sigma_gamma, vec_sigma_gamma)
    
    draw_stress(i_hat, j_hat, tau_alpha_beta, vec_tau_alpha_beta)
    draw_stress(j_hat, k_hat, tau_beta_gamma, vec_tau_beta_gamma)
    draw_stress(k_hat, i_hat, tau_gamma_alpha, vec_tau_gamma_alpha)
    
    draw_stress(j_hat, i_hat, tau_alpha_beta, vec_tau_beta_alpha)
    draw_stress(k_hat, j_hat, tau_beta_gamma, vec_tau_gamma_beta)
    draw_stress(i_hat, k_hat, tau_gamma_alpha, vec_tau_alpha_gamma)
    
    draw_traction(i_hat, j_hat, k_hat, sigma_alpha, tau_alpha_beta, tau_gamma_alpha, vec_traction_alpha)
    draw_traction(j_hat, k_hat, i_hat, sigma_beta, tau_beta_gamma, tau_alpha_beta, vec_traction_beta)
    draw_traction(k_hat, i_hat, j_hat, sigma_gamma, tau_gamma_alpha, tau_beta_gamma, vec_traction_gamma)

def draw_mohr_circle(sigma_u, sigma_v, tau_uv, mohr_circle_uv):
    C = (sigma_u + sigma_v)/2
    D = (sigma_u - sigma_v)/2
    R = np.sqrt(D**2 + tau_uv**2)
    circle_x = C + R*np.cos(ANGLE)
    circle_y = R*np.sin(ANGLE)
    mohr_circle_uv.set_data(circle_x, circle_y)
    
def draw_mohr_line(sigma_u, sigma_v, tau_uv, mohr_line_uv):
    mohr_line_uv.set_data([sigma_u, sigma_v], [-tau_uv, tau_uv])

def draw_mohr_letters(sigma_u, sigma_v, tau_uv, mohr_letter_uv_u, mohr_letter_uv_v, mohr_point_uv):
    mohr_letter_uv_u.set_position([sigma_u, -tau_uv]), mohr_letter_uv_v.set_position([sigma_v, tau_uv])
    mohr_point_uv.set_offsets([(sigma_u + sigma_v)/2, 0])
    if np.abs(sigma_u - sigma_v) <= MINIMUM_LETTER_DISTANCE and np.abs(tau_uv) <= MINIMUM_LETTER_DISTANCE/2:
        mohr_letter_uv_u.set_visible(False), mohr_letter_uv_v.set_visible(False)
        mohr_point_uv.set_visible(True)
    else:
        mohr_letter_uv_u.set_visible(True), mohr_letter_uv_v.set_visible(True)
        mohr_point_uv.set_visible(False)
        
def check_letters_overlap(sigma_u, tau_uv, tau_wu, mohr_letter_uv_u, mohr_letter_wu_u, mohr_point_u):
    mohr_point_u.set_offsets([sigma_u, (tau_uv - tau_wu)/2])
    if np.abs(tau_uv + tau_wu) <= MINIMUM_LETTER_DISTANCE:
        mohr_letter_uv_u.set_visible(False), mohr_letter_wu_u.set_visible(False)
        mohr_point_u.set_visible(True)
    else:
        mohr_letter_uv_u.set_visible(True), mohr_letter_wu_u.set_visible(True)
        mohr_point_u.set_visible(False)

def slider_update(val):
    yaw, pitch, roll = yaw_slider.val, pitch_slider.val, roll_slider.val
    sigma_x, sigma_y, sigma_z = sigma_x_slider.val, sigma_y_slider.val, sigma_z_slider.val
    tau_xy, tau_yz, tau_zx = tau_xy_slider.val, tau_yz_slider.val, tau_zx_slider.val
    R = get_rotation_matrix(yaw, pitch, roll)
    stress_tensor = get_stress_tensor(sigma_x, sigma_y, sigma_z, tau_xy, tau_zx, tau_yz)
    stress_tensor_prime, i_hat_prime, j_hat_prime, k_hat_prime = rotate_objects(stress_tensor, R)
    sigma_alpha, sigma_beta, sigma_gamma = stress_tensor_prime[0, 0], stress_tensor_prime[1, 1], stress_tensor_prime[2, 2]
    tau_alpha_beta, tau_beta_gamma, tau_gamma_alpha = stress_tensor_prime[0, 1], stress_tensor_prime[1, 2], stress_tensor_prime[2, 0]
    
    draw_half_cube(i_hat_prime, j_hat_prime, k_hat_prime)
    draw_stresses(i_hat_prime, j_hat_prime, k_hat_prime, stress_tensor_prime)
    draw_mohr_circle(sigma_x, sigma_y, tau_xy, mohr_circle_alpha_beta)
    draw_mohr_circle(sigma_y, sigma_z, tau_yz, mohr_circle_beta_gamma)
    draw_mohr_circle(sigma_z, sigma_x, tau_zx, mohr_circle_gamma_alpha)
    draw_mohr_line(sigma_alpha, sigma_beta, tau_alpha_beta, mohr_line_alpha_beta)
    draw_mohr_line(sigma_beta, sigma_gamma, tau_beta_gamma, mohr_line_beta_gamma)
    draw_mohr_line(sigma_gamma, sigma_alpha, tau_gamma_alpha, mohr_line_gamma_alpha)
    
    draw_mohr_letters(sigma_alpha, sigma_beta, tau_alpha_beta, mohr_letter_alpha_beta_alpha, mohr_letter_alpha_beta_beta, mohr_point_alpha_beta)
    draw_mohr_letters(sigma_beta, sigma_gamma, tau_beta_gamma, mohr_letter_beta_gamma_beta, mohr_letter_beta_gamma_gamma, mohr_point_beta_gamma)
    draw_mohr_letters(sigma_gamma, sigma_alpha, tau_gamma_alpha, mohr_letter_gamma_alpha_gamma, mohr_letter_gamma_alpha_alpha, mohr_point_gamma_alpha)
    
    check_letters_overlap(sigma_beta, tau_alpha_beta, tau_beta_gamma, mohr_letter_alpha_beta_beta, mohr_letter_beta_gamma_beta, mohr_point_beta)
    check_letters_overlap(sigma_gamma, tau_beta_gamma, tau_gamma_alpha, mohr_letter_beta_gamma_gamma, mohr_letter_gamma_alpha_gamma, mohr_point_gamma)
    check_letters_overlap(sigma_alpha, tau_gamma_alpha, tau_alpha_beta, mohr_letter_gamma_alpha_alpha, mohr_letter_alpha_beta_alpha, mohr_point_alpha)
    
    sigma_alpha_text.set_text(r'$\sigma_{\alpha}=$' + str(round(sigma_alpha, 2)))
    sigma_beta_text.set_text(r'$\sigma_{\beta}=$' + str(round(sigma_beta, 2)))
    sigma_gamma_text.set_text(r'$\sigma_{\gamma}=$' + str(round(sigma_gamma, 2)))
    tau_alpha_beta_text.set_text(r'$\tau_{\alpha\beta}=$' + str(round(tau_alpha_beta, 2)))
    tau_beta_gamma_text.set_text(r'$\tau_{\beta\gamma}=$' + str(round(tau_beta_gamma, 2)))
    tau_gamma_alpha_text.set_text(r'$\tau_{\gamma\alpha}=$' + str(round(tau_gamma_alpha, 2)))

# Plot initialization
plt.style.use('dark_background')
fig = plt.figure(figsize=(10,7.5))

# Half cube shaped stress element
ax_elmt = plt.subplot2grid((12,16), (0,0), colspan = 7, rowspan = 7, projection = '3d')
ax_elmt.set_title('Stress Element')
ax_elmt.set_xlim(-2, 2), ax_elmt.set_ylim(-2, 2), ax_elmt.set_zlim(-2, 2), ax_elmt.axis(False)
ax_elmt.set_box_aspect(aspect = (1, 1, 1))
ax_elmt.plot([-2, -1], [-2, -2], [-2, -2], color = 'w', linewidth = 1)
ax_elmt.plot([-2, -2], [-2, -1], [-2, -2], color = 'w', linewidth = 1)
ax_elmt.plot([-2, -2], [-2, -2], [-2, -1], color = 'w', linewidth = 1)
ax_elmt.text(-0.75, -2, -2, r'$x$', color = 'w', ha = 'center', va = 'center')
ax_elmt.text(-2, -0.75, -2, r'$y$', color = 'w', ha = 'center', va = 'center')
ax_elmt.text(-2, -2, -0.75, r'$z$', color = 'w', ha = 'center', va = 'center')

alpha_elmt_text = ax_elmt.text(None, None, None, '')
beta_elmt_text = ax_elmt.text(None, None, None, '')
gamma_elmt_text = ax_elmt.text(None, None, None, '')
half_cube1, = ax_elmt.plot([], [], [], color = 'w', alpha = 0.35, linestyle = 'dashed', linewidth = 1)
half_cube2, = ax_elmt.plot([], [], [], color = 'w', alpha = 0.35, linestyle = 'dashed', linewidth = 1)
vec_sigma_alpha, = ax_elmt.plot([], [], [], color = 'lime', linewidth = 1)
vec_sigma_beta, = ax_elmt.plot([], [], [], color = 'lime', linewidth = 1)
vec_sigma_gamma, = ax_elmt.plot([], [], [], color = 'lime', linewidth = 1)
vec_tau_alpha_beta, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_tau_beta_gamma, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_tau_gamma_alpha, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_tau_beta_alpha, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_tau_gamma_beta, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_tau_alpha_gamma, = ax_elmt.plot([], [], [], color = 'deepskyblue', linewidth = 1)
vec_traction_alpha, = ax_elmt.plot([], [], [], color = 'red', linewidth = 1)
vec_traction_beta, = ax_elmt.plot([], [], [], color = 'red', linewidth = 1)
vec_traction_gamma, = ax_elmt.plot([], [], [], color = 'red', linewidth = 1)

ax_mohr = plt.subplot2grid((12,16), (0,9), colspan = 7, rowspan = 7)
ax_mohr.set_xlim(-2, 2), ax_mohr.set_ylim(-2, 2)
ax_mohr.set_xlabel(r'$\sigma$'), ax_mohr.set_ylabel(r'$\tau$')
ax_mohr.set_xticks([-2, -1, 0, 1, 2]), ax_mohr.set_yticks([-2, -1, 0, 1, 2])
ax_mohr.set_title('Mohr Circle')
mohr_circle_alpha_beta, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linestyle = 'dashed', linewidth = 1)
mohr_circle_beta_gamma, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linestyle = 'dashed', linewidth = 1)
mohr_circle_gamma_alpha, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linestyle = 'dashed', linewidth = 1)

mohr_line_alpha_beta, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linewidth = 1)
mohr_line_beta_gamma, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linewidth = 1)
mohr_line_gamma_alpha, = ax_mohr.plot([], [], color = 'w', alpha = 0.35, linewidth = 1)

mohr_letter_alpha_beta_alpha = ax_mohr.text([], [], r'$\alpha$', color = 'y', ha = 'center', va = 'center', visible = False)
mohr_letter_alpha_beta_beta = ax_mohr.text([], [], r'$\beta$', color = 'y', ha = 'center', va = 'center', visible = False)
mohr_letter_beta_gamma_beta = ax_mohr.text([], [], r'$\beta$', color = 'y', ha = 'center', va = 'center', visible = False)
mohr_letter_beta_gamma_gamma = ax_mohr.text([], [], r'$\gamma$', color = 'y', ha = 'center', va = 'center', visible = False)
mohr_letter_gamma_alpha_gamma = ax_mohr.text([], [], r'$\gamma$', color = 'y', ha = 'center', va = 'center', visible = False)
mohr_letter_gamma_alpha_alpha = ax_mohr.text([], [], r'$\alpha$', color = 'y', ha = 'center', va = 'center', visible = False)

mohr_point_alpha_beta = ax_mohr.scatter([], [], color = 'y', s= 10, visible = False)
mohr_point_beta_gamma = ax_mohr.scatter([], [], color = 'y', s= 10, visible = False)
mohr_point_gamma_alpha = ax_mohr.scatter([], [], color = 'y', s= 10, visible = False)

mohr_point_alpha = ax_mohr.scatter([], [], color = 'y', s= 10, visible = True)
mohr_point_beta = ax_mohr.scatter([], [], color = 'y', s= 10, visible = True)
mohr_point_gamma = ax_mohr.scatter([], [], color = 'y', s= 10, visible = True)

ax_yaw = plt.subplot2grid((24, 18), (21, 0), colspan = 18, rowspan = 1)
yaw_slider = Slider(
    ax = ax_yaw,
    label = r'$Yaw$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'lightgoldenrodyellow',
)
yaw_slider.label.set_color('lightgoldenrodyellow')
yaw_slider.label.set_size(SLIDER_FONTSIZE)
yaw_slider.valtext.set_size(SLIDER_FONTSIZE)
yaw_slider.valtext.set_color('lightgoldenrodyellow')

ax_pitch = plt.subplot2grid((24, 18), (22, 0), colspan = 18, rowspan = 1)
pitch_slider = Slider(
    ax = ax_pitch,
    label = r'$Pitch$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'lightgoldenrodyellow',
)
pitch_slider.label.set_color('lightgoldenrodyellow')
pitch_slider.label.set_size(SLIDER_FONTSIZE)
pitch_slider.valtext.set_size(SLIDER_FONTSIZE)
pitch_slider.valtext.set_color('lightgoldenrodyellow')

ax_roll = plt.subplot2grid((24, 18), (23, 0), colspan = 18, rowspan = 1)
roll_slider = Slider(
    ax = ax_roll,
    label = r'$Roll$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = 'lightgoldenrodyellow',
)
roll_slider.label.set_color('lightgoldenrodyellow')
roll_slider.label.set_size(SLIDER_FONTSIZE)
roll_slider.valtext.set_size(SLIDER_FONTSIZE)
roll_slider.valtext.set_color('lightgoldenrodyellow')

ax_tau_xy = plt.subplot2grid((24,18), (20,0), colspan = 4)
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
tau_xy_slider.label.set_size(SLIDER_FONTSIZE)
tau_xy_slider.valtext.set_size(SLIDER_FONTSIZE)
tau_xy_slider.valtext.set_color('deepskyblue')

ax_tau_yz = plt.subplot2grid((24,18), (20,7), colspan = 4)
tau_yz_slider = Slider(
    ax = ax_tau_yz,
    label = r'$\tau_{yz}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'deepskyblue'
)
tau_yz_slider.label.set_color('deepskyblue')
tau_yz_slider.label.set_size(SLIDER_FONTSIZE)
tau_yz_slider.valtext.set_size(SLIDER_FONTSIZE)
tau_yz_slider.valtext.set_color('deepskyblue')

ax_tau_zx = plt.subplot2grid((24,18), (20,14), colspan = 4)
tau_zx_slider = Slider(
    ax = ax_tau_zx,
    label = r'$\tau_{zx}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'deepskyblue'
)
tau_zx_slider.label.set_color('deepskyblue')
tau_zx_slider.label.set_size(SLIDER_FONTSIZE)
tau_zx_slider.valtext.set_size(SLIDER_FONTSIZE)
tau_zx_slider.valtext.set_color('deepskyblue')

ax_sigma_x = plt.subplot2grid((24,18), (19,0), colspan = 4)
sigma_x_slider = Slider(
    ax = ax_sigma_x,
    label = r'$\sigma_{x}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'lime'
)
sigma_x_slider.label.set_color('lime')
sigma_x_slider.label.set_size(SLIDER_FONTSIZE)
sigma_x_slider.valtext.set_size(SLIDER_FONTSIZE)
sigma_x_slider.valtext.set_color('lime')

ax_sigma_y = plt.subplot2grid((24,18), (19,7), colspan = 4)
sigma_y_slider = Slider(
    ax = ax_sigma_y,
    label = r'$\sigma_{y}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'lime'
)
sigma_y_slider.label.set_color('lime')
sigma_y_slider.label.set_size(SLIDER_FONTSIZE)
sigma_y_slider.valtext.set_size(SLIDER_FONTSIZE)
sigma_y_slider.valtext.set_color('lime')

ax_sigma_z = plt.subplot2grid((24,18), (19,14), colspan = 4)
sigma_z_slider = Slider(
    ax = ax_sigma_z,
    label = r'$\sigma_{z}$',
    valmin = -1,
    valmax = 1,
    valstep = 0.05,
    valinit = 0,
    facecolor = 'lime'
)
sigma_z_slider.label.set_color('lime')
sigma_z_slider.label.set_size(SLIDER_FONTSIZE)
sigma_z_slider.valtext.set_size(SLIDER_FONTSIZE)
sigma_z_slider.valtext.set_color('lime')

ax_sigma_alpha = plt.subplot2grid((24,6), (16,0), rowspan = 2)
ax_sigma_alpha.set_xlim(-1, 1), ax_sigma_alpha.set_ylim(-1, 1), ax_sigma_alpha.set_xticks([]), ax_sigma_alpha.set_yticks([])
sigma_alpha_text = ax_sigma_alpha.text(0, 0, r'$\sigma_{\alpha}=0$', color = 'lime', ha = 'center', va = 'center', fontsize = 15)

ax_sigma_beta = plt.subplot2grid((24,6), (16,1), rowspan = 2)
ax_sigma_beta.set_xlim(-1, 1), ax_sigma_beta.set_ylim(-1, 1), ax_sigma_beta.set_xticks([]), ax_sigma_beta.set_yticks([])
sigma_beta_text = ax_sigma_beta.text(0, 0, r'$\sigma_{\beta}=0$', color = 'lime', ha = 'center', va = 'center', fontsize = 15)

ax_sigma_gamma = plt.subplot2grid((24,6), (16,2), rowspan = 2)
ax_sigma_gamma.set_xlim(-1, 1), ax_sigma_gamma.set_ylim(-1, 1), ax_sigma_gamma.set_xticks([]), ax_sigma_gamma.set_yticks([])
sigma_gamma_text = ax_sigma_gamma.text(0, 0, r'$\sigma_{\gamma}=0$', color = 'lime', ha = 'center', va = 'center', fontsize = 15)

ax_tau_alpha_beta = plt.subplot2grid((24,6), (16,3), rowspan = 2)
ax_tau_alpha_beta.set_xlim(-1, 1), ax_tau_alpha_beta.set_ylim(-1, 1), ax_tau_alpha_beta.set_xticks([]), ax_tau_alpha_beta.set_yticks([])
tau_alpha_beta_text = ax_tau_alpha_beta.text(0, 0, r'$\tau_{\alpha\beta}=0$', color = 'deepskyblue', ha = 'center', va = 'center', fontsize = 15)

ax_tau_beta_gamma = plt.subplot2grid((24,6), (16,5), rowspan = 2)
ax_tau_beta_gamma.set_xlim(-1, 1), ax_tau_beta_gamma.set_ylim(-1, 1), ax_tau_beta_gamma.set_xticks([]), ax_tau_beta_gamma.set_yticks([])
tau_beta_gamma_text = ax_tau_beta_gamma.text(0, 0, r'$\tau_{\beta\gamma}=0$', color = 'deepskyblue', ha = 'center', va = 'center', fontsize = 15)

ax_tau_gamma_alpha = plt.subplot2grid((24,6), (16,4), rowspan = 2)
ax_tau_gamma_alpha.set_xlim(-1, 1), ax_tau_gamma_alpha.set_ylim(-1, 1), ax_tau_gamma_alpha.set_xticks([]), ax_tau_gamma_alpha.set_yticks([])
tau_gamma_alpha_text = ax_tau_gamma_alpha.text(0, 0, r'$\tau_{\alpha\gamma}=0$', color = 'deepskyblue', ha = 'center', va = 'center', fontsize = 15)
    
draw_half_cube(I_HAT, J_HAT, K_HAT)
yaw_slider.on_changed(slider_update), pitch_slider.on_changed(slider_update), roll_slider.on_changed(slider_update)
sigma_x_slider.on_changed(slider_update), sigma_y_slider.on_changed(slider_update), sigma_z_slider.on_changed(slider_update)
tau_xy_slider.on_changed(slider_update), tau_yz_slider.on_changed(slider_update), tau_zx_slider.on_changed(slider_update)

plt.show()