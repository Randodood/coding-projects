import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

SLIDER_FONTSIZE = 12
SLIDER_COLOR = 'steelblue'
CUBE_COLOR = 'skyblue'
EDGE_WIDTH = 0.75

# Functions to initiate the cube elements
def initiate_cube_vertices(dimension, edge_length = 2):
    vertices = []
    for i in range(dimension):
        vertices.append(generate_vertice_column(2**(dimension - (i+1)), 2**(i+1)))
    vertices = edge_length*np.array(vertices).T/2
    return vertices

def generate_vertice_column(repetition, times):
    col = []
    group = np.ones(repetition)
    for i in range(times):
        col.append(group*(-1)**i)
    col = np.array(col)
    return col.flatten()

def initiate_cube_edges():
    edges = []
    for i in range(32):
        edge, = ax_cube.plot([], [], [], color = CUBE_COLOR, linewidth = EDGE_WIDTH)
        edges.append(edge)
    return edges

# Functions to draw the cube in the 3D space
def draw_edge(edge, vertices, row_index1, row_index2):
    edge.set_data([vertices[row_index1,0], vertices[row_index2,0]], [vertices[row_index1,1], vertices[row_index2,1]])
    edge.set_3d_properties([vertices[row_index1,2], vertices[row_index2,2]])

def draw_cube(edges, vertices):
    draw_edge(edges[0], vertices, 0, 1), draw_edge(edges[1], vertices, 0, 2), draw_edge(edges[2], vertices, 0, 4)
    draw_edge(edges[3], vertices, 0, 8), draw_edge(edges[4], vertices, 1, 3), draw_edge(edges[5], vertices, 1, 5)
    draw_edge(edges[6], vertices, 1, 9), draw_edge(edges[7], vertices, 2, 3), draw_edge(edges[8], vertices, 2, 6)
    draw_edge(edges[9], vertices, 2, 10), draw_edge(edges[10], vertices, 3, 7), draw_edge(edges[11], vertices, 3, 11)
    draw_edge(edges[12], vertices, 4, 5), draw_edge(edges[13], vertices, 4, 12), draw_edge(edges[14], vertices, 4, 6)
    draw_edge(edges[15], vertices, 5, 7), draw_edge(edges[16], vertices, 5, 13), draw_edge(edges[17], vertices, 6, 7)
    draw_edge(edges[18], vertices, 6, 14), draw_edge(edges[19], vertices, 7, 15), draw_edge(edges[20], vertices, 8, 9)
    draw_edge(edges[21], vertices, 8, 10), draw_edge(edges[22], vertices, 8, 12), draw_edge(edges[23], vertices, 9, 11)
    draw_edge(edges[24], vertices, 9, 13), draw_edge(edges[25], vertices, 10, 11), draw_edge(edges[26], vertices, 10, 14)
    draw_edge(edges[27], vertices, 11, 15), draw_edge(edges[28], vertices, 12, 13), draw_edge(edges[29], vertices, 12, 14)
    draw_edge(edges[30], vertices, 13, 15), draw_edge(edges[31], vertices, 14, 15)

# Function to get the rotation matrix
def get_rotation_matrix(xy, xz, yz, xw, yw, zw):
    R_xy = np.array([[np.cos(xy), -np.sin(xy), 0, 0], [np.sin(xy), np.cos(xy), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    R_yz = np.array([[1, 0, 0, 0], [0, np.cos(yz), -np.sin(yz), 0], [0, np.sin(yz), np.cos(yz), 0], [0, 0, 0, 1]])
    R_xz = np.array([[np.cos(xz), 0, -np.sin(xz), 0], [0, 1, 0, 0], [np.sin(xz), 0, np.cos(xz), 0], [0, 0, 0, 1]])
    R_xw = np.array([[np.cos(xw), 0, 0, -np.sin(xw)], [0, 1, 0, 0], [0, 0, 1, 0], [np.sin(xw), 0, 0, np.cos(xw)]])
    R_yw = np.array([[1, 0, 0, 0], [0, np.cos(yw), 0, -np.sin(yw)], [0, 0, 1, 0], [0, np.sin(yw), 0, np.cos(yw)]])
    R_zw = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, np.cos(zw), -np.sin(zw)], [0, 0, np.sin(zw), np.cos(zw)]])

    return np.matmul(R_zw, np.matmul(R_yw, np.matmul(R_xw, np.matmul(R_xz, np.matmul(R_yz, R_xy)))))

# Function called by all sliders
def rotate_cube(val):
    xy, xz, yz = np.radians(xy_slider.val), np.radians(xz_slider.val), np.radians(yz_slider.val)
    xw, yw, zw = np.radians(xw_slider.val), np.radians(yw_slider.val), np.radians(zw_slider.val)
    R = get_rotation_matrix(xy, xz, yz, xw, yw, zw)
    cube_vertices = np.matmul(initial_cube_vertices, R)
    draw_cube(cube_edges, cube_vertices)

# Initialization of plot
plt.style.use('dark_background')
fig = plt.figure(figsize=(10,7.5))

# 3D space to plot cube
ax_cube = plt.subplot2grid((15,10), (0,0), colspan = 10, rowspan = 13, aspect = 'equal', projection = '3d')
ax_cube.set_xlim(-1.5, 1.5), ax_cube.set_ylim(-1.5, 1.5), ax_cube.set_zlim(-1.5, 1.5), ax_cube.axis('off')
ax_cube.set_title('Projected Tesseract', color = SLIDER_COLOR, fontsize = 20, y = 1.1)

# xyz vectors
ax_cube.plot([0, 0.5], [0, 0], [0, 0], color = 'w', linewidth = 1)
ax_cube.plot([0, 0], [0, 0.5], [0, 0], color = 'w', linewidth = 1)
ax_cube.plot([0, 0], [0, 0], [0, 0.5], color = 'w', linewidth = 1)
ax_cube.text(0.6, 0, 0, r'$x$', color = 'w', ha = 'center', va = 'center')
ax_cube.text(0, 0.6, 0, r'$y$', color = 'w', ha = 'center', va = 'center')
ax_cube.text(0, 0, 0.6, r'$z$', color = 'w', ha = 'center', va = 'center')

# Sliders for angles
ax_xy = plt.subplot2grid((24, 18), (21, 0), colspan = 8)
xy_slider = Slider(
    ax = ax_xy,
    label = r'$\theta_{xy}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
xy_slider.label.set_color(SLIDER_COLOR)
xy_slider.label.set_size(SLIDER_FONTSIZE)
xy_slider.valtext.set_size(SLIDER_FONTSIZE)
xy_slider.valtext.set_color(SLIDER_COLOR)

ax_yz = plt.subplot2grid((24, 18), (22, 0), colspan = 8)
yz_slider = Slider(
    ax = ax_yz,
    label = r'$\theta_{yz}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
yz_slider.label.set_color(SLIDER_COLOR)
yz_slider.label.set_size(SLIDER_FONTSIZE)
yz_slider.valtext.set_size(SLIDER_FONTSIZE)
yz_slider.valtext.set_color(SLIDER_COLOR)

ax_xz = plt.subplot2grid((24, 18), (23, 0), colspan = 8)
xz_slider = Slider(
    ax = ax_xz,
    label = r'$\theta_{xz}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
xz_slider.label.set_color(SLIDER_COLOR)
xz_slider.label.set_size(SLIDER_FONTSIZE)
xz_slider.valtext.set_size(SLIDER_FONTSIZE)
xz_slider.valtext.set_color(SLIDER_COLOR)

ax_xw = plt.subplot2grid((24, 18), (21, 10), colspan = 8)
xw_slider = Slider(
    ax = ax_xw,
    label = r'$\theta_{xw}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
xw_slider.label.set_color(SLIDER_COLOR)
xw_slider.label.set_size(SLIDER_FONTSIZE)
xw_slider.valtext.set_size(SLIDER_FONTSIZE)
xw_slider.valtext.set_color(SLIDER_COLOR)

ax_yw = plt.subplot2grid((24, 18), (22, 10), colspan = 8)
yw_slider = Slider(
    ax = ax_yw,
    label = r'$\theta_{yw}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
yw_slider.label.set_color(SLIDER_COLOR)
yw_slider.label.set_size(SLIDER_FONTSIZE)
yw_slider.valtext.set_size(SLIDER_FONTSIZE)
yw_slider.valtext.set_color(SLIDER_COLOR)

ax_zw = plt.subplot2grid((24, 18), (23, 10), colspan = 8)
zw_slider = Slider(
    ax = ax_zw,
    label = r'$\theta_{zw}$',
    valmin = 0,
    valmax = 360,
    valstep = 0.5,
    valinit = 0,
    facecolor = SLIDER_COLOR,
)
zw_slider.label.set_color(SLIDER_COLOR)
zw_slider.label.set_size(SLIDER_FONTSIZE)
zw_slider.valtext.set_size(SLIDER_FONTSIZE)
zw_slider.valtext.set_color(SLIDER_COLOR)

# Initialize cube
initial_cube_vertices = initiate_cube_vertices(4)
cube_edges = initiate_cube_edges()
draw_cube(cube_edges, initial_cube_vertices)

# Update cube when sliders are interacted with
xy_slider.on_changed(rotate_cube), xz_slider.on_changed(rotate_cube), yz_slider.on_changed(rotate_cube)
xw_slider.on_changed(rotate_cube), yw_slider.on_changed(rotate_cube), zw_slider.on_changed(rotate_cube)

plt.show()