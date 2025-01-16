All files are independent of each other and can fully run without requiring an editor.

1D_Schrodinger_simulator.py

Simulates the behaviour of an initially gaussian wavicle. The wavefunction is replaced by a vector spanning the space with discrete steps. The Hamiltonian opertator is replaced by a matrix which includes a discrete second derivative. Eigenfunctions are then the eigenvectors of the Hamiltonian matrix. Built in potentials are: harmonic oscillator, potential barrier, and free particle.

![harmonic example](https://github.com/user-attachments/assets/735e4db1-00dc-42c6-9c80-896edf571ec3)

2D_Schrodinger_simulator.py

Like with the 1D case, simulates the behaviour of an initially gaussian wavicle. Extra graphs are shown to highlight the real and imaginary behaviour of the wavefunction as well as the potential function. Built in potentials are all the ones included in the 1D case plus the circular well and a potential simulating the double slit experiment.

![circular example](https://github.com/user-attachments/assets/5c0f9ebf-714d-47cf-b679-5e7537546148)

draw_chaotic_attractor.py

Generates points with a random position and models their time evolution in a known chaotic system. Uses the 4th order Runge-Kutta algorithm to compute the next point in the coupled set of differential equations. Built in attractors are: Lorenz, Dadras, Finance, Nosé-Hoover, Three-Scroll United, Wang-Sun, Lorenz83, Chen, and Thomas. Most attractors and default parameters are taken from: https://www.dynamicmath.xyz/strange-attractors/

![lorenz example](https://github.com/user-attachments/assets/95de40df-4473-44d0-932a-7b2f8003c6ee)

particle_in_a_sphere_eigenfunctions.py

Graphs the nodal surfaces of the eigenfunction with parameters n, l and m. The eigenfunctions are the solutions of the Laplace equation in spherical coordinates, precicely in a sphere of radius a. Even though these aren't technically orbitals, this file uses the Condon-Shortley phase convention.

![particle_in_a_sphere_eigenfuctions example; n3l2m1](https://github.com/user-attachments/assets/c682d462-d670-4aba-a232-305353cf3f4c)

hydrogen_orbitals.py

Generates points with random positions based on the electron probability functions of the hydrogen-like orbitals. The probablitity distribution for parameters n, l, and m, corresponds to the square magnitude of the eigenfunctions of the time-independant Schrödinger equation in spherical coordinates with a Coulomb potential function generated the atom nucleus. Given that the solutions are obtained through separation of variables, the probabilities spanning each variable can also be separated. This code computes the coordinates for each points separately based on one probablity function per coordinate, then puts them together in 3D space.

![hydrogen_orbitals example; n4l1m0_h](https://github.com/user-attachments/assets/2dd9a870-71ff-4144-a799-08cb1fed4374)

plane_stresses.py

Draws the normal stress, shear stress, and traction vectors onto all sides of a rotated wedge stress element based on xy stress input. The Mohr circle is shown on a second graph and uses the following sign convention: the shear stress axis points upwards and clockwise shear stress is considered positive. This convention is listed as convention #3 in the Wikipedia article on Mohr's circle: https://en.wikipedia.org/wiki/Mohr%27s_circle#Mohr-circle-space_sign_convention

![plane_stresses example](https://github.com/user-attachments/assets/a288126d-8792-4408-8a82-0ac4c13b44cc)

draw_geodesics.py

Draws a few geodesic curves, with set initial positions and velocities, in both extrinsic and intrinsic representations. This file only considers 3D manifolds that happen to be functions of x and y. Included functions with adjustable parameters are plane (linear), quadratic, gaussian, sinusoidal, and radial (spherical or hyperbolic).

![geodesics example](https://github.com/user-attachments/assets/af6a7954-c554-4467-84e5-753c5267c2d6)

3D_element_stresses.py

Draws the normal and the shear stresses as well as the traction vector on 3 faces of a cube. Faces can be rotated using Euler angles. This file uses a rotation matrix with the transformation order being yaw first, then pitch, and roll. The exact rotation matrix can be found on Wikipedia's article about rotation matrices as the second example under General 3D rotations: https://en.wikipedia.org/wiki/Rotation_matrix. The Mohr's circles are also drawn on a second graph.

beam_deflection.py

Graphs the the deflection and slope of a beam under a distributed load, a force, and a moment. Displays the reaction forces and moments located at the ends of the beam. All 4 possible support conditions, which are fixed-free, fixed-pinned, pinned-pinned, and fixed-fixed, are available.

![beam_deflection_example](https://github.com/user-attachments/assets/c5f75ae7-bff6-4431-862e-3c30506ec95a)