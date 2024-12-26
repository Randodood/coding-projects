All files are independent of each other and can fully run without requiring an editor.

1D_Schrodinger_simulator.py

Simulates the behaviour of an initially gaussian wavicle. The wavefunction is replaced by a vector spanning the space with discrete steps. The Hamiltonian opertator is replaced by a matrix which includes a discrete second derivative. Eigenfunctions are then the eigenvectors of the Hamiltonian matrix. Built in potentials are: harmonic oscillator, potential barrier, and free particle.

![harmonic_example](https://github.com/user-attachments/assets/735e4db1-00dc-42c6-9c80-896edf571ec3)

2D_Schrodinger_simulator.py

Like with the 1D case, simulates the behaviour of an initially gaussian wavicle. Extra graphs are shown to highlight the real and imaginary behaviour of the wavefunction as well as the potential function. Built in potentials are all the ones included in the 1D case plus the circular well and a potential simulating the double slit experiment.

![circular_example](https://github.com/user-attachments/assets/5c0f9ebf-714d-47cf-b679-5e7537546148)

draw_chaotic_attractor.py

Generates points with a random position and models their time evolution in a known chaotic system. Uses the 4th order Runge-Kutta algorithm to compute the next point in the coupled set of differential equations. Built in attractors are: Lorenz, Dadras, Finance, Nosé-Hoover, Three-Scroll United, Wang-Sun, Lorenz83, Chen, and Thomas. Most attractors and default parameters are taken from: https://www.dynamicmath.xyz/strange-attractors/

![lorenz_example](https://github.com/user-attachments/assets/95de40df-4473-44d0-932a-7b2f8003c6ee)

particle_in_a_sphere_eigenfunctions.py

Graphs the nodal surfaces of the eigenfunction with parameters n, l and m. The eigenfunctions are the solutions of the Laplace equation in spherical coordinates, precicely in a sphere of radius a. Even though these aren't technically orbitals, this file uses the Condon-Shortley phase convention.

![n3l2m1](https://github.com/user-attachments/assets/c682d462-d670-4aba-a232-305353cf3f4c)

hydrogen_orbitals.py

Generates points with random positions based on the electron probability functions of the hydrogen-like orbitals. The probablitity distribution for parameters n, l, and m, corresponds to the square magnitude of the eigenfunctions of the time-independant Schrödinger equation in spherical coordinates with a Coulomb potential function generated the atom nucleus. Given that the solutions are obtained through separation of variables, the probabilities spanning each variable can also be separated. This code computes the coordinates for each points separately based on one probablity function per coordinate, then puts them together in 3D space.

![n4l1m0_h](https://github.com/user-attachments/assets/2dd9a870-71ff-4144-a799-08cb1fed4374)

plane_stresses.py

Draws the normal stress, shear stress, and traction vectors onto all sides of a rotated wedge stress element based on xy stress input. The Mohr circle is shown on a second graph. This file uses the counter-clockwise positive convention for shear stress.

![plane_stresses_example2](https://github.com/user-attachments/assets/6ac91864-0631-4bd2-a87c-5fbf12af4cd2)
