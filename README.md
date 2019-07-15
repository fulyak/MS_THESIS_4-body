# MS_THESIS_4-body
This repo contains my Python codes for numerical simulations of the 4-body systems.

1)The REBOUND code:
This code simulates a generic 4-body system by using REBOUND package [1]. The number of bodies can be increased as desired. The simulation employs the WHFast integrator. The numerical analysis deals with exact Newton's equations of motion in the Jacobi coordinates. To minimize the error in the energy, one has to choose a time scale smaller than the orbital period of the bodies. As an output, one can obtain the position, velocity, orbital parameters, total energy and angular momentum, etc. of the bodies.

2)The REBOUND code with loop:
In addition to the previous code, this one can loop over the orbital elements, masses and time parameters to scan various parts of the parameter space. The initial data should be provided externally by creating a type of batch file named "input". The results are appended to one another for each loop in a single output file.

3)Gauss-TPO code:
I prepared this code to simulate the secular evolution of 4-body systems for my thesis. The system  of our interest is made up of a cental star, two planets and a companion star. In this code, we use the combination of the Gaussian ring algorithm [2] and the test particle octupole approximation. For the interaction between the two planets, which can be closely separated, we use the Gaussian ring algorhtim with a softening parameter b. For the effect of the companion star, we use the Hamiltionan perturbation theory, where the ratio of the semi-major axes is expanded up to octupole order in the disturbing function of the companion star. 
Also, since the mass of each planet is much smaller than that of the companion star, the test particle approximation is used, i.e., the planets do not affect the companion star and its orbital plane is approximated to be fixed. The double averaged equations of motion are integrated using Bulirsch-Stoer integrator [3] with much larger scales than orbital periods. In the Gaussian ring algorithm, first a single average is performed over the orbit of a perturbing body analytically. Then, another average is taken over the orbit of the perturbed body using the discrete Fourier transform. In this numerical analysis, the orbit of the perturbed body is divided into N equally spaced points in eccentric anamoly. One has to choose a suitable value of N to minimize the energy error. Also, the BS integrator introduces its own tolerance. As an output, one can obtain the orbital parameters and angular momenta of the bodies.

References:
[1] see the main webpage of REBOUND package: https://rebound.readthedocs.io/en/latest/
[2] arXiv:0811.2812
[3] Computational Physics by Mark E. J. Newman, Createspace Independent Pub, 2012.
