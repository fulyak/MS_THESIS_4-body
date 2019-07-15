#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

# Import the rebound module
import rebound
import numpy as np
from numpy import arccos,sin,cos,pi,sqrt
#We simulate Kozai cycles of the two planets perturbed by a distant star 
#Objects: 1:inner planet, 2:outer planet, 3:companion star

#REBOUND uses Jacobi coordinates: orbital elements describe the particle’s orbit around 
#the centre of mass of all particles added previously

#Any components not passed automatically default to 0


M=1 #mass
m1=1*0.000954
m2=0.03
m3=1
a1=4 #semi-major axis
a2=50
a3=950
e1=0.01 #eccentricity
e2=0.01
e3=0.01
i1=65*np.pi/180. #inclination
i2=120*np.pi/180.
i3=0*np.pi/180.

# Create Simulation object
sim = rebound.Simulation()

#Units
sim.units = ('yr', 'AU', 'Msun')

# Add particle to rebound
sim.add(m = M, hash = "primary star")
sim.add(m = m1, a =a1 , e =e1, inc = i1, hash = "planet")
sim.add(m = m2, a = a2, e = e2, inc=i2, hash = "planet 2")
sim.add(m = m3, a = a3, e = e3, inc=i3, hash = "binary star")
ps=sim.particles

sim.integrator = "whfast" 
sim.ri_whfast.corrector = 11
sim.ri_whfast.safe_mode = 0
sim.dt=ps[1].P/500 #time step
tmax = 10*10**6 #total integration time
N=1000 #number of data

times = np.linspace(0.,tmax,N)

E0=sim.calculate_energy() #initial total energy

#integrate
for t in times:	
	sim.integrate(t,exact_finish_time = 0) 
	Irel=arccos(cos(ps[1].inc)*cos(ps[2].inc)+sin(ps[1].inc)*sin(ps[2].inc)*cos(ps[1].Omega-ps[2].Omega))*180./np.pi
	print(t,ps[1].e,ps[1].inc*180./np.pi,ps[1].omega*180./np.pi,ps[1].Omega*180./np.pi,ps[2].e,180./np.pi*ps[2].inc,180./np.pi*ps[2].omega,180./np.pi*ps[2].Omega,abs((sim.calculate_energy()-E0)/E0),Irel,ps[3].e,180./np.pi*ps[3].inc,180./np.pi*ps[3].omega,180./np.pi*ps[3].Omega)
 
