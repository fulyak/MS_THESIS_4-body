#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import rebound
import numpy as np
from numpy import cos,sin

from input import get_orbitpar1
from input import get_orbitpar2
from input import get_orbitpar3
from input import get_integrate

#This code is used to scan various regions of the parameter space
#Initial parameters are taken from an external file named "input"

totaltime,numberofdata,integrator,init,final,step=get_integrate()
variable_domain = np.linspace(init,final,step)

M_J=0.0009542 #Jupiter mass

for variable in variable_domain:

	m1,a1 ,e1,i1,omega1,Omega1=get_orbitpar1(variable)
	m2,a2 ,e2,i2,omega2,Omega2=get_orbitpar2(variable)
	m3,a3,e3,dt=get_orbitpar3(variable)

	sim = rebound.Simulation()
	sim.units = ('yr', 'AU', 'Msun')
	ps=sim.particles	
	sim.add(m = 1.0, hash = "primary star")
	sim.add(m = m1*M_J, a = a1, e = e1, inc = i1*np.pi/180.,omega=omega1*np.pi/180,Omega=Omega1*np.pi/180, hash = "planet1")
	sim.add(m = m2*M_J, a = a2, e = e2, inc = i2*np.pi/180.,omega=omega2*np.pi/180,Omega=Omega2*np.pi/180, hash = "planet2")
	sim.add(m = m3, a = a3, e = e3,  hash = "binary star")
	
	sim.integrator =integrator
	sim.dt = ps[1].P/dt #time step
	sim.ri_whfast.safe_mode = 0
	sim.ri_whfast.corrector = 11
	N = numberofdata #store the orbital elements at N times during the interval
	E0=sim.calculate_energy() #total energy of the system
	tmax = totaltime #total integration time 
	times = np.linspace(0.,tmax,N)

	for t in times:
		sim.integrate(t,exact_finish_time = 0) 
		Irel=np.arccos((cos(ps[1].inc)*cos(ps[2].inc)+sin(ps[1].inc)*sin(ps[2].inc)*cos(ps[1].Omega-ps[2].Omega)))*180./np.pi #relative inclination
		if variable==init:
			if t==0:
				print('#m1=',m1,' ','a1=',a1,' ','e1=',e1,' ','omega1=',omega1,' ','Omega1=',Omega1)
				print('#m2=',m2,' ','a2=',a2,' ','e2=',e2,' ','omega2=',omega2,' ','Omega1=',Omega2)
				print('#m3=',m3,' ','a3=',a3,' ','e3=',e3)			
				print('#dt=',dt,' ','tmax=',totaltime,' ','N=',numberofdata,' ','init=',init,' ','final=',final)
				
		print(t,ps[1].e,ps[1].inc*180./np.pi,ps[1].omega*180./np.pi,ps[1].Omega*180./np.pi,ps[2].e,180./np.pi*ps[2].inc,180./np.pi*ps[2].omega,180./np.pi*ps[2].Omega,abs((sim.calculate_energy()-E0)/E0),Irel)






