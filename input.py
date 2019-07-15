
def get_orbitpar1(variable):

	m1=1
	a1=4
	e1=0.01
	i1=65 
	omega1=0 
	Omega1=0 
	return m1,a1,e1,i1,omega1,Omega1

def get_orbitpar2(variable):
	m2=0.03/0.0009542
	a2=50
	e2=variable 
	i2=120 
	omega2=0 
	Omega2=0
	return m2,a2,e2,i2,omega2,Omega2

def get_orbitpar3(variable):

	m3=1 
	a3=950 
	e3=0.01 
	dt=5000
	return m3,a3,e3,dt

def get_integrate():
 
	totaltime=5*10**6 
	numberofdata=1000
	integrator="whfast"
	init=0.01
	final=0.8
	step=8
	return totaltime,numberofdata,integrator,init,final,step


