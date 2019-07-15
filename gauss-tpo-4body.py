#!/usr/bin/python -u 
# -*- coding: utf-8 -*-
from __future__ import print_function
from pylab import plot,show,xlim,ylim,xlabel
from numpy import array,zeros,abs,sqrt,pi,cos,sin,arccos,multiply,dot,cross,arange,empty
from numpy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy

#We simulate Kozai cycles of the two planets perturbed by a distant star using the Gauss's method.
#Primary star(M), Binary star(a,e,m), Planet1(a1,e1,E1,m1,omega1,Omega1,i1),Planet2(a2,e2,E2,m2,omega2,Omega2,i2)
#Units: (yr, AU, Msun)


###Calculation of the planetary interactions:

#Perform analytical average over perturbing orbit and numerical average (discrete Fourier transfrom) over perturbed orbit
#by using the work of Tremaine,Touma,Kazandjian: "Gauss’s Method for Secular Dynamics, Softened".

###Calculation of the interactions between the companion star and each planet:
#Use the test particle octupole approximation



#Since the orbit of the binary star carries the majority of the angular momentum of the system,
#it is taken as the reference plane (I = 0) with coordinates (x,y,z).



G = 4.* pi**2 #gravitaional constant

#masses
m1 =1*0.0009542 
m2 =0.03
m=1
M = 1

#initial orbital parameters

a1 = 4 #semimajor axis
a2=50 
a = 950
e = 0.01 #eccentricity
e10= 0.01
e20= 0.01

i10 = 65. *pi/180. #inclination relative to orbit of the companion star 
i20 = 120. *pi/180.
omega10 = 0. *pi/180. #argument of periapsis
omega20 = 0. *pi/180. 
Omega10 = 0* pi/180. #longitude of the ascending node
Omega20 = 0* pi/180.

n1 = sqrt(G *(M + m1)/a1**3) #mean motion
n2 = sqrt(G *(M + m2)/a2**3)
b = 0.01*a1 #softening parameter

#fixed unit vectors
xhat = [1, 0, 0]
yhat = [0, 1, 0]
zhat = [0, 0, 1]

#periods
P1 = sqrt(a1**3/M)
P2 = sqrt(a2**3/M)

#discrete FT parameters
K1 = 300
N1=K1+1 #number of equally spaced points in eccentric anomaly
K2 = 200
N2=K2+1

totaltime = 3*10**6 #number of big steps
H=5*10**3 #size of big steps 
numofsteps=totaltime/H

delta=10**(-5) #tolerance of BS integator
error=2*delta #integration error

#In the calculation of elliptic integrals, we use Chebyshev polynomial expansion method.
#Table is taken from "Chebyshev Polynomial Expansions of Complete Elliptic Integrals", W. J. Cody.

coa=zeros((30,1),float)
cop=zeros((28,1),float)

#Table Ib ,page:252
coa[0] = 1.081544554559937186096290 *pi/2
coa[1] = .044573221432776369067285 *pi
coa[2] = .004245760369819504775625 *pi
coa[3] = .000502612797966246046695 *pi
coa[4] = .000065770168092913324847 *pi
coa[5] = .000009117111066723701032 *pi
coa[6] = .000001312028529098578893 *pi
coa[7] = .000000193836613347125696 *pi
coa[8] = .000000029199279326654288 *pi
coa[9] = .000000004464889073184542 *pi
coa[10] = .000000000690892119069053* pi
coa[11] = .000000000107945779201563 *pi
coa[12] = .000000000017001189289808 *pi
coa[13] = .000000000002695771410000* pi
coa[14] = .000000000000429921254075 *pi
coa[15] = .000000000000068905869287 *pi
coa[16] = .000000000000011091989146 *pi
coa[17] = .000000000000001792348951 *pi
coa[18] = .000000000000000290608238 *pi
coa[19] = .000000000000000047261548 *pi
coa[20] = .000000000000000007707066 *pi
coa[21] = .000000000000000001259903 *pi
coa[22] = .000000000000000000206421 *pi
coa[23] = .000000000000000000033889 *pi
coa[24] = .000000000000000000005574 *pi
coa[25] = .000000000000000000000918 *pi
coa[26] = .000000000000000000000152 *pi
coa[27] = .000000000000000000000025 *pi
coa[28] = .000000000000000000000004 *pi
coa[29] = .000000000000000000000001 *pi

#Table IIa ,page:255
cop[0] = 2.928225850405146882999545/2
cop[1] = -.109838557243451911762083
cop[2] = -.003370779633972361482362
cop[3] = -.000235300858731369414039
cop[4] = -.000021764144792006684306
cop[5] = -.000002330164928439468235
cop[6] = -.000000272992738839219275
cop[7] = -.000000033998892039790023
cop[8] = -.000000004425755444003036
cop[9] = -.000000000595739318488316
cop[10] = -.000000000082323461496100
cop[11] = -.000000000011618876971255
cop[12] = -.000000000001668577566166
cop[13] = -.000000000000243130062812
cop[14] = -.000000000000035866435645
cop[15] = -.000000000000005347403815
cop[16] = -.000000000000000804634863
cop[17] = -.000000000000000122057159
cop[18] = -.000000000000000018647874
cop[19] = -.000000000000000002867194
cop[20] = -.000000000000000000443363
cop[21] = -.000000000000000000068912
cop[22] = -.000000000000000000010761
cop[23] = -.000000000000000000001687
cop[24] = -.000000000000000000000266
cop[25] = -.000000000000000000000042
cop[26] = -.000000000000000000000007
cop[27] = -.000000000000000000000001


#define the complete elliptic integrals
def EllK(k):
	sum = 0
	N=29
	for nn in range(N+1):
		Tnstar= cos(nn* arccos(4 *k**2 - 1))
		sum+= coa[nn]*Tnstar
	return sum

def EllE(k):
	sum = 0
	N=27
	for nn in range(N+1):
		Tnstar= cos(nn *arccos(4 *k**2 - 1))
		sum+= cop[nn]*Tnstar
	return sum

#initiate vectors for each body

x1hat=zeros(3,float)
y1hat=zeros(3,float)
z1hat=zeros(3,float)
r1hat=zeros(3,float)
t1hat=zeros(3,float)
acceleration1=zeros(3,float)

x2hat=zeros(3,float)
y2hat=zeros(3,float)
z2hat=zeros(3,float)
r2hat=zeros(3,float)
t2hat=zeros(3,float)
acceleration2=zeros(3,float)

acceleration3=zeros(3,float)
rphat=zeros(3,float)
tphat=zeros(3,float)

#notation:
#xp-> x', p:primed  is used for the perturbed body in the 3 body problem
#Here,in the 4 body problem, we use
#for m1: xp-> x1 and for m2: xp-> x2

#def3 function calculates the rates of change of orbital parameters of m1 and m2 due to the 
#effect of the companion star using the double averaged interaction potential (quadrupole and octupole level) 

#if we want to obtain the effect of the companion star on m1:f3(i,omega,Omega,ep,ap -> i1,omega1,Omega1,e1,a1) 
#if we want to obtain the effect of the companion star on m2:f3(i,omega,Omega,ep,ap -> i2,omega2,Omega2,e2,a2) 
def f3(i,omega,Omega,ep,ap): 

	#transform the orbit of m1 or m2 from the fixed coordinate system to the orbital (perifocal) coordinate system
	xi=cos(Omega)* cos(omega) - sin(Omega) *cos(i)* sin(omega)
	xj=sin(Omega)* cos(omega) + cos(Omega) *cos(i) *sin(omega)      	
	xk=sin(i)*sin(omega)
	yi=-cos(Omega) *sin(omega) - sin(Omega)* cos(i) *cos(omega)
	yj=-sin(Omega) *sin(omega) + cos(Omega) *cos(i) *cos(omega)
	yk=sin(i)* cos(omega)
	zi=sin(i)* sin(Omega)
	zj=-sin(i) *cos(Omega)
	zk=cos(i)

	xphat= [xi,xj,xk] #unit vectors e.g x1^=(xi x^+xj y^+xk z^)	
	yphat= [yi,yj,yk] 
	zphat= [zi,zj,zk] 

	#potential terms
	epsoct=(ap/a)*e/(1-e**2) 

	fay0=(G*m*ap**2)/(a**3*(1-e**2)**1.5)

	j=sqrt(1-ep**2) #angular momentum of the perturbed body
	#components of the angular momentum vector in the fixed coordinate system
	jx=j*zi 
	jy=j*zj
	jz=j*zk
	jvec=array([jx,jy,jz],float)
	
	#components of the eccenticity vector in the fixed coordinate system
	ex=ep*xi
	ey=ep*xj 
	ez=ep*xk
	evec=array([ex,ey,ez],float) #eccentricity vector of the pertubed body


	#The expressions below are obtained using Mathematica from the octupole order potential
	# in Eqs. (3.57-3.58) in my thesis
	delj=array([-(75./32)* epsoct* ez* fay0* jz, 0,fay0* ((3* jz)/4. + (75/64.) *epsoct* (-2* ez* jx - 2* ex* jz))],float)
	dele=array([fay0* ((3* ex)/2. + (75./64)* epsoct* (1/5 - (16* ex**2)/5. + 7* ez**2 - (8/5.)* (ex**2 + ey**2 + ez**2) - jz**2)), ((3* ey)/2. - (15* epsoct* ex* ey)/4.)* fay0, fay0* (-((9* ez)/4.) + (75/64.)* epsoct* ((54* ex* ez)/5. - 2* jx* jz))],float)
	
	#The rates of change of the angular momentum and eccentricity vectors:
	jdot = (cross(jvec, delj) + cross(evec, dele))/sqrt(G* M*ap)
	edot = (cross(jvec, dele) + cross(evec, delj))/sqrt(G* M*ap)
	
	#rate of change of orbital parameters of the perturbed body due to effect of the companion star
	fi3=(-sin(omega)*dot(xphat,jdot)-cos(omega)*dot(yphat,jdot))/norm(jvec)
	fOmega3=(cos(omega)*dot(xphat,jdot)-sin(omega)*dot(yphat,jdot))/(norm(jvec)*sin(i))
	fomega3=dot(edot,yphat)/ep-cos(i)*fOmega3
	fe3=dot(edot,xphat)
	

	
	return array([fi3,fOmega3,fomega3,fe3],float)


#Calculate the effect of m2 on m1
#First, perform analytical average over perturbing orbit, then numerical average(discrete FT) over perturbed orbit

#accelation2 function calculates the averaged acceleration of m1 analytically due the effect of m2 
def acceleration2(i1,omega1,Omega1,e1,E1,a1,i2,omega2,Omega2,e2,a2): 
	
	x1i=cos(Omega1)* cos(omega1) - sin(Omega1) *cos(i1)* sin(omega1)
	x1j=sin(Omega1)* cos(omega1) + cos(Omega1) *cos(i1) *sin(omega1)      	
	x1k=sin(i1)*sin(omega1)
	y1i=-cos(Omega1) *sin(omega1) - sin(Omega1)* cos(i1) *cos(omega1)
	y1j=-sin(Omega1) *sin(omega1) + cos(Omega1) *cos(i1) *cos(omega1)
	y1k=sin(i1)* cos(omega1)
	z1i=sin(i1)* sin(Omega1)
	z1j=-sin(i1) *cos(Omega1)
	z1k=cos(i1)
	
	x2i=cos(Omega2)* cos(omega2) - sin(Omega2) *cos(i2)* sin(omega2)
	x2j=sin(Omega2)* cos(omega2) + cos(Omega2) *cos(i2) *sin(omega2)      	
	x2k=sin(i2)*sin(omega2)
	y2i=-cos(Omega2) *sin(omega2) - sin(Omega2)* cos(i2) *cos(omega2)
	y2j=-sin(Omega2) *sin(omega2) + cos(Omega2) *cos(i2) *cos(omega2)
	y2k=sin(i2)* cos(omega2)
	z2i=sin(i2)* sin(Omega2)
	z2j=-sin(i2) *cos(Omega2)
	z2k=cos(i2)
	
	x1hat= [x1i,x1j,x1k]	
	y1hat= [y1i,y1j,y1k] 
	z1hat= [z1i,z1j,z1k]

	x2hat= [x2i,x2j,x2k]	
	y2hat= [y2i,y2j,y2k] 
	z2hat= [z2i,z2j,z2k]

	#radial unit vector
	r1hat= multiply(x1hat,(cos(E1) - e1)/(1. -e1 *cos(E1))) +  multiply(y1hat,sqrt(1. - e1**2)* sin(E1)/(1. -e1 *cos(E1)))
	#tangential unit vector		
	t1hat= multiply(y1hat,(cos(E1) - e1)/(1. -e1 *cos(E1))) +  multiply(x1hat,-sqrt(1. - e1**2)* sin(E1)/(1. -e1 *cos(E1)))

	r1=a1*(1- e1*cos(E1))
	r1vec=multiply(r1hat, r1) #position vector of m1

	F0=-r1vec-multiply(x2hat,a2*e2)
	F1=multiply(y2hat,a2 *sqrt(1.-e2**2))
	F2 =multiply(x2hat,a2) 
	
	Ab=dot(r1vec,r1vec) + a2**2 + b**2 + 2* a2* e2* dot(r1vec,x2hat) 
	Bsineps= a2* sqrt(1. - e2**2) *dot(r1vec,y2hat);
	Bcoseps= a2*dot(r1vec,x2hat) + a2**2 *e2
	B=sqrt(Bsineps**2+Bcoseps**2)
	c = a2**2* e2**2
	
	Q= 1/9.* (c - Ab)**2 - 1/3. *(B**2 - Ab*c)
	R= 1/27. *(c - Ab)**3 - 1/6. *(c - Ab)* (B**2 - Ab* c) + 1/2.* c *Bsineps**2
	theta = arccos(R/(sqrt(Q**3)))
	
	lam0 = (Ab - c)/3. - 2.*sqrt(Q)*cos((2.*pi)/3. + theta/3.)
	lam1 = (Ab - c)/3. - 2.*sqrt(Q)*cos((-2.*pi)/3. + theta/3.)
	lam2 = (Ab - c)/3. - 2.*sqrt(Q)*cos(theta/3.)
	#Define the terms Qij  by taking the square of the definitions in the paper
	#If we do not define in this way, sqrt(Q01) would give imaginary result.
	Q00 = (lam0*(c + lam0))/((lam0 - lam1)*(lam0 - lam2)) #these are square of  Qij
	Q01 = (lam1*(c + lam1))/((lam0 - lam1)*(lam1 - lam2))
	Q02 = (c + lam2)/((lam0 - lam2)*(lam1 - lam2))*abs(lam2)
	Q10 = (c + lam0)/(lam0*(lam0 - lam1)*(lam0 - lam2)) #Bsineps
	Q11 =(c + lam1)/((lam0 - lam1)*lam1*(lam1 - lam2)) #Bsineps
	Q12 = (c + lam2)/((lam0 - lam2)*(lam1 - lam2)*abs(lam2)) #-Bsineps
	Q20 = lam0/((c + lam0)*(lam0 - lam1)*(lam0 - lam2)) #Bcoseps
	Q21 = lam1/((lam0 - lam1)*(c + lam1)*(lam1 - lam2)) #Bcoseps
	Q22 = abs(lam2)/((lam0 - lam2)*(lam1 - lam2)*(c + lam2))#Bcoseps
	
	#Now,their multiplication e.g Q02*Q22 will be positive, so sqrt(Q02*Q22) will be positive, as well.
	U0=Q00 + Q02 - e2*sqrt(Q00*Q20)*Bcoseps - e2*sqrt(Q02*Q22)*Bcoseps
	U1=sqrt(Q00*Q10)*Bsineps - sqrt(Q02*Q12)*Bsineps - e2*sqrt(Q10*Q20)*Bsineps*Bcoseps + e2*sqrt(Q12*Q22)*Bcoseps*Bsineps
	U2=sqrt(Q00*Q20)*Bcoseps - e2*Q20*Bcoseps**2 + sqrt(Q02*Q22)*Bcoseps - e2*Q22*Bcoseps**2
	
	V0=Q01 - Q02 - e2*sqrt(Q01*Q21)*Bcoseps + e2*sqrt(Q02*Q22)*Bcoseps
	V1=sqrt(Q01*Q11)*Bsineps + sqrt(Q02*Q12)*Bsineps - e2*sqrt(Q11*Q21)*Bsineps*Bcoseps - e2*sqrt(Q12*Q22)*Bsineps*Bcoseps
	V2=sqrt(Q01*Q21)*Bcoseps - e2*Q21*Bcoseps*Bcoseps - sqrt(Q02*Q22)*Bcoseps + e2*Q22*Bcoseps**2
	
	FU=multiply(F0,U0)+multiply(F1,U1)+multiply(F2,U2)
	FV=multiply(F0,V0)+multiply(F1,V1)+multiply(F2,V2)
	
	k=sqrt((lam1 - lam2)/(lam0 - lam2))
	#the acceleration vector is
	acc2=(2.*G*sqrt(lam0 - lam2)*m2/((lam0 - lam1)*(lam1 - lam2)*pi))*((FV + multiply(FU,k**2))*EllE(k) - FV*(1 - k**2)*EllK(k))
	#radial,tangential and normal components of the acceleration vector:
	R2=dot(acc2,r1hat)
	W2=dot(acc2,z1hat)	
	S2=dot(acc2,t1hat)
	return R2,W2,S2

#accelation1 function calculates the averaged acceleration of m2 analytically due the effect of m1
def acceleration1(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,E2,a2):
	
	x1i=cos(Omega1)* cos(omega1) - sin(Omega1) *cos(i1)* sin(omega1)
	x1j=sin(Omega1)* cos(omega1) + cos(Omega1) *cos(i1) *sin(omega1)      	
	x1k=sin(i1)*sin(omega1)
	y1i=-cos(Omega1) *sin(omega1) - sin(Omega1)* cos(i1) *cos(omega1)
	y1j=-sin(Omega1) *sin(omega1) + cos(Omega1) *cos(i1) *cos(omega1)
	y1k=sin(i1)* cos(omega1)
	z1i=sin(i1)* sin(Omega1)
	z1j=-sin(i1) *cos(Omega1)
	z1k=cos(i1)
	
	x2i=cos(Omega2)* cos(omega2) - sin(Omega2) *cos(i2)* sin(omega2)
	x2j=sin(Omega2)* cos(omega2) + cos(Omega2) *cos(i2) *sin(omega2)      	
	x2k=sin(i2)*sin(omega2)
	y2i=-cos(Omega2) *sin(omega2) - sin(Omega2)* cos(i2) *cos(omega2)
	y2j=-sin(Omega2) *sin(omega2) + cos(Omega2) *cos(i2) *cos(omega2)
	y2k=sin(i2)* cos(omega2)
	z2i=sin(i2)* sin(Omega2)
	z2j=-sin(i2) *cos(Omega2)
	z2k=cos(i2)

	x1hat= [x1i,x1j,x1k]	
	y1hat= [y1i,y1j,y1k] 
	z1hat= [z1i,z1j,z1k]

	x2hat= [x2i,x2j,x2k]	
	y2hat= [y2i,y2j,y2k] 
	z2hat= [z2i,z2j,z2k]

	r2hat= multiply(x2hat,(cos(E2) - e2)/(1. -e2 *cos(E2))) +  multiply(y2hat,sqrt(1. - e2**2)* sin(E2)/(1. -e2 *cos(E2)))
	t2hat= multiply(y2hat,(cos(E2) - e2)/(1. -e2 *cos(E2))) +  multiply(x2hat,-sqrt(1. - e2**2)* sin(E2)/(1. -e2 *cos(E2)))

	r2=a2*(1- e2*cos(E2)) 
	r2vec=multiply(r2hat, r2)

	F0=-r2vec-multiply(x1hat,a1*e1)
	F1=multiply(y1hat,a1 *sqrt(1.-e1**2))
	F2 =multiply(x1hat,a1) 
	
	Ab=dot(r2vec,r2vec) + a1**2 + b**2 + 2* a1* e1* dot(r2vec,x1hat) 
	Bsineps= a1* sqrt(1. - e1**2) *dot(r2vec,y1hat);
	Bcoseps= a1*dot(r2vec,x1hat) + a1**2 *e1
	B=sqrt(Bsineps**2+Bcoseps**2)
	c = a1**2* e1**2
	
	Q= 1/9.* (c - Ab)**2 - 1/3. *(B**2 - Ab*c)
	R= 1/27. *(c - Ab)**3 - 1/6. *(c - Ab)* (B**2 - Ab* c) + 1/2.* c *Bsineps**2
	theta = arccos(R/(sqrt(Q**3)))
	
	lam0 = (Ab - c)/3. - 2.*sqrt(Q)*cos((2.*pi)/3. + theta/3.)
	lam1 = (Ab - c)/3. - 2.*sqrt(Q)*cos((-2.*pi)/3. + theta/3.)
	lam2 = (Ab - c)/3. - 2.*sqrt(Q)*cos(theta/3.)
	
	Q00 = (lam0*(c + lam0))/((lam0 - lam1)*(lam0 - lam2)) #these are square of  Qij
	Q01 = (lam1*(c + lam1))/((lam0 - lam1)*(lam1 - lam2))
	Q02 = (c + lam2)/((lam0 - lam2)*(lam1 - lam2))*abs(lam2)
	Q10 = (c + lam0)/(lam0*(lam0 - lam1)*(lam0 - lam2)) #Bsineps
	Q11 =(c + lam1)/((lam0 - lam1)*lam1*(lam1 - lam2)) #Bsineps
	Q12 = (c + lam2)/((lam0 - lam2)*(lam1 - lam2)*abs(lam2)) #-Bsineps
	Q20 = lam0/((c + lam0)*(lam0 - lam1)*(lam0 - lam2)) #Bcoseps
	Q21 = lam1/((lam0 - lam1)*(c + lam1)*(lam1 - lam2)) #Bcoseps
	Q22 = abs(lam2)/((lam0 - lam2)*(lam1 - lam2)*(c + lam2))#Bcoseps
	
	
	U0=Q00 + Q02 - e1*sqrt(Q00*Q20)*Bcoseps - e1*sqrt(Q02*Q22)*Bcoseps
	U1=sqrt(Q00*Q10)*Bsineps - sqrt(Q02*Q12)*Bsineps - e1*sqrt(Q10*Q20)*Bsineps*Bcoseps + e1*sqrt(Q12*Q22)*Bcoseps*Bsineps
	U2=sqrt(Q00*Q20)*Bcoseps - e1*Q20*Bcoseps**2 + sqrt(Q02*Q22)*Bcoseps - e1*Q22*Bcoseps**2
	
	V0=Q01 - Q02 - e1*sqrt(Q01*Q21)*Bcoseps + e1*sqrt(Q02*Q22)*Bcoseps
	V1=sqrt(Q01*Q11)*Bsineps + sqrt(Q02*Q12)*Bsineps - e1*sqrt(Q11*Q21)*Bsineps*Bcoseps - e1*sqrt(Q12*Q22)*Bsineps*Bcoseps
	V2=sqrt(Q01*Q21)*Bcoseps - e1*Q21*Bcoseps*Bcoseps - sqrt(Q02*Q22)*Bcoseps + e1*Q22*Bcoseps**2
	
	FU=multiply(F0,U0)+multiply(F1,U1)+multiply(F2,U2)
	FV=multiply(F0,V0)+multiply(F1,V1)+multiply(F2,V2)
	
	k=sqrt((lam1 - lam2)/(lam0 - lam2))
	
	acc1=(2.*G*sqrt(lam0 - lam2)*m1/((lam0 - lam1)*(lam1 - lam2)*pi))*((FV + multiply(FU,k**2))*EllE(k) - FV*(1 - k**2)*EllK(k))
		
	R1=dot(acc1,r2hat)
	W1=dot(acc1,z2hat)	
	S1=dot(acc1,t2hat)
	return R1,W1,S1


##Calculate the components of the averaged acceleration vector of m1 due to m2
def RSW1(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2):
	# initial Fourier coefficients
	rcn0 = 0
	rcn1 = 0
	rcn2 = 0
	rsn1 = 0
	rsn2 = 0
	scn0 = 0
	scn1 = 0
	scn2 = 0
	ssn1 = 0
	ssn2 = 0
	wcn0 = 0
	wcn1 = 0
	wcn2 = 0
	wsn1 = 0
	wsn2 = 0

	
	
      	#discrete FT
	for j in range(N1):
               
		E1=2.*pi*j/N1
		R2,W2,S2=acceleration2(i1,omega1,Omega1,e1,E1,a1,i2,omega2,Omega2,e2,a2)	
		
		rcn0+= (1./N1)*R2* cos(0*2*pi*j/N1) 
		rcn1+= (1./N1)*R2* cos(1*2*pi*j/N1) 
		rcn2+= (1./N1)*R2* cos(2*2*pi*j/N1) 
		rsn1+= (1./N1)*R2* sin(1*2*pi*j/N1) 
		rsn2+= (1./N1)*R2* sin(2*2*pi*j/N1) 
		scn0+= (1./N1)*S2* cos(0*2*pi*j/N1) 
		scn1+= (1./N1)*S2* cos(1*2*pi*j/N1) 
		scn2+= (1./N1)*S2* cos(2*2*pi*j/N1) 
		ssn1+= (1./N1)*S2* sin(1*2*pi*j/N1) 
		ssn2+= (1./N1)*S2* sin(2*2*pi*j/N1) 
		wcn0+= (1./N1)*W2* cos(0*2*pi*j/N1) 
		wcn1+= (1./N1)*W2* cos(1*2*pi*j/N1) 
		wcn2+= (1./N1)*W2* cos(2*2*pi*j/N1) 
		wsn1+= (1./N1)*W2* sin(1*2*pi*j/N1) 
		wsn2+= (1./N1)*W2* sin(2*2*pi*j/N1) 
	return rcn0,rcn1,rcn2,rsn1,rsn2,scn0,scn1,scn2,ssn1,ssn2,wcn0,wcn1,wcn2,wsn1,wsn2

#Calculate the components of the averaged acceleration vector of m2 due to m1
def RSW2(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2):
	
	rcn0 = 0
	rcn1 = 0
	rcn2 = 0
	rsn1 = 0
	rsn2 = 0
	scn0 = 0
	scn1 = 0
	scn2 = 0
	ssn1 = 0
	ssn2 = 0
	wcn0 = 0
	wcn1 = 0
	wcn2 = 0
	wsn1 = 0
	wsn2 = 0

	
	
      	
	for j in range(N2):
               
		E2=2.*pi*j/N2
		R1,W1,S1=acceleration1(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,E2,a2)	
		
	
		rcn0+= (1./N2)*R1* cos(0*2*pi*j/N2) 
		rcn1+= (1./N2)*R1* cos(1*2*pi*j/N2) 
		rcn2+= (1./N2)*R1* cos(2*2*pi*j/N2) 
		rsn1+= (1./N2)*R1* sin(1*2*pi*j/N2) 
		rsn2+= (1./N2)*R1* sin(2*2*pi*j/N2) 
		scn0+= (1./N2)*S1* cos(0*2*pi*j/N2) 
		scn1+= (1./N2)*S1* cos(1*2*pi*j/N2) 
		scn2+= (1./N2)*S1* cos(2*2*pi*j/N2) 
		ssn1+= (1./N2)*S1* sin(1*2*pi*j/N2) 
		ssn2+= (1./N2)*S1* sin(2*2*pi*j/N2) 
		wcn0+= (1./N2)*W1* cos(0*2*pi*j/N2) 
		wcn1+= (1./N2)*W1* cos(1*2*pi*j/N2) 
		wcn2+= (1./N2)*W1* cos(2*2*pi*j/N2) 
		wsn1+= (1./N2)*W1* sin(1*2*pi*j/N2) 
		wsn2+= (1./N2)*W1* sin(2*2*pi*j/N2) 
	return rcn0,rcn1,rcn2,rsn1,rsn2,scn0,scn1,scn2,ssn1,ssn2,wcn0,wcn1,wcn2,wsn1,wsn2

#Calculate the rates of change of the angular momentum and eccentricity vector of m1

#dL/dt=N, dA/dt=Adat
def NA1dat(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2):
      
	
	rcn0,rcn1,rcn2,rsn1,rsn2,scn0,scn1,scn2,ssn1,ssn2,wcn0,wcn1,wcn2,wsn1,wsn2=RSW1(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2)
	#components of the torque in m1's orbital frame (e.g N1=N1xp x1^+ N1yp y1^ + N1zp z1^)		
	N1xp= (1./(n1*a1 ** 2))*a1* sqrt(1 - e1 ** 2)*(wsn1 - 1./2*e1 *wsn2)
	N1yp= - (1./(n1*a1 ** 2))*a1* ((1 + e1 ** 2)*wcn1 - 3./2* e1* wcn0 - 1./2*e1* wcn2)
	N1zp= a1*(1./(n1*a1 ** 2))  *((1 + 1./2*e1 ** 2)*scn0 - 2*e1 *scn1 + 1./2 *e1 ** 2 *scn2)

	Adat1xp=(sqrt(1 - e1 ** 2))/(2*n1*a1)*((4 *scn1 - e1 *scn2 - 3*e1 *scn0) + 2*sqrt (1 - e1 ** 2)*rsn1)
	Adat1yp=(1./(2*n1*a1))*(2*(2 - e1 ** 2) *ssn1 - e1* ssn2 - 2*sqrt(1 - e1 ** 2)*(rcn1 - e1* rcn0))
	Adat1zp=(-e1/(n1*a1))*(wsn1 - 1./2*wsn2)

	
	return N1xp,N1yp,N1zp,Adat1xp,Adat1yp,Adat1zp

#Calculate the rate of change of the angular momentum and eccentricity vector of m2
def NA2dat(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2):
      
	
	rcn0,rcn1,rcn2,rsn1,rsn2,scn0,scn1,scn2,ssn1,ssn2,wcn0,wcn1,wcn2,wsn1,wsn2=RSW2(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2)
	N2xp= (1./(n2*a2 ** 2))*a2* sqrt(1 - e2 ** 2)*(wsn1 - 1./2*e2 *wsn2)
	N2yp= - (1./(n2*a2 ** 2))*a2* ((1 + e2 ** 2)*wcn1 - 3./2* e2* wcn0 - 1./2*e2* wcn2)
	N2zp= a2*(1./(n2*a2 ** 2))  *((1 + 1./2*e2 ** 2)*scn0 - 2*e2 *scn1 + 1./2 *e2 ** 2 *scn2)

	Adat2xp=(sqrt(1 - e2 ** 2))/(2*n2*a2)*((4 *scn1 - e2 *scn2 - 3*e2 *scn0) + 2*sqrt (1 - e2 ** 2)*rsn1)
	Adat2yp=(1./(2*n2*a2))*(2*(2 - e2 ** 2) *ssn1 - e2* ssn2 - 2*sqrt(1 - e2 ** 2)*(rcn1 - e2* rcn0))
	Adat2zp=(-e2/(n2*a2))*(wsn1 - 1./2*wsn2)

	
	return N2xp,N2yp,N2zp,Adat2xp,Adat2yp,Adat2zp

#Angular momentum function
def angularmom(i,Omega,ep):

	zi=sin(i)* sin(Omega)
	zj=-sin(i) *cos(Omega)
	zk=cos(i)
	j=sqrt(1-ep**2)
	jx=j*zi 
	jy=j*zj
	jz=j*zk
	return jx,jy,jz

#Calculate the rates of change of the orbital parameters of m1 and m2
def f(r):
	i1=r[0]
	Omega1=r[1]
	omega1=r[2]
	e1=r[3]
	i2=r[4]
	Omega2=r[5]
	omega2=r[6]
	e2=r[7]
	#for m1:
	N1xp,N1yp,N1zp,Adat1xp,Adat1yp,Adat1zp=NA1dat(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2)
	
	N1vec=array([N1xp,N1yp,N1zp],float)
	Adat1vec=array([Adat1xp,Adat1yp,Adat1zp],float)	
	L1vec=array([0,0,sqrt(1-e1**2)],float) #angular momentum vector of m1

	fi3,fOmega3,fomega3,fe3=f3(i1,omega1,Omega1,e1,a1)

	#the rate of change of the orbital parameters of m1 due to m2 and m
	fi1=fi3+(-sin(omega1)*dot(xhat,N1vec)-cos(omega1)*dot(yhat,N1vec))/norm(L1vec)
	fOmega1=fOmega3+(cos(omega1)*dot(xhat,N1vec)-sin(omega1)*dot(yhat,N1vec))/(norm(L1vec)*sin(i1))
	fomega1=fomega3+dot(Adat1vec,yhat)/e1-cos(i1)*(cos(omega1)*dot(xhat,N1vec)-sin(omega1)*dot(yhat,N1vec))/(norm(L1vec)*sin(i1))
	fe1=fe3+dot(Adat1vec,xhat)

	#for m2:
	N2xp,N2yp,N2zp,Adat2xp,Adat2yp,Adat2zp=NA2dat(i1,omega1,Omega1,e1,a1,i2,omega2,Omega2,e2,a2)
	
	N2vec=array([N2xp,N2yp,N2zp],float)
	Adat2vec=array([Adat2xp,Adat2yp,Adat2zp],float)	
	L2vec=array([0,0,sqrt(1-e2**2)],float) #angular momentum vector of m2

	fi3,fOmega3,fomega3,fe3=f3(i2,omega2,Omega2,e2,a2)

	#the rate of change of the orbital parameters of m2 due to m1 and m
	fi2=fi3+(-sin(omega2)*dot(xhat,N2vec)-cos(omega2)*dot(yhat,N2vec))/norm(L2vec)
	fOmega2=fOmega3+(cos(omega2)*dot(xhat,N2vec)-sin(omega2)*dot(yhat,N2vec))/(norm(L2vec)*sin(i2))
	fomega2=fomega3+dot(Adat2vec,yhat)/e2-cos(i2)*(cos(omega2)*dot(xhat,N2vec)-sin(omega2)*dot(yhat,N2vec))/(norm(L2vec)*sin(i2))
	fe2=fe3+dot(Adat2vec,xhat)
	
	
	
	return array([fi1,fOmega1,fomega1,fe1,fi2,fOmega2,fomega2,fe2],float)



tpoints=arange(0,totaltime+H,H)


i1points=[]
Omega1points=[]
omega1points=[]
e1points=[]


i2points=[]
Omega2points=[]
omega2points=[]
e2points=[]

r=array([i10,Omega10,omega10,e10,i20,Omega20,omega20,e20],float)

print('#t\t e1\t i1\t omega1\t Omega1\t error \t e2\t i2\t omega2\t Omega2\t')

#Evolve the rings in time using a Bulirsch-Stoer (BS) integrator:
for t in tpoints:
        
	i1points.append(r[0])
	Omega1points.append(r[1])
	omega1points.append(r[2])
	e1points.append(r[3])

	i2points.append(r[4])
	Omega2points.append(r[5])
	omega2points.append(r[6])
	e2points.append(r[7])
	Irel=arccos((cos(r[0])*cos(r[4])+sin(r[0])*sin(r[4])*cos(r[1]-r[5])))*180./pi #relative inclination of m1 and m2
	
	j1x,j1y,j1z=angularmom(r[0],r[1],r[3]) #components of the angular momentum vector of m1 
	j2x,j2y,j2z=angularmom(r[4],r[5],r[7]) #components of the angular momentum vector of m2

	print(t,r[3],r[0]*180./pi,r[2]*180./pi,r[1]*180./pi,r[7],r[4]*180./pi,r[6]*180./pi,r[5]*180./pi,error,Irel,j1x,j1y,j1z,j2x,j2y,j2z)
	
	error=2*delta #update the error

#do one modified midpoint step of size H	
	r1=r+0.5*H*f(r)
		
	
       
	r2=r+H*f(r1)
	
	R1=empty([1,8],float)
	R1[0]=0.5*(r1+r2+0.5*H*f(r2))
#print(R1[0])
#Now increase n untill the required accuracy is reached
	n=1
	while error>delta:
		n+=1
		h=H/n
		#Modified midpoint method
		
		r1=r+0.5*h*f(r)
		
		r2=r+h*f(r1)

		for ii in range(n-1):
			

			r1+=h*f(r2)


			r2+=h*f(r1)
	

		#calculate extrapolation estimates
		
		R2=R1
		R1=empty([n,8],float)
		R1[0]=0.5*(r1+r2+0.5*h*f(r2))
		#print(R1[0])
		for mm in range(1,n):
			divide=((n/(n-1.))**(2*mm)-1.)
			epsilon=(R1[mm-1]-R2[mm-1])/divide
			R1[mm]=R1[mm-1]+epsilon
			#print(R1[m])
		error=abs(epsilon[0])
		#print(error)
#set r equal to the most accurate estimate we have before moving on to the next step
	r=R1[n-1]
	



