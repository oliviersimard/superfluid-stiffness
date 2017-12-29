import numpy as np
import pylab as plt
import sys
import getopt
import scipy.integrate as sci

def get_option():
	try:
		name = ''
		wmin = -6
		wmax = 6
		nbrw = -1
		option, test = getopt.getopt(sys.argv[1:], 'f:l:m:n:')
	except getopt.GetoptError:
		print('options for analytical continuation:\n-f <inputFile>\n-l <lowest real frequency>\n-m <largest real frequency>\n-n <number of frequencies between wmin and wmax')
		sys.exit()
	for opt, arg in option:
		if opt =='-f':
			name = arg
		elif opt =='-l':
			wmin = float(arg)
		elif opt == '-m':
			wmax = float(arg)
		elif opt== '-n':
			nbrw = int(arg)
	if nbrw <0:
		print('number of real frequency not found. Standard value taken: nbr=100')
		number =100	
	return name, wmin, wmax, nbrw

name, wmin, wmax, nbr_w = get_option()
R = np.loadtxt(name)
omega_n = R[:,0]
Rf_n = R[:, 1]
If_n = R[:,2]
eta = 0.001


f_n = Rf_n +1j*If_n
N = len(omega_n)-7
g = np.zeros((N, N), dtype = complex)
omega_n = omega_n[0:N]

f_n = f_n[0:N]
g[0,:] = f_n[:]
g[1,:] = (f_n[0] - f_n)/((1j*omega_n -1j*omega_n[0])*f_n)
for k in range(2, N):
	g[k, :] = (g[k-1, k-1] - g[k-1, :])/((1j*omega_n - 1j*omega_n[k-1])*g[k-1, :])

a = np.diagonal(g)

A = np.zeros((N, ), dtype = complex)
B = np.zeros((N, ), dtype = complex)
P = np.zeros((N, ), dtype = complex)
fw = []
omega = [ wmin+i*(wmax-wmin)/(nbr_w-1.) for i in range(nbr_w)]

for k  in range (len(omega)):
	z = omega[k] +1j*eta
	P[0] = a[0]
	A[1] = 1.
	B[1]= 1.+a[1]*(z - 1j*omega_n[0]);
	P[1] = P[0]*A[1]/B[1]
	for c in range(2, N):
		A[c] = 1. + a[c]*(z -1j*omega_n[c-1])/A[c-1]
		B[c] = 1. + a[c]*(z - 1j*omega_n[c-1])/B[c-1]
		P[c] = P[c-1]*A[c]/B[c]
	#print P[-1][0]
	fw.append(P[-1])
f = open("analytic_continuation_"+name, "w")
Spec_func = [fw[i].imag for i in range(len(fw))]

print("norm of the spectral function = ", -(1./np.pi)*sci.simps(Spec_func, omega))


for i in range(len(fw)):
	f.write(str(omega[i]) + "\t\t" +str(fw[i].real) + "\t\t" + str(fw[i].imag)+"\n")

f.close()
