import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

tol = 1e-5
filename = "stiffness_b70_w_200_coex_int_K_per_tperp_U12.dat"
filename_to_write = "stiffness/stiffness_b70_w_200_coex_int_K_per_tperp_U12__test_.dat"

def fitting_function(x:float, a:float, b:float, c:float, d:float, e:float, f:float, g:float, h:float) -> float:
    return a*x + b*x**2 + c*x**3 + d*x**4 + e*x**5 + f*x**6 + g*x**7 + h

NOCOEX_data = np.loadtxt("stiffness/stiffness_b70_w_200_nocoex_int_K_AFM_per_tperp_U12.dat",usecols=(1,2))
print("NOCOEX_data: ",NOCOEX_data)
AFM_SC_NOCOEX = np.loadtxt("stiffness/stiffness_b70_w_200_coex_AFM_SC_1_int_K_per_tperp_U12.dat",usecols=(1,2))
print("AFM_SC_NOCOEX: ",AFM_SC_NOCOEX)
dop_M_data = np.loadtxt("COEX/U12/Loop_COEX_tmp.dat", skiprows=1, usecols=(2,4))
print("dop_M_data: ",dop_M_data)
COEX_data = np.loadtxt(filename, usecols=(1,2))
print("COEX_data: ",COEX_data)
COEX_data_cnc = np.loadtxt(filename_to_write, usecols=(1,2))

print("COEX_data_cnc",COEX_data_cnc)

#dop_M_data = dop_M_data.tolist()
#print(len(dop_M_data), "\n\n")
dop_M_data_new = []
for element in dop_M_data:
    element[1] = float(element[1])
    print(element[0],element[1])
    if element[1] >= tol:
        dop_M_data_new.append(element)
dop_M_data_new = np.asarray(dop_M_data_new)
min_dop_AFM = min(dop_M_data_new[:,0])
max_dop_AFM = max(dop_M_data_new[:,0])

print(min_dop_AFM,max_dop_AFM)

comp = np.hstack((AFM_SC_NOCOEX,NOCOEX_data))

comp_new = []
for element in comp:
    if min_dop_AFM <= element[0] <= max_dop_AFM:
        comp_new.append(element[0:4])
comp_new = np.asarray(comp_new)

popt, pcov = curve_fit(fitting_function, comp_new[:,0], comp_new[:,1])
popt2, pcov2 = curve_fit(fitting_function, comp_new[:,2], comp_new[:,3])

COEX_data_new  = []
for element in COEX_data:
    if min_dop_AFM <= element[0] <= max_dop_AFM:
        COEX_data_new.append([element[0],element[1]*(fitting_function(element[0],*popt2)/fitting_function(element[0],*popt))])
COEX_data_new = np.asarray(COEX_data_new)

for element in COEX_data_cnc:
    print("element_cnc = ", element)
    for element2 in COEX_data_new:
        print("element_new = ",element2)
        if element[0] == element2[0]:
            print("element = ",element)
            element[1] = element2[1]
print("COEX_data_cnc = ",COEX_data_cnc)

#print("comp_new = ",comp_new,"\n\n\n")
#print("data = ", AFM_SC_NOCOEX,"\t", NOCOEX_data)
print("COEX_data_new = ",COEX_data_new)
#print(dop_M_data_new, "\n\n")
#print(dop_M_data[:,0])
with open("COEX_data_new_"+filename,'w') as fi:
    for i in range(len(COEX_data_new)):
        fi.write("{0:5.5f}\t\t{1}\n".format(COEX_data_new[i,0],COEX_data_new[i,1]))
fi.close()

with open(filename_to_write,'w') as fi:
    for i in range(len(COEX_data_cnc)):
        fi.write("{0:5.5f}\t\t{1}\n".format(COEX_data_cnc[i,0],COEX_data_cnc[i,1]))
fi.close()

plt.plot(comp_new[:,2], fitting_function(comp_new[:,2], *popt2), '.', label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f, f=%5.3f, g=%5.3f, h=%5.3f' % tuple(popt2))
plt.plot(comp_new[:,0], fitting_function(comp_new[:,0], *popt), '.',label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f, f=%5.3f, g=%5.3f, h=%5.3f' % tuple(popt))
plt.plot(COEX_data_new[:,0],COEX_data_new[:,1],'.',label='COEX_data_new')
plt.plot(comp_new[:,0], comp_new[:,1],'.',label='data_AFM_SC_NOCOEX')
plt.plot(comp_new[:,2], comp_new[:,3], '.', label='data_nocoex')
plt.plot(COEX_data_cnc[:,0],COEX_data_cnc[:,1],'.',label='full_data')
plt.xlabel(r"Density of particles $n$")
plt.ylabel(r"Superfluid stiffness $\rho_s$")
plt.legend()

plt.show()

