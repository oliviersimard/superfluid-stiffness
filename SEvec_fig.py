import numpy as np
import matplotlib.pyplot as plt
import os
import re
from fnmatch import fnmatch

paths = './'
pattern = "*.npy"
file_array = []
list_temp = []
file_dict = {}
#tt=1500 #Pay attention to these values
tt=400
beta=50
key_fold = "./NOCOEX/U8/SEvec/SEvec_b50_SC/SEvec_b50_SC" #folder in which the program will go fetch the data
MU=30
MU2=31
MU3=32
MU4=33
w_list = [(2*n+1)*np.pi/beta for n in range(tt)]

for path, subdirs, files in os.walk(paths):
    for name in files:
        if key_fold in path:
            if fnmatch(name, pattern):
                file_array.append(os.path.join(path,name))
print(file_array)
regex = re.compile(r'\d+')

for el in file_array:
    list_temp.append(int(regex.findall(el)[-1]))

list_im_list = []
list_re_list = []

sorted_list = sorted(zip(list_temp,file_array))
mu_list = np.loadtxt('NOCOEX/U8/stiffness_b50.dat',skiprows=0,usecols=(1,))
for i, path_to_files in enumerate(sorted_list):
    file_dict[i] = np.load(path_to_files[1])
    imag_SEpart = []
    real_SEpart = []
    for j in range(file_dict[i].shape[1]):
        imag_SEpart.append(np.imag(file_dict[i][0,j,0,0]).tolist()) #Change the element taken from self-energy with numbers after j
        real_SEpart.append(np.real(file_dict[i][0,j,0,0]).tolist())
    #print(list(zip(w_list,real_SEpart)))
    list_im_list.append(zip(w_list,imag_SEpart))
    list_re_list.append(zip(w_list,real_SEpart))


#print(list(list_im_list[1]))


filename = "IM_vs_RE_00_GF_NOCOEX_U8"#_%.5f"%(-1.0*(mu_list[MU]-1.0)) #Change the chemical potential here

fig = plt.figure()

ax = plt.subplot(211)

#xlabel = (r"$\omega$")
ylabel = (r"Im$\Sigma_{00}$"r"$\left(\omega\right)$")

#plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

self_array_im = []
for elem in list_im_list[MU]:
    self_array_im.append(elem[1])

self_array_im2 = []
for elem in list_im_list[MU2]:
    self_array_im2.append(elem[1])

self_array_im3 = []
for elem in list_im_list[MU3]:
    self_array_im3.append(elem[1])

self_array_im4 = []
for elem in list_im_list[MU4]:
    self_array_im4.append(elem[1])

bbox_props = dict(boxstyle="square,pad=0.3", fc="yellow", ec="b", lw=2) #used for the box of text in plot
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9]) #WHEN ONLY PLOTTING ONE GRAPH' UNCOMMENT THESE LINES
#ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)

#ax.plot(w_list, self_array_re, linestyle="-", marker='o', markersize=3, linewidth=2, color='green', label='Re')
ax.plot(w_list, self_array_im, linestyle="-", marker='o', markersize=3, linewidth=2, color='red')
ax.plot(w_list, self_array_im2, linestyle="-", marker='o', markersize=3, linewidth=2, color='green')
ax.plot(w_list, self_array_im3, linestyle="-", marker='o', markersize=3, linewidth=2, color='blue')
ax.plot(w_list, self_array_im4, linestyle="-", marker='o', markersize=3, linewidth=2, color='cyan')
#plt.text(0.8, 0.8, r"$\delta=%.4f$"%(-1.0*(mu_list[MU]-1.0)), transform=ax.transAxes, bbox=bbox_props)
#ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), fancybox=True, shadow=True, ncol=3)
plt.xlim([-0.01,30.01])
plt.xticks(np.arange(0,31,1))#, rotation='vertical')
plt.title(r"Self-energy $\omega$-dependence at fixed doping for $\beta = 50$") # at  " r"$\mu=%.2f$"%mu[0][0])

##############################################################
ax2 = plt.subplot(212) #ELIMINATE THIS SUBPLOT TO HAVE ONLY ONE GRAPH

xlabel = (r"$i\omega$")
ylabel = (r"Re$\Sigma_{00}$"r"$\left(\omega\right)$")

plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)


self_array_re = []
for elem in list_re_list[MU]:
    self_array_re.append(elem[1])

self_array_re2 = []
for elem in list_re_list[MU2]:
    self_array_re2.append(elem[1])

self_array_re3 = []
for elem in list_re_list[MU3]:
    self_array_re3.append(elem[1])

self_array_re4 = []
for elem in list_re_list[MU4]:
    self_array_re4.append(elem[1])

box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

ax2.plot(w_list, self_array_re, linestyle="-", marker='o', markersize=3, linewidth=2, color='red',label=r"$\delta=%.5f$"%(mu_list[MU]))#, label="Re")
ax2.plot(w_list, self_array_re2, linestyle="-", marker='o', markersize=3, linewidth=2, color='green',label=r"$\delta=%.5f$"%(mu_list[MU2]))#, label="Re")
ax2.plot(w_list, self_array_re3, linestyle="-", marker='o', markersize=3, linewidth=2, color='blue',label=r"$\delta=%.5f$"%(mu_list[MU3]))
ax2.plot(w_list, self_array_re4, linestyle="-", marker='o', markersize=3, linewidth=2, color='cyan',label=r"$\delta=%.5f$"%(mu_list[MU4]))

ax2.legend(loc='upper center', bbox_to_anchor=(0.5,-0.2), fancybox=True, shadow=True, ncol=3)

#plt.text(0.8, 0.8, r"$\delta=%.4f$"%(-1.0*(mu_list[MU]-1.0)), transform=ax2.transAxes, bbox=bbox_props)
plt.xlim([-0.01,30.01])
plt.xticks(np.arange(0,31,1))#, rotation='vertical')

plt.show()
fig.savefig(filename + ".pdf", format='pdf')

################################################################
##Printing Imag and real parts 
################################################################
with open('self_energy.dat', 'w') as l:
    for i in range(len(w_list)):
        l.write("{0:6.5f}\t\t{1:6.5f}\t\t{2:6.5f}\t\t{3:3.2f}\n".format(w_list[i],self_array_re[i],self_array_im[i],mu_list[MU]))

l.close()


