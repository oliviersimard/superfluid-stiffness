import numpy as np
import matplotlib.pyplot as plt
import os
import re
from fnmatch import fnmatch
import itertools
import copy
import warnings

paths = './'
pattern = "*.npy"
file_dict = {}
#tt=1500 #Pay attention to these values
tt=200
U=8
beta=60
key_fold1 = "./COEX/U8/SEvec/w_400/SEvec_b60_SC_AFM/SEvec_b60_SC_AFM"#/SEvec_b60_SC_AFM"#/SEvec_b70_SC" #folder in which the program will go fetch the data
key_fold2 = "./COEX/U8/SEvec/w_400/SEvec_b60_SC_AFM/SEvec_b60_SC_AFM_s"
same_repository = False ##Very important to change for file loading
dop_min = 0.8
dop_max = 1.2
decide_curves = 0 ## Input 1 to state that you want certain number equally-spaced self-energy curves to be plotted
num_curves = 1 ## Number of equally-spaced self-energy curves to be plotted
list_dop_curves = [0,5,10,20] ## If decide_curves = 1, type in relevant curves you want to plot
w_list = [(2*n+1)*np.pi/beta for n in range(tt)]
#w_list = [x for x in np.linspace(-10,10,tt)] ###For DOS plots

file_array=[]
for path, subdirs, files in os.walk(paths):
    for name in files:
        if key_fold1 in path:
            if fnmatch(name, pattern):
                file_array.append(os.path.join(path,name))

print("file_array : ", file_array, "\n\n")
regex1 = re.compile(r'(.*?)(?=_)')
regex2 = re.compile(r'(.*?)(?=/)')
print("regex2",regex2.findall(file_array[0]))
file_array_s = []
file_array_d = file_array.copy()
if same_repository:
    if "DOS" in regex2.findall(file_array[0]):
        #print(regex2.findall(file_array[0]))
        key_word = "DOS"
        boool = True

        warnings.warn("All the files share the same repository")
        regex = re.compile(r'\d+')
        list_temp=[]
        for el in file_array:
            list_temp.append(int(regex.findall(el)[-1]))
            
        list_temp_tot = [list(a) for a in zip(list_temp,file_array)]
        list_temp_tot.sort()

        list_temp_tot = np.asarray(list_temp_tot)
    #print("list_temp_tot = ", list_temp_tot, type(list_temp_tot))
    elif "SEvec" in regex2.findall(file_array[0]):
        key_word = "Not_DOS"
        boool = True

        warnings.warn("All the files share the same repository")
        regex = re.compile(r'\d+')
        list_temp=[]
        for el in file_array:
            list_temp.append(int(regex.findall(el)[-1]))
            
        list_temp_tot = [list(a) for a in zip(list_temp,file_array)]
        list_temp_tot.sort()
        list_temp_tot = np.asarray(list_temp_tot)
        print("list_temp_tot",list_temp_tot)
else:
    for el in file_array:
        if "s/SEvec" in regex1.findall(el):   ######Change if looking after Gvec or SEvec
            file_array_s.append(file_array_d.pop(file_array_d.index(el)))
            boool = True
            key_word = "Not_DOS"
        if "s/DOS" in regex1.findall(el):
            file_array_s.append(file_array_d.pop(file_array_d.index(el)))
            boool = True
            key_word = "DOS"
        elif "s/Gvec" in regex1.findall(el):
            file_array_s.append(file_array_d.pop(file_array_d.index(el)))
            boool = False
            key_word = "Not_DOS"

    regex = re.compile(r'\d+')
    list_temp_s=[]
    for el in file_array_s:
        list_temp_s.append(int(regex.findall(el)[-1]))

    list_temp_d=[]
    for el in file_array_d:
        list_temp_d.append(int(regex.findall(el)[-1]))

    list_temp_tot_d = [list(a) for a in zip(list_temp_d,file_array_d)]
    list_temp_tot_d.sort()
    #print("list_temp_tot_d = ", list_temp_tot_d)
    list_temp_tot_s = [list(a) for a in zip(list_temp_s,file_array_s)]
    list_temp_tot_s.sort()
    #print("list_temp_tot_s = ", list_temp_tot_s)
    list_temp_tot = np.vstack((list_temp_tot_d, list_temp_tot_s))
    print("list_temp_tot = ",list_temp_tot)

list_im_list = []
list_re_list = []
list_DoS_1_list = []
list_DoS_2_list = []
list_im_winf = []
list_re_winf = []

#sorted_list = sorted(zip(list_temp,file_array))
mu_list = np.loadtxt('COEX/U8/SEvec/w_400/SEvec_b60_SC_AFM/stiffness.dat',skiprows=0,usecols=(1,))
i = 0
dop_range = []
index_range = []
while i<len(mu_list):
    if mu_list[i]>dop_min and mu_list[i]<dop_max:
        dop_range.append(mu_list[i])
        index_range.append(i)
    i+=1

#print("list_temp_tot = ", list_temp_tot)
if key_word=="DOS":
    for i, path_to_files in enumerate(list_temp_tot):
        if i in index_range:
            file_dict[i] = np.load(path_to_files[1])
            DoS_1 = []
            DoS_2 = []
            for j in range(file_dict[i].shape[1]):
                DoS_1.append(np.real(file_dict[i][0,j,0]).tolist())
                DoS_2.append(np.real(file_dict[i][0,j,1]).tolist())
            list_DoS_1_list.append([list(a) for a in zip(w_list,DoS_1)])
            list_DoS_2_list.append([list(a) for a in zip(w_list,DoS_2)])
else:
    for i, path_to_files in enumerate(list_temp_tot):
        if i in index_range:
            file_dict[i] = np.load(path_to_files[1])
            imag_SEpart = []
            real_SEpart = []
            if 'COEX' in regex2.findall(file_array[0]):
                #print("Length list", len(file_dict[i][0,:,0,0].tolist()))
                list_im_winf.append(np.imag(np.trace(file_dict[i][0,-1,0:2,0:2])).tolist())
                list_re_winf.append(np.real(np.trace(file_dict[i][0,-1,0:2,0:2])).tolist())
            elif 'NOCOEX' in regex2.findall(file_array[0]):
                list_im_winf.append(np.imag(np.trace(file_dict[i][0,-1,0:2,0:2])).tolist())
                list_re_winf.append(np.real(np.trace(file_dict[i][0,-1,0:2,0:2])).tolist())
            #print("list_im_list_winf and list_re_list_winf :", len(list_im_list_winf),"\t",len(list_re_list_winf))
            #print("file_dict shape = ", file_dict[i].shape)
            for j in range(file_dict[i].shape[1]):
                if 'COEX' in regex2.findall(file_array[0]):
                    imag_SEpart.append(np.imag(np.trace(file_dict[i][0,j,0:2,0:2])).tolist()) #Change the element taken from self-energy with numbers after j
                    real_SEpart.append(np.real(np.trace(file_dict[i][0,j,0:2,0:2])).tolist())
                elif 'NOCOEX' in regex2.findall(file_array[0]):
                    imag_SEpart.append(np.imag(file_dict[i][0,j,0,0]).tolist()) #Change the element taken from self-energy with numbers after j
                    real_SEpart.append(np.real(file_dict[i][0,j,0,0]).tolist())
            #print(list(zip(w_list,real_SEpart)))
            list_im_list.append([list(a) for a in zip(w_list,imag_SEpart)])
            list_re_list.append([list(a) for a in zip(w_list,real_SEpart)])

print("winf lists:",list_im_winf,"\t",list_re_winf)
list_im_winf = np.asarray(list_im_winf)
list_re_winf = np.asarray(list_re_winf)
list_im_list = np.asarray(list_im_list)
list_re_list = np.asarray(list_re_list)
list_DoS_1_list = np.asarray(list_DoS_1_list)
list_DoS_2_list = np.asarray(list_DoS_2_list)

print("list_im_list : ", list_im_list.shape, "\t", "list_re_list : ", list_re_list.shape)
#print(np.trace(np.asarray(np.imag(file_dict[15][0,0,:,:]).tolist())))
#print("lengths of list_im_list and list_re_list = ", list_im_list[0][:,1]," and \n\n",list_re_list[0][:,1])

#filename = "IM_vs_RE_00_GF_NOCOEX_U8_beta_70_w_200"#_%.5f"%(-1.0*(mu_list[MU]-1.0)) #Change the chemical potential here
filename = "SEvec_COEX_U8_beta_60_w_400_cluster_trace_0_2"
filename2 = "SEvec_COEX_U8_beta_60_w_400_cluster_winf_trace_0_2"

fig = plt.figure(1)

ax = plt.subplot(211)

#xlabel = (r"$\omega$")
if boool and key_word=="Not_DOS":
    ylabel = (r"Im$\Sigma_{00}$"r"$\left(\omega\right)$")
elif boool and key_word=="DOS":
    ylabel = (r"DOS")
else:
    ylabel = (r"ImG$_{00}$"r"$\left(\omega\right)$")

bbox_props = dict(boxstyle="square,pad=0.3", fc="yellow", ec="b", lw=2) #used for the box of text in plot
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9]) #WHEN ONLY PLOTTING ONE GRAPH' UNCOMMENT THESE LINES
#ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
#plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

dop_list_im_list = list(sorted(zip(dop_range,list_im_list), key=lambda x: x[0]))
dop_list_im_winf = list(sorted(zip(dop_range,list_im_winf),key=lambda x: x[0]))
dop_list_DoS_1_list = list(sorted(zip(dop_range,list_DoS_1_list), key=lambda x: x[0]))
print("dop_list_DoS_1_list : ", len(dop_list_DoS_1_list))
#dop_list_im_list = dop_list_im_list.sort()
#print(dop_list_im_list, " vs \n\n ", list_im_list)
if key_word=="DOS":
    list_DoS_1 = []
    if decide_curves==0:
        eq_spaced_dop = np.arange(0,len(list_DoS_1_list),num_curves)
        for dop in eq_spaced_dop:
            list_DoS_1.append(dop_list_DoS_1_list[dop][1])
    elif decide_curves==1:
        for dop in list_dop_curves:
            list_DoS_1.append(dop_list_DoS_1_list[dop][1])
else:
    list_self_array_im = []
    list_self_array_im_winf = []
    if decide_curves==0:
        eq_spaced_dop = np.arange(0,len(list_im_list),num_curves)
        for dop in eq_spaced_dop:
            list_self_array_im.append(dop_list_im_list[dop][1])
            list_self_array_im_winf.append(dop_list_im_winf[dop][1])
    elif decide_curves==1:
        for dop in list_dop_curves:
            list_self_array_im.append(dop_list_im_list[dop][1])
            list_self_array_im_winf.append(dop_list_im_winf[dop][1])

#print("list_DoS_1 : ", list_DoS_1,"\n")
#print("list_self_array_im = ", list_self_array_im, "\n")
list_self_array_im_winf = np.asarray(list_self_array_im_winf)
#print("list_self_array_im_winf",list_self_array_im_winf)

color_im = iter(plt.cm.rainbow(np.linspace(0,1,200)))
color_re = iter(plt.cm.rainbow(np.linspace(0,1,200)))
color_2 = iter(plt.cm.rainbow(np.linspace(0,1,10)))
if key_word=="DOS":
    for DoS_1 in list_DoS_1:
        if decide_curves==0:
            ax.plot(w_list,DoS_1[:,1],linestyle="-",marker='o',markersize=3,linewidth=2,color=next(color_im))
        elif decide_curves==1:
            ax.plot(w_list,DoS_1[:,1],linestyle="-",marker='o',markersize=3,linewidth=2,color=next(color_2))
else:
    for self_array_im in list_self_array_im:
        #print("self_array_im = ", self_array_im, "\n\n")
        #ax.plot(w_list, self_array_re, linestyle="-", marker='o', markersize=3, linewidth=2, color='green', label='Re')
        if decide_curves==0:
            ax.plot(w_list, self_array_im[:,1], linestyle="-", marker='o', markersize=3, linewidth=2, color=next(color_im))
        elif decide_curves==1:
            ax.plot(w_list, self_array_im[:,1], linestyle="-", marker='o', markersize=3, linewidth=2, color=next(color_2))

if key_word=="Not_DOS":
    if tt == 400:
        plt.xlim([-0.01,41.01])
        plt.xticks(np.arange(0,42,1))#, rotation='vertical')
    elif tt == 200:
        plt.xlim([-0.01,21.01])
        plt.xticks(np.arange(0,22,1))
    else:
        raise ValueError('Should only have a grid of 200 or 400 Matsubara frequencies')
elif key_word=="DOS":
    plt.xlim([-10.01,10.01])
    plt.xticks(np.arange(-10,10,1))

if boool and key_word=="Not_DOS":
    plt.title(r"Self-energy $\omega$-dependence at different dopings for $\beta = {0:2.2f}$".format(beta)) # at  " r"$\mu=%.2f$"%mu[0][0])
elif boool and key_word=="DOS":
    plt.title(r"DoS at different dopings for $\beta = {0:2.2f}$".format(beta))    
else:
    plt.title(r"Green's function $\omega$-dependence at different dopings for $\beta = {0:2.2f}$".format(beta))

##############################################################
ax2 = plt.subplot(212) #ELIMINATE THIS SUBPLOT TO HAVE ONLY ONE GRAPH

xlabel = (r"$i\omega$")

if boool and  key_word=="Not_DOS":
    ylabel = (r"Re$\Sigma_{00}$"r"$\left(\omega\right)$")
elif boool and key_word=="DOS":
    ylabel = (r"DOS")
else:
    ylabel = (r"ReG$_{00}$"r"$\left(\omega\right)$")

plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

dop_list_re_list = list(sorted(zip(dop_range,list_re_list), key=lambda x : x[0]))
dop_list_re_winf = list(sorted(zip(dop_range,list_re_winf),key=lambda x: x[0]))
dop_list_DoS_2_list = list(sorted(zip(dop_range,list_DoS_2_list), key=lambda x: x[0]))
if key_word=="DOS":
    list_DoS_2 = []
    if decide_curves==0:
        eq_spaced_dop = np.arange(0,len(list_DoS_2_list),num_curves)
        for dop in eq_spaced_dop:
            list_DoS_2.append(dop_list_DoS_2_list[dop])
    elif decide_curves==1:
        for dop in list_dop_curves:
            list_DoS_2.append(dop_list_DoS_2_list[dop])
else:
    list_self_array_re = []
    list_self_array_re_winf = []
    if decide_curves==0:
        eq_spaced_dop = np.arange(0,len(list_re_list),num_curves)
        for dop in eq_spaced_dop:
            list_self_array_re.append(dop_list_re_list[dop]) ###Not dop_list_im_list[dop][1] to include doping in list_self_array_re
            list_self_array_re_winf.append(dop_list_re_winf[dop])
    elif decide_curves==1:
        for dop in list_dop_curves:
            list_self_array_re.append(dop_list_re_list[dop])  ###Not dop_list_im_list[dop][1] to include doping in list_self_array_re
            list_self_array_re_winf.append(dop_list_re_winf[dop])

#print("list_DoS_2 : ", list_DoS_2)
#print("list_self_array_re = ", list_self_array_re, "\n")
list_self_array_re_winf = np.asarray(list_self_array_re_winf)
print("list_self_array_re_winf",list_self_array_re_winf)

box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

if key_word=="DOS":
    for dop,DoS_2 in list_DoS_2:
        if decide_curves==0:
            ax2.plot(w_list,DoS_2[:,1],linestyle="-",marker='o',markersize=3,linewidth=2,color=next(color_re),label=r"$\delta=%.5f$"%(dop))
        elif decide_curves==1:
            ax2.plot(w_list,DoS_2[:,1],linestyle="-",marker='o',markersize=3,linewidth=2,color=next(color_2),label=r"$\delta=%.5f$"%(dop))
else:
    for dop,self_array_re in list_self_array_re:
        if decide_curves==0:
            ax2.plot(w_list, self_array_re[:,1], linestyle="-", marker='o', markersize=3, linewidth=2, color=next(color_re),label=r"$\delta=%.5f$"%(dop))#, label="Re") 
        elif decide_curves==1:
            ax2.plot(w_list, self_array_re[:,1], linestyle="-", marker='o', markersize=3, linewidth=2, color=next(color_2),label=r"$\delta=%.5f$"%(dop))

#plt.text(0.8, 0.8, r"$\delta=%.4f$"%(-1.0*(mu_list[MU]-1.0)), transform=ax2.transAxes, bbox=bbox_props)

ax2.legend(loc='upper center', bbox_to_anchor=(0.5,-0.15), fancybox=True, shadow=True, ncol=12)

if key_word=="Not_DOS":
    if tt == 400:
        plt.xlim([-0.01,41.01])
        plt.xticks(np.arange(0,42,1))#, rotation='vertical')
    elif tt == 200:
        plt.xlim([-0.01,21.01])
        plt.xticks(np.arange(0,22,1))
    else:
        raise ValueError('Should only have a grid of 200 or 400 Matsubara frequencies')
elif key_word=="DOS":
    plt.xlim([-10.01,10.01])
    plt.xticks(np.arange(-10,10,1))

fig2 = plt.figure(2)

axfig2 = plt.subplot(111)

ylabel = (r"Re$\Sigma_{00}$"r"$\left(\omega\to\infty\right)$")
xlabel = (r"doping $\delta$")
plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)
plt.title(r"Re$\Sigma\left(\omega\to\infty\right)$ at different dopings for $\beta = {0:2.2f}$".format(beta))
axfig2.plot(list_self_array_re_winf[:,0],list_self_array_re_winf[:,1]/U*1,linestyle="-",marker='o',markersize=3,linewidth=2,color='black')

plt.show()
fig.savefig(filename + ".pdf", format='pdf')
fig2.savefig(filename2 + ".pdf", format='pdf')


if key_word=="DOS":
    dop_list = [dop for dop,DoS in list_DoS_2]
    DOS_list = [DoS for dop,DoS in list_DoS_2]
else:
    dop_list = [dop for dop,self_array_re in list_self_array_re]
    self_array_re_list = [self_array_re for dop,self_array_re in list_self_array_re]

#dop_self_array_re_list = np.vstack((dop_list,self_array_re_list))
#print(len(w_list),len(DOS_list[3][:,0]))
################################################################
##Printing Imag and real parts for pade.py
################################################################
doping_indice = 30 ##Doping value for coexistence
#doping_indice = 68  ###Doping value for no coexistence
print("Doping is {0:1.6f}".format(dop_list[doping_indice]))
if boool and key_word=="Not_DOS":
    with open('self_energy.dat', 'w') as l:
        for i in range(len(w_list)):
            l.write("{0:6.5f}\t\t{1:6.5f}\t\t{2:6.5f}\t\t{3:3.2f}\n".format(self_array_re_list[doping_indice][i,0],self_array_re_list[doping_indice][i,1],list_self_array_im[doping_indice][i,1],dop_list[doping_indice]))
elif boool and key_word=="DOS":
    with open('DOS.dat', 'w') as l:
        for i in range(len(w_list)):
            l.write("{0:6.5f}\t\t{1:6.5f}\t\t{2:6.5f}\n".format(DOS_list[doping_indice][i,0],DOS_list[doping_indice][i,1],dop_list[doping_indice]))
else:
    with open('green_function.dat', 'w') as l:
        for i in range(len(w_list)):
            l.write("{0:6.5f}\t\t{1:6.5f}\t\t{2:6.5f}\t\t{3:3.2f}\n".format(self_array_re_list[doping_indice][i,0],self_array_re_list[doping_indice][i,1],list_self_array_im[doping_indice][i,1],dop_list[doping_indice]))

l.close()
################################################################
##Quasiparticle weight computation
################################################################

def Quasiparticle_weigth(w_0:float,self_im_w_0) -> float:
    im_sigma_w_0 = 1/4*np.trace(self_im_w_0)
    Z = 1/(1-im_sigma_w_0/w_0)
    return Z
if key_word=="Not_DOS":
    print("Quasiparticle weigth = ",Quasiparticle_weigth(w_list[0],np.asarray(np.imag(file_dict[doping_indice][0,0,:,:]))))
"""
anal_cont_data = np.loadtxt('analytic_continuation_self_energy.dat')
def eps(kx: float, ky: float, t: float = 1.0, tp: float = -0.3, tpp: float = 0.2) -> float:
    coskx = np.cos(kx)
    cosky = np.cos(ky)
    return(-2.0*t*(coskx + cosky) - 2.0*tp*(np.cos(kx+ky) + np.cos(kx-ky)) - 2.0*tpp*(np.cos(2.0*kx)+np.cos(2.0*ky)))

def E_rel(omega: float, kx: float, ky: float, ReSigma: float) -> float:
    return(omega - eps(kx,ky) - ReSigma)

grid = 200
kx = np.linspace(-np.pi,np.pi,grid)
ky = np.linspace(np.pi,-np.pi,grid)
k_grid = np.array(list(itertools.product(kx,ky)),dtype='float,float').reshape(grid,grid)

N = len(anal_cont_data[:,0])
peak_E = np.zeros((N,),dtype=float)

"""
