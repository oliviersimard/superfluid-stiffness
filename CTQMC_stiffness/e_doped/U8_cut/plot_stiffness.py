import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from plot_utilities import container

############################# Important parameters for the plots ######################################
path_to_file_ = "PER_periodized_unitary_transformation_yy.dat"
beta_to_plot = 1.0/50
tp_arr = np.array([-0.1,-0.2,-0.3],dtype=float)
tp_cut = -0.3
data_full = np.genfromtxt(path_to_file_,dtype=float,names=True,delimiter="\t\t")
list_n = data_full["n"]; list_tp = data_full["tp"] 
list_mu = data_full["mu"]; list_beta = 1.0/data_full["beta"]
list_stiff = data_full["Stiff"]; list_SC_order = data_full["SC_Gap"]
list_zipped_params = list(zip(list_n,list_beta,list_tp,list_stiff,list_SC_order))

######################## Part of the code showing the 3D plot of superfluid stiffness #########################
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(list_n,list_tp,list_beta,cmap=plt.cm.jet,c=list_stiff)
plt.title("Color plot of the superfluid stiffness")
max_tp = max(list_tp)+.05; min_tp = min(list_tp)-.05
max_n = max(list_n)+.01; min_n = min(list_n)-.01
max_beta = max(list_beta)+(1/200); min_beta = min(list_beta)-(1/200)
ax.set_xlim([min_n, max_n])
ax.set_xticks(np.arange(min_n,max_n,.04))
ax.set_ylim([min_tp, max_tp])
ax.set_yticks(np.arange(min_tp,max_tp,.05))
ax.set_zlim([min_beta, max_beta]) 
ax.set_xlabel(r"$n$") 
ax.set_ylabel(r"$tp$")
ax.set_zlabel(r"$T$")
fig.colorbar(p)

# Plotting colormap of order parameter amplitude
fig_2 = plt.figure(2)
ax_2 = fig_2.add_subplot(111, projection='3d')
ax_2.set_ylim([min_tp, max_tp])
ax_2.set_yticks(np.arange(min_tp,max_tp,.05))
ax_2.set_xlabel(r"$n$"); ax_2.set_ylabel(r"$tp$"); ax_2.set_zlabel(r"$T$")
err = ax_2.scatter(list_n,list_tp,list_beta,cmap=plt.cm.jet,c=list_SC_order)
plt.title("Color plot of the order parameter amplitude")
fig_2.colorbar(err)

##################### Part of the code showing superfluid stiffness and order parameter as a function of n at any desired T ##########################
dict_vals = {}
max_dop_list = []
for tp_val in tp_arr:
    list_n_tmp = []; list_stiff_tmp = []; list_SC_order_tmp = []
    for params_tup in list_zipped_params:
        beta_tmp = params_tup[1]; tp_tmp = params_tup[2]
        if beta_tmp == beta_to_plot and tp_tmp == tp_val:
            list_n_tmp.append(params_tup[0])
            list_stiff_tmp.append(params_tup[3])
            list_SC_order_tmp.append(params_tup[4])
    dict_vals[tp_val] = (list_n_tmp,list_stiff_tmp,list_SC_order_tmp)
    max_dop_list.append(max(list_n_tmp))
# Plot for superfluid stiffness
fig_3 = plt.figure(3)
ax_3 = fig_3.add_subplot(111)
col = iter(plt.cm.rainbow(np.linspace(0,1,len(tp_arr))))
ax_3.set_xlabel(r"Particle density $n$", fontsize=15)
ax_3.set_ylabel(r"Superfluid stiffness $\rho_s$", fontsize=15)
plt.title(r"$\rho_s$ vs $n$ for all tp's at temperature $\beta$ = %.1f"%(1.0/beta_to_plot))
max_dop_list_stiff = []; max_stiff_list = []
for tp_val in dict_vals.keys():
    max_rho = max(dict_vals[tp_val][1])
    max_stiff_list.append(max_rho)
    max_doping = dict_vals[tp_val][0][dict_vals[tp_val][1].index(max_rho)]
    max_dop_list_stiff.append(max_doping)
    ax_3.plot(dict_vals[tp_val][0],dict_vals[tp_val][1],linestyle="-",marker='.',markersize=2,linewidth=1,c=next(col),label="tp = %.1f"%(tp_val))
    ax_3.annotate(r"max $\rho_s=%.3f$"%(max_rho), xy=(max_doping, max_rho), xytext=(max_doping-0.023, max_rho-0.01),
            arrowprops=dict(facecolor='black', shrink=0.05, width=0.05, headwidth=3, headlength=3),
            )

# Plot for order parameter
fig_4 = plt.figure(4)
ax_4 = fig_4.add_subplot(111)
ax_4.set_xlabel(r"Particle density $n$", fontsize=15)
ax_4.set_ylabel(r"Order parameter amplitude $\left|\phi\right|$", fontsize=15)
plt.title(r"$\left|\phi\right|$ vs $n$ for all tp's at temperature $\beta$ = %.1f"%(1.0/beta_to_plot))
max_dop_list_order = []; max_order_param_list = []
coll = iter(plt.cm.rainbow(np.linspace(0,1,len(tp_arr))))
for i,tp_val in enumerate(dict_vals.keys()):
    max_order_param = max(dict_vals[tp_val][2])
    max_order_param_list.append(max_order_param)
    max_doping = dict_vals[tp_val][0][dict_vals[tp_val][2].index(max_order_param)]
    max_dop_list_order.append(max_doping)
    ax_4.plot(dict_vals[tp_val][0],dict_vals[tp_val][2],linestyle="-",marker='.',markersize=2,linewidth=1,c=next(coll),label="tp = %.1f"%(tp_val))
    #ax_4.annotate(r"max $\left|\phi\right|=%.3f$"%(max_order_param), xy=(max_doping, max_order_param), xytext=(max_doping-0.03, max_order_param-0.01),
    #        arrowprops=dict(facecolor='black', shrink=0.05, width=0.05, headwidth=3, headlength=3),
    #        )

########## Plot showing the behavior of both order parameter maxima and superfluid stiffness maxima for each tp as a function of doping ###########

fig_5 = plt.figure(5)
ax_5 = fig_5.add_subplot(111)
ax_5.set_xlabel(r"Particle density $n$", fontsize=15)
ax_5.set_ylabel(r"$t^{\prime}$", fontsize=15)
plt.title(r"Maximum of $\left|\phi\right|$ and maximum of $\rho_s$ for each tp vs $n$ at temperature $\beta$ = %.1f"%(1.0/beta_to_plot))
ax_5.plot(max_dop_list_stiff,tp_arr,linestyle="-",marker='.',markersize=2,linewidth=1,color='red',label=r"max $\rho_s$")
ax_5.plot(max_dop_list_order,tp_arr,linestyle="-",marker='.',markersize=2,linewidth=1,color='blue',label=r"max $\left|\phi\right|$")
ax_5.plot(max_dop_list,tp_arr,linestyle="-",marker='.',markersize=2,linewidth=1,color='green',label=r"max $n$")

######### Plot showing difference between maximum doping and doping at which the superfluid stiffness is maximum, as a function of T #############

beta_diff = []
dict_vals_tp_cut = {}
for beta in list_beta:
    container(beta_diff,beta)
for beta in beta_diff:
    list_n_tmp = []; list_stiff_tmp = []; list_SC_order_tmp = []
    for params in list_zipped_params:
        beta_tmp = params[1]; tp_tmp = params[2]
        if beta_tmp == beta and tp_tmp == tp_cut:
            list_n_tmp.append(params[0])
            list_stiff_tmp.append(params[3])
            list_SC_order_tmp.append(params[4])
    try:
        max_n = max(list_n_tmp) 
        max_stiff_dop = list_n_tmp[list_stiff_tmp.index(max(list_stiff_tmp))]
        max_SC_dop = list_n_tmp[list_SC_order_tmp.index(max(list_SC_order_tmp))]
    except:
        max_n = max_stiff_dop = max_SC_dop = 0.0
    if max_n != 0.0 and max_stiff_dop != 0.0 and max_SC_dop != 0.0: 
        dict_vals_tp_cut[beta] = (max_n,max_stiff_dop,max_SC_dop)

fig_6 = plt.figure()
ax_6 = fig_6.add_subplot(111)
plt.title(r"$\delta_c-\delta_{max}(\rho_s)$ vs $\beta$ for $tp$ = %.2f"%(tp_cut))
plt.xlabel(r"Inverse temperature $\beta$",fontsize=15)
for beta in beta_diff:
    try:
        ax_6.plot(beta,dict_vals_tp_cut[beta][0]-dict_vals_tp_cut[beta][1],marker='.',markersize=3,color='red')
    except:
        print("Empty lists encountered!")

plt.legend()

plt.show()
