import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

path_to_file_ = "CUM_periodized_unitary_transformation.dat"
data_full = np.genfromtxt(path_to_file_,dtype=float,names=True,delimiter="\t\t")
list_n = data_full["n"]; list_tp = data_full["tp"]; 
list_mu = data_full["mu"]; list_beta = 1.0/data_full["beta"]
list_stiff = data_full["Stiff"]; list_SC_order = data_full["SC_Gap"]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(list_n,list_tp,list_beta,cmap=plt.cm.jet,c=list_stiff)
plt.title("Color plot of the superfluid stiffness")
max_tp = max(list_tp)+.1; min_tp = min(list_tp)-.1
max_mu = max(list_n)+.1; min_mu = min(list_n)-.1
max_beta = max(list_beta)+(1/200); min_beta = min(list_beta)-(1/200)

ax.set_xlim([min_mu, max_mu])
ax.set_ylim([min_tp, max_tp])
ax.set_yticks(np.arange(min_tp,max_tp,.1))
ax.set_zlim([min_beta, max_beta]) 
#ax.plot(filling, beta, 'k+', zdir='y', zs=max_U) # Useful if one wants to project data on plane
#ax.plot(U, beta, 'k+', zdir='x', zs=min_dop)
#ax.plot(filling, U, 'k+', zdir='z', zs=min_beta)

ax.set_xlabel(r"$n$") 
ax.set_ylabel(r"$tp$")
ax.set_zlabel(r"$T$")
#ax.plot(list_error_order_param, list_tp, list_beta, marker="_") #To eventually plot the error bars on SC order param. amplitude
fig.colorbar(p)
# Plotting colormap of errors
fig_e = plt.figure()
ax_e = fig_e.add_subplot(111, projection='3d')
ax_e.set_xlabel(r"$n$"); ax_e.set_ylabel(r"$tp$"); ax_e.set_zlabel(r"$T$")
err = ax_e.scatter(list_n,list_tp,list_beta,cmap=plt.cm.jet,c=list_SC_order)
plt.title("Color plot of the order parameter amplitude")
fig_e.colorbar(err) 
plt.show()
