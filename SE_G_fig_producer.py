from SEvec_fig_util import Cdmft_data, Correction

######################Extract Self-Energy or Green's function#############################

#Parameters to specify in order to fetch the data to be plotted
params_val = ['./','*.npy',2000,8,120,0.8,1.3,0]
params_key = ["paths","pattern","w_grid","U","beta","dop_min","dop_max","verbose"] #paths is for the current folder from where to search, 
                                                                         #pattern for the files to load,
                                                                         #w_grid is the grid of frequency 
                                                                         #U is the Hubbard interaction
                                                                         #beta is the fictitious T to smooth out the integration process over imaginary axis
                                                                         #dop_min is the minimal doping to be plotted
                                                                         #dop_max is the maximum doping to be plotted
params_dict = {key : val for key, val in zip(params_key,params_val)}

#The key folders into which the files with the desired pattern are fetched
key_folders = ["./COEX/U8/SEvec/w_2000/SEvec_b120_SC_AFM","./COEX/U8/SEvec/w_400/SEvec_b60_SC_AFM/SEvec_b60_SC_AFM_s"]

#Instantiating cdmft
cdmft = Cdmft_data(**params_dict)

#Building the list of paths to the files containing the cdmft data
opt = "HF" #options implemented are 'HF' for plotting Hartree-Fock behavior of either G or SE with respect to doping and 'freq' for plotting either G or SE with respect to Matsubara frequency at different dopings

#path to the chemical-potential-vs-doping file
mu_list = './COEX/U8/SEvec/w_2000/SEvec_b120_SC_AFM/stiffness.dat'

#Generate the file path list to check if looking inside right folder
file_path_test = cdmft.gen_file_path(key_folders[0])

#Write down the name of the file for the plot to be stored
filename = "SEvec_COEX_U8_beta_120_w_2000_cluster_trace_0_8_"

#Code line to generate the plots
tr_indices = ":".join(("0","8"))
cdmft.gen_plot(False, key_folders[0], tr_indices, opt, mu_list, filename) #True if the folders holding the data are sparse.
#cdmft.get_quasi_weight(False, key_folders[0], tr_indices, opt, mu_list)
###########################Correction part of the code######################################

#Files to load for the correction operation
tol = 1e-3
file_COEX_data = "stiffness/stiffness_U655/COEX/Per/w_2000/tperp/stiffness_b90_w_2000_coex_int_K_per_tperp_U655_1e-3.dat"
file_NOCOEX_data = "stiffness/stiffness_U655/NOCOEX/Per/w_2000/tperp/stiffness_b90_w_2000_nocoex_int_K_per_tperp_U655_1e-3.dat"
file_NOCOEX_AFM_per = "stiffness/stiffness_U655/NOCOEX/Per/w_2000/tperp/AFM_SC_1/stiffness_b90_w_2000_nocoex_AFM_SC_1_int_K_per_tperp_U655_1e-3.dat"
var_para_info_file = "COEX/U655/Loop_COEX_tmp.dat"
filename_comparison = "stiffness_b90_w_2000_coex_int_K_AFM_per_1e-3_tperp_U655.dat.corr" #Filename of the figure showing the corrected coex data

#Loading files' info
correct_info = Correction(tol, file_COEX_data, file_NOCOEX_data, file_NOCOEX_AFM_per, var_para_info_file)
correct_info.gen_plot_comparison("polyfit",6,filename_comparison)

#print(type(correct_info), Correction.__dict__)
