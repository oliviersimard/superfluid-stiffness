from SEvec_fig_util import Cdmft_data
import numpy as np
import re
from fnmatch import fnmatch
import warnings

param_dict = {}

#Parameters to specify in order to fetch the data to be plotted
params_val = ['./','*.npy',400,8,60,0.8,1.4]
params_key = ["paths","pattern","w_grid","U","beta","dop_min","dop_max"] #paths is for the current folder from where to search, 
                                                                         #pattern for the files to load,
                                                                         #w_grid is the grid of frequency 
                                                                         #U is the Hubbard interaction
                                                                         #beta is the fictitious T to smooth out the integration process over imaginary axis
                                                                         #dop_min is the minimal doping to be plotted
                                                                         #dop_max is the maximum doping to be plotted
params_dict = {key : val for key, val in zip(params_key,params_val)}

#The key folders into which the files with the desired pattern are fetched
key_folders = ["./NOCOEX/U8/SEvec/w_400/SEvec_b60_SC/SEvec_b60_SC","./COEX/U8/SEvec/w_400/SEvec_b60_SC_AFM/SEvec_b60_SC_AFM_s"]

#Instantiating cdmft
cdmft = Cdmft_data(**params_dict)

#Building the list of paths to the files containing the cdmft data
opt = "freq" #options implemented are 'HF' for plotting Hartree-Fock behavior of either G or SE with respect to doping and 'freq' for plotting either G or SE with respect to Matsubara frequency at different dopings

#path to the chemical-potential-vs-doping file
mu_list = './NOCOEX/U8/SEvec/w_400/SEvec_b60_SC/stiffness.dat'

#Generate the file path list to check if looking inside right folder
file_path_test = cdmft.gen_file_path(key_folders[0])

#Write down the name of the file for the plot to be stored
filename = "SEvec_NOCOEX_U8_beta_60_w_400_cluster_trace_0_2"

#Code line to generate the plots
cdmft.gen_plot(True, key_folders[0], 2, opt, mu_list, filename)
