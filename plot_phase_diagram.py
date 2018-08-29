import numpy as np
import os 
from fnmatch import fnmatch
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import argparse

License = """
        plot_phase_diagram: A program nobody cares about
        copyleft (C, but facing left) 2018  Olivier Simard <olivier.simard2@usherbrooke.ca> 

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""


global _verbose
_verbose = 0
_NcN = 8
_plotting_option = 'pg'

#-------------------------------------------------------------------------------Setting the command line parser------------------------------------------------------------------------------------#
####################################################################################################################################################################################################
def argparser() -> tuple:
    """Simple command line argument parser for more options"""
    parser = argparse.ArgumentParser(description='Plot the phase diagram or the self-energies/Green\'s functions. The program saves the data in binary files for later usage.')
    parser.add_argument('-p', '--plotting',dest='plot',help='Option to plot either the phase diagram or the self-energies (SE)/Green\'s functions (GE). Default is plotting the phase diagram. Choose \'s\' to plot only SE, \'g\' to plot only  GE or \'all\' to plot all.',default = argparse.SUPPRESS)
    parser.add_argument('-f', '--filename',dest='filename',help='Option to specify the file for the SE and/or the GE to be plotted. Must follow \'-p\'.',default = argparse.SUPPRESS)
    parser.add_argument('-v', '--verbose',dest='verbose',help='Option to set the verbose. Default value is 0.')
    parser.add_argument('-a','--about',dest='about',help='About this program',nargs='?',const=True,default=False)
    
    args = parser.parse_args()

    if args.about:
        print(License)
        return None, None, None
    else:
        try:
            return args.plot, args.verbose, args.filename
        except:
            return args.plot, args.verbose, None

#---------------------------------------------------------------------------------Getting the converged data---------------------------------------------------------------------------------------#
####################################################################################################################################################################################################

def convert_converged_green(file_converged_g: dict) -> None:
    """Function to save converted Green's function (also self-energy) ready for superfluid stiffness calculations"""
    data = np.genfromtxt(file_converged_g,names=True)
    len_iwn = len(data["w_matsubara"])
    Gvec = np.zeros((len_iwn,_NcN,_NcN),dtype=complex)

    for w_ind in range(len_iwn):
        g11 = data["g11_Re"][w_ind]+1j*data["g11_Im"][w_ind]; g12 = data["g12_Re"][w_ind]+1j*data["g12_Im"][w_ind];
        g13 = data["g13_Re"][w_ind]+1j*data["g13_Im"][w_ind]; g12_N = data["g12_nambu_Re"][w_ind]+1j*data["g12_nambu_Im"][w_ind]
        U_trans = np.array([[1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 0, 1],
                            [0, 0, 1, 0]],dtype=int) # Matrix transformation to recover qcm way to label cluster sites
        G_up = np.array([[g11, g12, g13, g12],
                         [g12, g11, g12, g13],
	                 [g13, g12, g11, g12],
	                 [g12, g13, g12, g11]],dtype=complex)
        G_up = U_trans*G_up*U_trans
        G_down = -G_up.copy()
        F = np.array([[0, g12_N, 0, g12_N],
                      [g12_N, 0, g12_N,0],
                      [0, g12_N, 0, g12_N],
                      [g12_N, 0, g12_N, 0]],dtype=complex)
        F = U_trans*F*U_trans
        F_dag = F.copy()
        G_up_F = np.concatenate((G_up,F),axis=1); F_dag_G_down = np.concatenate((F_dag,G_down),axis=1)
        Gvec[w_ind] = np.concatenate((G_up_F,F_dag_G_down),axis=0)

    fi = os.path.basename(os.path.realpath(file_converged_g)).split(".")[0]
    path = os.path.dirname(os.path.realpath(file_converged_g))
    np.savetxt(path+"/w.dat",data["w_matsubara"],fmt='%4.8f')
    full_path = path+"/"+fi
    #if os.path.exists(full_path+".npy"):
    #    os.remove(full_path+".npy")
    if _verbose > 0:
        print("full path: ", full_path)
    np.save(full_path, Gvec)
    return None


def get_dict(list_of_files: str, list_of_subdirs: str) -> dict:
    """Function generating dictionnary containing subdirs and associated Green's function files"""
    dir_vals = {}
    list_subdirs_and_files = []
    for subdirs in list_of_subdirs:
        for files in list_of_files:
            if subdirs in files:
                if _verbose > 0:
                    print(subdirs, "inside", files)
                list_subdirs_and_files.append(files)    
        dir_vals[subdirs] = list_subdirs_and_files
        list_subdirs_and_files = []
    return dir_vals

def get_converged_Green(list_of_files: str, list_of_subdirs: str) -> dict:
    """Function returning converged Green's function (last iteration) for each subdirs"""
    dic = get_dict(list_of_files, list_of_subdirs)
    list_int = []
    dict_converged_file = {}
    for keys in dic.keys():
        max_val = 0
        if _verbose > 0:
            print("keys: ", keys, "\n")
        for path_files in dic[keys]:
            if _verbose > 0:
                print("path_files: ", path_files, "\n")
            curr_int = int(re.findall(r'(?<=\w)(\d+\.\d+|\d+)',path_files)[-1])
            list_int.append(curr_int)
            if curr_int > max_val:
                max_val = curr_int
        dict_converged_file[keys] = [path_files for path_files in dic[keys] if "greenR"+str(max_val) in path_files][0]
        assert max_val !=0, "Check if folder is empty or if the files are properly named! Problem occured in %s"%(keys)
        if _verbose > 0:
            print("list of int: ", list_int, "\n")
        list_int = []
    return dict_converged_file


def get_SC_order_parameter(list_of_files: str, list_of_subdirs:str) -> list:
    """Function returning the relevant physical parameters to plot the phase diagram"""
    relevant_param_list = []
    dict_converged_file = get_converged_Green(list_of_files, list_of_subdirs)
    for keys in dict_converged_file.keys():
        info_from_keys = re.findall(r'(?<=\w)(\d+\.\d+|\d+)',keys)
        filling = float(info_from_keys[1]); U = float(info_from_keys[2]); beta = float(info_from_keys[3])
        loading_Nambu_G = np.loadtxt(dict_converged_file[keys],dtype=float,comments='#',skiprows=1,usecols=(0,7))
        self_dirname = os.path.dirname(dict_converged_file[keys])
        self_basename = os.path.basename(dict_converged_file[keys])
        num_self = re.search(r"\d+",self_basename).group()
        string_for_self = self_dirname+"/"+"selfR"+num_self+".dat"
        convert_converged_green(dict_converged_file[keys]) # Converting Green's functions
        convert_converged_green(string_for_self) # Converting self-energies
        sum_G_Nambu = sum(loading_Nambu_G[:,1])
        abs_SC_order_param = (1/beta)*abs(sum_G_Nambu)
        print(abs_SC_order_param,filling,U,1/beta)
        relevant_param_list.append([abs_SC_order_param,filling,U,beta])
    return relevant_param_list

def plot_phase_diagram(list_of_files: str,list_of_subdirs: str) -> None:
    relevant_params = get_SC_order_parameter(list_of_files,list_of_subdirs)
    abs_SC_order_param = np.array([el[0] for el in relevant_params],dtype=float); filling = np.array([el[1] for el in relevant_params],dtype=float)
    U = np.array([el[2] for el in relevant_params],dtype=float); beta = np.array([1/el[3] for el in relevant_params],dtype=float)
    print("Order parameter: ", abs_SC_order_param, "\n", "Filling: ", filling, "\n", "U: ", U, "\n", "Temperature: ", beta)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    p = ax.scatter(filling,U,beta,cmap=cm.jet,c=abs_SC_order_param)

    max_U = max(U)+1; min_U = min(U)-1
    max_dop = max(filling)+0.01; min_dop = min(filling)-0.01
    max_beta = max(beta)+(1/200); min_beta = min(beta)-(1/200)

    ax.set_xlim([min_dop, max_dop])
    ax.set_ylim([min_U, max_U])
    ax.set_yticks(np.arange(min_U,max_U,1))
    ax.set_zlim([min_beta, max_beta]) 
    ax.plot(filling, beta, 'k+', zdir='y', zs=max_U)
    ax.plot(U, beta, 'k+', zdir='x', zs=min_dop)
    ax.plot(filling, U, 'k+', zdir='z', zs=min_beta)

    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$U$")
    ax.set_zlabel(r"$T$")

    fig.colorbar(p) 
    
    plt.show()
    return None

def plot_SE_and_or_GE(list_of_files: str,list_of_subdirs: str,particular_file: str, args: str) -> None:
    fig = plt.figure()
    axIm = fig.add_subplot(211)
    x_lim = 10
    if len(particular_file) != 0:
        params_arr = re.findall(r'(?<=\w)(\d+\.\d+|\d+)',particular_file)
        beta = float(params_arr[2]); n = float(params_arr[0]); U = float(params_arr[1])
        file_to_be_loaded = ''
        dict_converged_file = get_converged_Green(list_of_files, list_of_subdirs)
        tmp_file_str = dict_converged_file[particular_file]
        f_dirname = os.path.dirname(tmp_file_str)
        w_list = f_dirname+"/w.dat"
        if args == 'g':
            f_without_ext = os.path.basename(os.path.realpath(tmp_file_str)).split(".")[0]
            file_to_be_loaded = f_dirname+"/"+f_without_ext+".npy"
        elif args == 's':
            s_basename = os.path.basename(tmp_file_str)
            num_self = re.search(r"\d+",s_basename).group()
            string_to_self = f_dirname+"/"+"selfR"+num_self+".npy"
            file_to_be_loaded = string_to_self
        else:
            raise ValueError
        data_file_loaded = np.load(file_to_be_loaded)
        w_list = np.loadtxt(f_dirname+"/w.dat")
        SorG_el_im = [el[0,0].imag for el in data_file_loaded]
        SorG_el_re = [el[0,0].real for el in data_file_loaded]
        if args == 's':
            plt.title(r"Self-energy $\omega$-dependence for $\beta$ = %.1f, U = %.1f and n = %.1f"%(beta,U,n))
        elif args == 'g':
            plt.title(r"Green's function $\omega$-dependence for $\beta$ = %.1f, U = %.1f and n = %.1f"%(beta,U,n))

        axIm.plot(w_list,SorG_el_im,linestyle="-",marker='o',markersize=3,linewidth=2,color='red')
        plt.ylabel("Im", fontsize=15)
        if args == 's':
            plt.xlim([-0.01,x_lim])
            plt.ylim([-2.5,0])

        axRe = fig.add_subplot(212)
        axRe.plot(w_list,SorG_el_re,linestyle="-",marker='o',markersize=3,linewidth=2,color='green')
        plt.ylabel("Re", fontsize=15)
        plt.xlabel(r"Matsubara frequency $i\omega_n$",fontsize=15)
        if args == 's':
            plt.xlim([-0.01,x_lim])
            plt.ylim([4,2])
        plt.show()
    else:
        relevant_params = get_SC_order_parameter(list_of_files,list_of_subdirs)
    return None


#___________________________________________________________________________________MAIN________________________________________________________________________________________________#
#########################################################################################################################################################################################
if __name__ == '__main__':
    args = argparser()
    print("ARGS:", args)
    Filename = ''
    if not (args[2] is None):
        Filename = args[2]
    if not (args[1] is None):
        _verbose = int(args[1])
    if not (args[0] is None):
        _plotting_option = args[0]
    else:
        pass
    pattern_subdirs = "plaquette*from_self*"
    pattern_file = "greenR*.dat"
    re_p = re.compile("(?<=R)[0-9]+")
    list_of_files = []; list_of_subdirs = []; list_of_params = []

    for path, subdirs, files in os.walk(os.getcwd()):
        for name in subdirs:
            print("name of subdirs: ", name)
            if fnmatch(name,pattern_subdirs):
                list_of_subdirs.append(name)
                list_of_params.append(re.findall(r'[-+]?\d*\.\d+|\d+',name))
        for filename in files:
            if fnmatch(filename,pattern_file):
                list_of_files.append(os.path.join(path,filename))

    #m = re.findall(r'[-+]?\d*\.\d+|\d+', list_of_files[0]) # m and p options give the same result
    #n = re.findall(r'(?<=R)[0-9]+', list_of_files[12])
    #o = re.findall(r'[\d+\.\d+|\d+](?=(_.*_))', list_of_files[0]) ## Lookahead expression
    #p = re.findall(r'(?<=\w)(\d+\.\d+|\d+)', list_of_files[0]) ## Lookbehind expression to look for numbers (whether decimals or not) following a letter

    d = get_dict(list_of_files,list_of_subdirs)
    dd = get_converged_Green(list_of_files,list_of_subdirs)
    
    if _verbose > 0:
        for k in d.keys():
            print(k, d[k], "\n")

#---------------------------------------------------------------------------------Plotting the data---------------------------------------------------------------------------------------#
###########################################################################################################################################################################################
    if _plotting_option == 'pg':
        plot_phase_diagram(list_of_files,list_of_subdirs)
    elif _plotting_option == 's' or 'g' or 'all':
        try:
            plot_SE_and_or_GE(list_of_files, list_of_subdirs, Filename, _plotting_option)
        except ValueError as e:
            print("Invalid input in command line after \'-p\'.")
            raise e
    else:
        raise 

    exit()

