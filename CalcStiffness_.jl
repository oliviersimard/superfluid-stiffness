using StiffSquare2.PeriodizeSC
using StiffSquare2.Stiffness


macro assertion(ex, text)
    :($ex ? nothing : error("Assertion failed: ", $(text)))
end

# Reading parameters from JSON file
paramsfile = "params__.json"
params = JSON.parsefile(paramsfile)
println(typeof(params))
# Defining own dictionnary for convenience
block_params = Dict{String,Any}(
	"path_to_files" => params["path_to_files"],
	"data_loop" => params["data_loop"],
	"pattern" => params["pattern"], # Pattern of the files to be loaded
	"beta" => params["beta"],
	"print_mu_dop" => params["Print_mu_dop"], # Parameter used to create stiffness.dat file (file used with SE_G_fig_producer.py)
	"w_discretization" => params["w_discretization"], # Matsubara frequency grid
	"AFM_SC_NOCOEX" => params["AFM_SC_NOCOEX"],
	"DOS" => params["DOS"],
	"cumulant" => params["cumulant"],
	"Periodization" => params["Periodization"],
	"fout_name" => params["fout_name"],
               "fout_name_DOS" => params["fout_name_DOS"],
               "abc" => params["abc"]
	)

# Check if appropriate input was entered
if block_params["cumulant"] in [0,1]; nothing else; throw(ErrorException("Cumulant takes in 0 or 1 as arguments")) end
if block_params["DOS"] in [0,1]; nothing else; throw(ErrorException("DOS takes in 0 or 1 as arguments")) end
if params["Save_gf_only"] in [0,1]; nothing else; throw(ErrorException("Save_gf_only takes in 0 or 1 as arguments")) end

# Check out if appriopriate inputs entered in data_loop_sparse
match_data_loop = match(r"^(.*?)(?=/)",block_params["data_loop"]).match
match_data_loop_arr = Array{String,1}([match_data_loop])
M_tol = 1e-3

# Decide the name of the file to be created
if block_params["DOS"] in [0,1]; nothing else; throw(ErrorException("DOS takes in either 0 or 1")) end
fout_name = block_params["DOS"] == 1 ? block_params["fout_name_DOS"] : block_params["fout_name"] # DOS can take 1 or 0

super_data_M, super_datap, list_of_files, datap = Stiffness.read_data_loop(block_params,block_params["print_mu_dop"],block_params["pattern"])
list_of_files_and_datap = datap != nothing ? collect(zip(list_of_files,datap)) : collect(zip(list_of_files,super_datap))

@assertion(length(list_of_files) != 0, "Length of list_of_files must differ from 0. Check if the path_to_files variable has been given the right path!")
@assertion(length(super_datap) != 0, "Length of super_datap must differ from 0. Check your file containing the chemical potential and the doping and make sure they are respectively marked \'mu\' and \'ave_mu\' in the header")

if Stiffness.verbose > 1
    println("Length list_of_files is : ", length(list_of_files),"\n")
    println("Length super_datap is: ", length(super_datap),"\n")
end

# Building list_modulevec defining the instance of modulevec
list_list_modulevec = Array{Array{Any,1},1}()
for (list_file,datap_mu) in list_of_files_and_datap
    length(list_file) != length(datap_mu[:,1]) && throw(File_Input_Error(:list_file))
    if Stiffness.verbose >=0
        println("Length list_file: ", length(list_file), "\tLength datap_mu[:,1]: ", length(datap_mu[:,1]))
    end
    list_modulevec = Stiffness.gen_modulevec_args(block_params,list_file,datap_mu[:,1],params)
    push!(list_list_modulevec, list_modulevec)
    if Stiffness.verbose >= 0
        println("typeof list_modulevec: ",typeof(list_modulevec))
    end
end

if length(list_list_modulevec) > 1 
    global super_list_modulevec
    super_list_modulevec = vcat(list_list_modulevec...)
else
    global super_list_modulevec
    super_list_modulevec = list_list_modulevec[1]
end

println(Stiffness.ts(super_list_modulevec[1]))
if Stiffness.verbose >= 0
    println("Size of sEvec_c in StiffnessArray: ", size(super_list_modulevec[1].data_), " | typeof StiffnessArray: ", typeof(super_list_modulevec[1]))
end

########################################### Saving Green's functions ###########################################
if params["Save_gf_only"] == 1
    save_gf(super_list_modulevec) # Green's functions are integrated over the Brillouin zone
    println("Saved the Green's functions (see function save_gf() for further informations)")  
    exit(0)
elseif params["Save_gf_only"] == 0
    nothing
end

########################################### Beginning of the superfluid stiffness calculations ###########################################
# Computing the superfluid stiffness
list_kIntegral_stiff, list_kIntegral_DOS, list_kIntegral_DOS_nocoex = Stiffness.calc_stiff_funct(super_data_M, super_list_modulevec, match_data_loop_arr, M_tol)

# Summing over the Matsubara frequencies
beta = block_params["beta"]
println("Summation in frequency beta = $(@sprintf("%.2f", beta))")
stiffness = (block_params["DOS"] == 1 && block_params["Periodization"] == 1) ? Array{Tuple{Int64,Float64},1}() : Array{Float64,1}()

if block_params["DOS"] == 1 && block_params["Periodization"] == 1
    if match_data_loop_arr[1] == "COEX"
        for j in 1:size(list_kIntegral_DOS)[1]
            stiffness_tmp = 0.0
            ind = list_kIntegral_DOS[j][1] #Tuples of indice to keep track of the order and of w vs nk(w)    
            nk_w = list_kIntegral_DOS[j][2]                                                                
            for ii in 1:size(nk_w)[1]
                stiffness_tmp+=nk_w[ii,2]                                                            
            end                                                                                                  
            push!(stiffness,(ind,1/2+(2.0/beta)*real(stiffness_tmp))) #Added 1/2 to avoid double counting                                
        end
    end
    for j in 1:size(list_kIntegral_DOS_nocoex)[1]                                                    
        stiffness_tmp = 0.0                                                         
        ind = list_kIntegral_DOS_nocoex[j][1]                                                 
        nk_w = list_kIntegral_DOS_nocoex[j][2]                                         
        for ii in 1:size(nk_w)[1]                                                     
            stiffness_tmp+=nk_w[ii,2]
        end
        push!(stiffness,(ind,1.0+(4.0/beta)*real(stiffness_tmp)))
    end
elseif block_params["DOS"] == 1 && block_params["Periodization"] == 0
    println("DOS == 1 and Periodization == 0; 4.0")
    for j in 1:size(list_kIntegral_DOS)[1]
        stiffness_tmp = 0.0
        for ii in 1:size(list_kIntegral_DOS[j])[1]
            stiffness_tmp+=list_kIntegral_DOS[j][ii,2]
        end
        push!(stiffness,(1.0 + (4.0/beta)*real(stiffness_tmp)))
        println(1.0 + (4.0/beta)*imag(stiffness_tmp))
    end

elseif block_params["DOS"] == 0
    for j in 1:size(list_kIntegral_stiff)[1]#-1
        stiffness_tmp = 0.0
        for ii in 1:size(list_kIntegral_stiff[j])[1]#-1
            stiffness_tmp+=list_kIntegral_stiff[j][ii, 2]
        end
        push!(stiffness,(1/beta)*real(stiffness_tmp))
        println((1/beta)*imag(stiffness_tmp))
    end
end

if block_params["Periodization"] == 1 && block_params["DOS"] == 1
    stiffness_sort = sort(stiffness, by=x->x[1])
    stiffness = [d[2] for d in stiffness_sort]
end

mu_dop_stiff = hcat(super_datap[1][1:length(stiffness),1],super_datap[1][1:length(stiffness),2],stiffness)

# Printing the file which contains the superfluid stiffness results 
println("Stiffness = ", mu_dop_stiff)
writedlm(fout_name, mu_dop_stiff, "\t\t")

println("Program exited")
