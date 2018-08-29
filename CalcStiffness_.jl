using JSON
using NPZ
using Glob
using Memento
using StiffSquare2.PeriodizeSC

###########################################Functions and structures###########################################
global pwd_ = pwd()
# Set verbose here
global verbose = 0
global logger

# Verify if Julia.log exists. Removes it if case.
if isfile(pwd_*"/Julia.log")
    rm(pwd_*"/Julia.log")
end

# Memento logger definition
print_with_color(:green, "Setting parameters for log file Julia.log\n")

logger = getlogger(current_module())
push!(logger, DefaultHandler("Julia.log"))
Memento.config!("debug";fmt="[{level}|{name}]:{msg}")

macro assertion(ex, text)
    :($ex ? nothing : error("Assertion failed: ", $(text)))
end

struct File_Input_Error <: Exception
    var::Symbol
end

# Exception constructor
Base.showerror(io::IO, e::File_Input_Error) = println(io, e.var, " Keys are not balanced (See function read_data_loop)")

# Static type object containing all important parameters
mutable struct StiffnessArray{T, N, S <: AbstractArray, Y <: Associative{String,Any}} <: AbstractArray{T, N}
    # Contains the data array
    data_::S
    # Can have maximum three non-zero hopping parameters
    hopping_val_::Base.RefValue{Tuple{T, T, T, T}}
    # Contains informations about the calculations to do
    vals_::Y
    #addContent::Function
    #setHeader::Function
    #getHeaders::Function

    #function StiffnessArray{T,N,S,Y}() where {T,N,S,Y}
    #    this = new()
    #    this.vals_ = Dict{String,Any}()

        #this.addContent = function (append::String)
        #    this.data = this.data * append
        #end

    #    this.setHeader = function(header::String, value::Any)
    #        this.vals_[header] = value
    #    end

        #this.getHeaders = function ()
        #    headers = ""
        #    for (header, value) in this.vals_
        #        vals_ = vals_ * header
        #	 if length(value) > 0
        #	     vals_ = vals_ * ": " * value
        #        end
        #	 vals_ = vals_ * "\n"
        #    end
        #    return vals_
        #end
    #    return this
    #end
end

# Outer constructor:
StiffnessArray{T,N}(data_::AbstractArray{Complex{T},N}, hopping_val_::Tuple{T,T,T,T}, vals_::Associative{String,Any}) = StiffnessArray{T,N,typeof(data_),typeof(vals_)}(data_, Ref(hopping_val_), vals_)
ts(t::StiffnessArray) = t.hopping_val_[] # Getter
Base.size(A::StiffnessArray) = size(A.data_)
Base.IndexStyle{T<:StiffnessArray}(::Type{T}) = Base.IndexLinear()

function Base.getindex(A::StiffnessArray, I...)
    checkbounds(Bool, A.data_, I...) && return A.data_[I...]
    I... < 1 ? throw(ErrorException("Out of bounds (lower bound)")) : throw(ErrorException("Out of bounds (upper bound)"))
end

"""
Function used to build main dictionnary containing file related information

#Argument(s):

- block_params: Dict-valued argument produced form JSON parsing.
- print_mu_dop: Binary-valued argument (1 or 0) to decide to print (or not) a file
  describing chemical potential vs doping
- pattern: String-valued argument representing the extension common to the python binary files (*.npy)

#Returns:

- super_datap: 2d array containing chemical potential on first column and doping on second one (Array{Float64,2})
- stiffness.dat: If print_mu_dop = 1, returns a file and exits the program
- list_of_files:
"""
function read_data_loop(block_params::Dict{String,Any}, print_mu_dop::Int64, pattern::String)
    dict_data_file = Dict{String,Any}()
    if block_params["not_sparse"] == 0 # <----------------------------------------------------
        len_block_data = length(block_params["data_loop_sparse"])
        for l in 1:len_block_data
            dict_data_file["data_file_s$(l)"] = open(readdlm, block_params["data_loop_sparse"][l])

            dict_data_file["data_file_header_s$(l)"] = dict_data_file["data_file_s$(l)"][1,:]

            dict_data_file["indmu_s$(l)"] = find(x->x=="mu",dict_data_file["data_file_header_s$(l)"])
            dict_data_file["inddop_s$(l)"] = find(x->x=="ave_mu",dict_data_file["data_file_header_s$(l)"])
            dict_data_file["indM_s$(l)"] = find(x->x=="ave_M",dict_data_file["data_file_header_s$(l)"])

            dict_data_file["data_file_mu_s$(l)"] = dict_data_file["data_file_s$(l)"][:,dict_data_file["indmu_s$(l)"]]
            dict_data_file["data_file_dop_s$(l)"] = dict_data_file["data_file_s$(l)"][:,dict_data_file["inddop_s$(l)"]]
            dict_data_file["data_file_M_s$(l)"] = dict_data_file["data_file_s$(l)"][:,dict_data_file["indM_s$(l)"]]

            dict_data_file["data_file_mu_s$(l)"] = convert(Array{Float64},filter(x->x!="mu",dict_data_file["data_file_mu_s$(l)"]))
            dict_data_file["data_file_dop_s$(l)"] = filter(x->x!="ave_mu",dict_data_file["data_file_dop_s$(l)"])
            dict_data_file["data_file_M_s$(l)"] = data_file_M_s1 = filter(x->x!="ave_M",dict_data_file["data_file_M_s$(l)"])

            dict_data_file["datap_h_s$(l)"] = hcat(dict_data_file["data_file_mu_s$(l)"],dict_data_file["data_file_dop_s$(l)"])

            dict_data_file["list_of_files_s$(l)"] = glob(string(block_params["path_to_files_sparse"][l],pattern))
        end

        # Check for the keys with name "datap_h_s" in dict_data_file
        extracted_keys_M = [dict_data_file[l] for l in keys(dict_data_file) if ismatch(r"data[\w]*_M_s",l)]
        extracted_keys_datap_h = [dict_data_file[l] for l in keys(dict_data_file) if ismatch(r"datap[\w]*_s",l)]
        extracted_keys_list_files = [dict_data_file[l] for l in keys(dict_data_file) if ismatch(r"list[\w]*_s",l)]
        if verbose >= 0
            # The order should be *_s2 and *_s1 to refer to previous code where *_s2 = * and *_s1 = *_s
            println("Order of datap_h (for vcat): ", [l for l in keys(dict_data_file) if ismatch(r"datap[\w]*_s",l)])
            println("Order of list_of_files (for vcat): ", [l for l in keys(dict_data_file) if ismatch(r"list[\w]*_s",l)])
            println("Order of data_M_s: ", [l for l in keys(dict_data_file) if ismatch(r"data[\w]*_M_s",l)])
        end
        len_datap_h = length(extracted_keys_datap_h)
        println(len_datap_h)
        len_list_files = length(extracted_keys_list_files)
        println(len_list_files)
        len_data_M = length(extracted_keys_M)
        println(len_data_M)
        # Assertion that len_datap_h == len_list_files
        (len_datap_h != len_list_files || len_datap_h != len_data_M) && throw(File_Input_Error(:len_datap_h))
        if verbose >= 0
            println("Length of files extracted: ", len_list_files, "\ttype: ", typeof(extracted_keys_list_files))
            println("Number of different (sparse) mu vs doping tables: ", len_datap_h, "\ttype: ", typeof(extracted_keys_datap_h))
            println("Length of data_M: ", len_data_M, "\ttype: ", typeof(extracted_keys_M))
            println("Size(s) of files extracted: ", size.(extracted_keys_list_files))
            println("Size(s) of mu-vs-doping tables extracted: ", size.(extracted_keys_datap_h))
            println("Size(s) of set(s) of data_M: ", size.(extracted_keys_M))
        end
        super_datap = vcat(extracted_keys_datap_h...)
        super_data_M = vcat(extracted_keys_M...)
        super_datap = convert(Array{Float64,2},super_datap)
        println(super_data_M)
        list_super_datap = Array{Array{Float64,2},1}()
        push!(list_super_datap,super_datap)
        # Exits the program if print_mu_dop = 1
        if print_mu_dop == 1
            writedlm("stiffness_.dat", super_datap, "\t\t")
            print_with_color(:red, "Printed stiffness_.dat for later use (SE_G_fig_producer.py)\n")
            exit(0)
        else
            nothing
        end
        return super_data_M, list_super_datap, extracted_keys_list_files, extracted_keys_datap_h

    elseif block_params["not_sparse"] == 1 # <----------------------------------------------------
        data_file = open(readdlm, block_params["data_loop"])
        data_file_header = data_file[1,:]

        indmu = find(x->x=="mu",data_file_header)
        inddop = find(x->x=="ave_mu",data_file_header)
        indM = find(x->x=="ave_M",data_file_header)

        data_file_mu = data_file[:,indmu]
        data_file_dop = data_file[:,inddop]
        data_file_M = data_file[:,indM]
        #println(length(data_file_mu),"\n",length(data_file_dop),"\n",length(data_file_M))

        data_file_mu = filter(x->x!="mu",data_file_mu)
        data_file_dop = filter(x->x!="ave_mu",data_file_dop)
        data_file_M = filter(x->x!="ave_M",data_file_M)
        #println(length(data_file_mu),"\n",length(data_file_dop),"\n",length(data_file_M))

        datap_h = hcat(data_file_mu,data_file_dop)
        datap_h = convert(Array{Float64,2},datap_h)
        data_file_M = convert(Array{Float64,1},data_file_M)
        #println(length(datap_h),"\n",length(data_file_M))
        # Having an Array{Array{String,1},1} and an Array{Array{Float64,2},1} object to facilitate the task
        list_of_files = Array{Array{String,1},1}()
        list_datap_h = Array{Array{Float64,2},1}()
        push!(list_of_files,glob(string(block_params["path_to_files"],pattern)))
        push!(list_datap_h,datap_h)
        # Exits the program if print_mu_dop = 1
        if print_mu_dop == 1
            writedlm("stiffness.dat", super_datap, "\t\t")
            print_with_color(:red, "Printed stiffness.dat for later use (SE_G_fig_producer.py)\n")
            exit(0)
        else
            nothing
        end
        #println(length(data_file_M), "\n", length(list_datap_h), "\n", length(list_of_files))

        return data_file_M, list_datap_h, list_of_files, nothing
    end
end


"""
Function extracting the last number of the binary filename

#Argument(s):

- list_of_files: 1d array-valued argument containing the filenames

#Returns:

- list_num: 1d array-valued output containing the last number of the filenames
"""
function gen_file_num(list_of_files::Array{String,1})
    list_num = Array{Int64,1}(length(list_of_files))
    for (i,files) in enumerate(list_of_files)
        m = [x for x in eachmatch(r"\d+",files)][end]
        list_num[i] = parse(Int64,m.match)
    end
    if verbose >= 0
        println("Length list_num: ", length(list_num))
    end
    return list_num
end


"""
Function producing the arguments defining modulevec (Initiate the instance)

#Argument(s):

- list_of_files: 1d array-valued argument containing the filenames

#Returns:

- list__mu: 1d array containing chemical potential
- list_SEvec_c: 1d array containing the self-energies taken from qcm output files
"""
function gen_modulevec_args(block_params::Dict{String,Any},list_of_files::Array{String,1},data_mu::Array{Float64,1},params::Dict{String,Any})
    list_modulevec = Array{Any,1}() #; list_t_mu = Array{Tuple{Float64,Float64,Float64,Float64},1}()
    list_num = gen_file_num(list_of_files)
    zip_num_files = collect(zip(list_num,list_of_files))
    sort!(zip_num_files)
    for (num,files) in zip_num_files
        sEvec_c = npzread(files)
        if verbose >= 0
            println(num, " : ", (params["t"],params["tp"],params["tpp"],data_mu[num+1]))
        end
        push!(list_modulevec, StiffnessArray(sEvec_c, (params["t"],params["tp"],params["tpp"],data_mu[num+1]), block_params))
    end

    return list_modulevec
end

"""
Function calculating the superfluid stiffness calling the module StiffSquare2.PeriodizeSC

#Argument(s):

- list_of_files: 1d array-valued argument containing the filenames

#Returns:

- list__mu: 1d array containing chemical potential
- list_SEvec_c: 1d array containing the self-energies taken from qcm output files
"""
function calc_stiff_funct(super_data_M::Array{Float64,1}, super_list_modulevec::Array{Any,1}, match_data_loop_arr::Array{String,1}, M_tol::Float64)
    # Array definitions
    list_kIntegral_stiff = Array{Any,1}()
    list_kIntegral_DOS = Array{Any,1}()
    list_kIntegral_DOS_nocoex = Array{Any,1}()
    # The Matsubara frequency grid is defined here
    w_discretization = super_list_modulevec[1].vals_["w_discretization"]
    beta = super_list_modulevec[1].vals_["beta"]
    zvec = 1.0im*[(2*n+1)*pi/beta for n in 1:w_discretization]
    if super_list_modulevec[1].vals_["not_sparse"] == 1
        AFM_SC = super_list_modulevec[1].vals_["AFM_SC_NOCOEX"]
        cumulant = super_list_modulevec[1].vals_["cumulant"]
        DOS = super_list_modulevec[1].vals_["DOS"]
        notice(logger, "not_sparse option set to 1")
        notice(logger, "AFM_SC_NOCOEX option set to $(AFM_SC)")
        notice(logger, "cumulant option set to $(cumulant)")
        notice(logger, "DOS option set to $(DOS)")
        CHECK = match_data_loop_arr[1]
        notice(logger, "Treating $(CHECK) data")
        if super_list_modulevec[1].vals_["Periodization"] == 1
            notice(logger, "Periodization option set to 1")
            for (l,modulevec_el) in enumerate(super_list_modulevec)
	 #if l % 10 == 0
                modelvec = PeriodizeSC.ModelVector(ts(modulevec_el)[1],ts(modulevec_el)[2],ts(modulevec_el)[3],ts(modulevec_el)[4],zvec[1:w_discretization], modulevec_el.data_[1,1:w_discretization,:,:])
                if (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 0 && cumulant == 0 && DOS == 0
                    push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_SC))
                elseif (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 0 && cumulant == 1 && DOS == 0
	     println("CUM")
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_SC))
                elseif (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 1 && cumulant == 0 && DOS == 0
	     println("Per AFM_SC_1")
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_test))
                elseif (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 1 && cumulant == 1 && DOS == 0
	     println("CUM AFM_SC_1")
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_AFM_SC))
                elseif (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 0 && cumulant == 0 && DOS == 1
                    cond = false # Is there coexistence? : false
                    push!(list_kIntegral_DOS_nocoex,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_nocoex,cond)))
	 elseif (CHECK == "NOCOEX" || CHECK == "AFM") && AFM_SC == 0 && cumulant == 1 && DOS == 1
                    cond = false # Is there coexistence? : false 
                    push!(list_kIntegral_DOS_nocoex,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_nocoex_cum,cond)))
                elseif CHECK == "COEX" && DOS == 1 && cumulant == 0
	     cond = abs(super_data_M[l]) > M_tol
                    if cond
                        push!(list_kIntegral_DOS,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_coex,cond)))
                    else
                        push!(list_kIntegral_DOS_nocoex,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_nocoex,cond)))
                    end
                elseif CHECK == "COEX" && DOS == 1 && cumulant == 1
                    cond = abs(super_data_M[l]) > M_tol
                    if cond
                        push!(list_kIntegral_DOS,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_coex_cum,cond)))
                    else
                        push!(list_kIntegral_DOS_nocoex,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_k_nocoex_cum,cond)))
                    end
                elseif CHECK == "COEX" && DOS == 0 && cumulant == 0
                    cond = abs(super_data_M[l]) > M_tol
                    if cond
                        push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_test))
                    else
                        push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_SC))
                    end
                elseif CHECK == "COEX" && DOS == 0 && cumulant == 1
                    cond = abs(super_data_M[l]) > M_tol
                    if cond
                        push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_AFM_SC))
                    else
                        push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_SC))
                    end
                #end
                end
            end
        elseif super_list_modulevec[1].vals_["Periodization"] == 0
            notice(logger, "Periodization option set to 0")
            for (l,modulevec_el) in enumerate(super_list_modulevec)
	 modelvec = PeriodizeSC.ModelVector(ts(modulevec_el)[1],ts(modulevec_el)[2],ts(modulevec_el)[3],ts(modulevec_el)[4],zvec[1:w_discretization], modulevec_el.data_[1,1:w_discretization,:,:])
                if (CHECK == "NOCOEX" || CHECK == "AFM") && (AFM_SC == 0 || AFM_SC == 1) && DOS == 0
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_trace_G_kintegrand))
                elseif CHECK == "COEX" && DOS == 0
                    push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_cluster_G_kintegrand))
	 elseif (CHECK == "NOCOEX" || CHECK == "COEX" || CHECK == "AFM") && DOS == 1
	     println("yo")
                    cond = false # There is no periodization here, so it doesn't matter -> no coexistence : false
                    push!(list_kIntegral_DOS_nocoex,(l,PeriodizeSC.build_k_precomputed(modelvec,PeriodizeSC.DOS_cluster_trace,cond)))
                end
            end
        end
    elseif super_list_modulevec[1].vals_["not_sparse"] == 0
        print_with_color(:red, "Will do this tomorrow or today")
    end

    return list_kIntegral_stiff, list_kIntegral_DOS, list_kIntegral_DOS_nocoex
end


function save_gf(super_list_modulevec::Array{Any,1})
    notice(logger, "Only saving the Green's functions")
    println("type: ", typeof(super_list_modulevec))
    w_discretization = super_list_modulevec[1].vals_["w_discretization"]
    beta = super_list_modulevec[1].vals_["beta"]
    println("beta: ", beta)
    zvec = 1.0im*[(2*n+1)*pi/beta for n in 1:w_discretization]
    for (l,modulevec_el) in enumerate(super_list_modulevec)
        modelvec = PeriodizeSC.ModelVector(ts(modulevec_el)[1],ts(modulevec_el)[2],ts(modulevec_el)[3],ts(modulevec_el)[4],zvec[1:w_discretization],modulevec_el.data_[1,1:w_discretization,:,:])
        if verbose >= 0
            println("mu: ", modelvec.mu_)
        end
        for n in 1:length(zvec)
            model = PeriodizeSC.Model(modelvec, n)
            if verbose >= 0
                println("w: ", model.w_)
            end
            PeriodizeSC.save_green_function(model,beta)
        end
    end

    return nothing
end


########################################### End of functions and structures ###########################################

# Reading parameters from JSON file
paramsfile = "params_.json"
params = JSON.parsefile(paramsfile)
println(typeof(params))
# Defining own dictionnary for convenience
block_params = Dict{String,Any}(
	"path_to_files_sparse" => params["path_to_files_sparse"],
	"path_to_files" => params["path_to_files"],
	"data_loop_sparse" => params["data_loop_sparse"],
	"data_loop" => params["data_loop"],
	"pattern" => params["pattern"], # Pattern of the files to be loaded
	"beta" => params["beta"],
	"print_mu_dop" => params["Print_mu_dop"], # Parameter used to create stiffness.dat file (file used with SE_G_fig_producer.py)
	"w_discretization" => params["w_discretization"], # Matsubara frequency grid
	"AFM_SC_NOCOEX" => params["AFM_SC_NOCOEX"],
	"DOS" => params["DOS"],
	"cumulant" => params["cumulant"],
	"Periodization" => params["Periodization"],
	"not_sparse" => params["not_sparse"],
	"fout_name" => params["fout_name"],
	"fout_name_DOS" => params["fout_name_DOS"]
	)

# Check if appropriate input was entered
setlevel!(logger,"info")
if block_params["cumulant"] in [0,1]; nothing else; throw(ErrorException("Cumulant takes in 0 or 1 as arguments")) end
if block_params["DOS"] in [0,1]; nothing else; throw(ErrorException("DOS takes in 0 or 1 as arguments")) end
if params["Save_gf_only"] in [0,1]; nothing else; throw(ErrorException("Save_gf_only takes in 0 or 1 as arguments")) end

# Check out if appriopriate inputs entered in data_loop_sparse
match_data_loop = match(r"^(.*?)(?=/)",block_params["data_loop"]).match
match_data_loop_sparse1 = match(r"^(.*?)(?=/)",block_params["data_loop_sparse"][1]).match
match_data_loop_sparse2 = match(r"^(.*?)(?=/)",block_params["data_loop_sparse"][2]).match
match_data_loop_sparse1 == match_data_loop_sparse2 ? nothing : throw(ErrorException("Must enter same file directory (ie, NOCOEX or COEX)"))
match_data_loop_arr = Array{String,1}([match_data_loop,match_data_loop_sparse1,match_data_loop_sparse2])
M_tol = 1e-3

# Decide the name of the file to be created
if block_params["DOS"] in [0,1]; nothing else; throw(ErrorException("DOS takes in either 0 or 1")) end
fout_name = block_params["DOS"] == 1 ? block_params["fout_name_DOS"] : block_params["fout_name"] # DOS can take 1 or 0

# Opening files to read
info(logger, "Reading relevant data from averages*.dat or Loop*.dat files")

super_data_M, super_datap, list_of_files, datap = read_data_loop(block_params,block_params["print_mu_dop"],block_params["pattern"])
list_of_files_and_datap = datap != nothing ? collect(zip(list_of_files,datap)) : collect(zip(list_of_files,super_datap))

# Building list_modulevec defining the instance of modulevec
list_list_modulevec = Array{Array{Any,1},1}()
for (list_file,datap_mu) in list_of_files_and_datap
    length(list_file) != length(datap_mu[:,1]) && throw(File_Input_Error(:list_file))
    if verbose >=0
        println("Length list_file: ", length(list_file), "\tLength datap_mu[:,1]: ", length(datap_mu[:,1]))
    end
    list_modulevec = gen_modulevec_args(block_params,list_file,datap_mu[:,1],params)
    push!(list_list_modulevec, list_modulevec)
    if verbose >= 0
        println("typeof list_modulevec: ",typeof(list_modulevec))
    end
end

if length(list_list_modulevec) > 1 ##<------------------------------ To verify with sparse data
    global super_list_modulevec
    super_list_modulevec = vcat(list_list_modulevec...)
else
    global super_list_modulevec
    super_list_modulevec = list_list_modulevec[1]
end

println(ts(super_list_modulevec[1]))
if verbose >= 0
    println("Size of sEvec_c in StiffnessArray: ", size(super_list_modulevec[1].data_), " | typeof StiffnessArray: ", typeof(super_list_modulevec[1]))
end

########################################### Saving Green's functions ###########################################
setlevel!(logger, "notice")
if params["Save_gf_only"] == 1
    save_gf(super_list_modulevec) # Green's functions are integrated over the Brillouin zone
    notice(logger, "Saved the Green's functions (see function save_gf() for further informations)")  
    exit(0)
elseif params["Save_gf_only"] == 0
    nothing
end

########################################### Beginning of the superfluid stiffness calculations ###########################################
setlevel!(logger, "info")

info(logger, "Entering decision loop regarding if NOCOEX or COEX. Integration over BZ or RBZ")

setlevel!(logger, "notice")

# Computing the superfluid stiffness
list_kIntegral_stiff, list_kIntegral_DOS, list_kIntegral_DOS_nocoex = calc_stiff_funct(super_data_M, super_list_modulevec, match_data_loop_arr, M_tol)

# Summing over the Matsubara frequencies
setlevel!(logger, "info")
beta = block_params["beta"]
info(logger, "Summation in frequency beta = $(@sprintf("%.2f", beta))")
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
    for j in 1:size(list_kIntegral_stiff)[1]-1
        stiffness_tmp = 0.0
        for ii in 1:size(list_kIntegral_stiff[j])[1]-1
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

setlevel!(logger, "warn")
warn(logger, "Program exited")
