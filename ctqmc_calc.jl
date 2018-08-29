using StiffSquare2.PeriodizeSC
using NPZ
using Glob

# You have to run plot_phase_diagram.py before running this program properly

_Nc = 4
t = 1.0; tp = -0.3; tpp = 0.2

macro assertion(ex, text)
    :($ex ? nothing : error("Assertion failed: ", $(text)))
end

function grasp_mu(mu_file_path::String)
    data_mu = readdlm(mu_file_path)[:,2][end]
    return data_mu
end

function load_iwn_list(w_file_path)
    data_iwn = readdlm(w_file_path)[:,1]
    return convert(Array{Complex{Float64},1},data_iwn)
end

pattern = r"plaquetteC2vSc"
list_subdirs = readdir()
list_of_npy_SE_files = Array{String,1}()
list_of_npy_GE_files = Array{String,1}()
list_list_iwn = Array{Array{Complex{Float64},1},1}()
list_mu = Array{Float64,1}()
list_params = Array{Tuple{Float64,Float64,Float64},1}()
for subdirs in list_subdirs
    if ismatch(pattern,subdirs)
        matches_params = [el for el in eachmatch(r"[-+]?\d*\.\d+|\d+",subdirs)]
        U = parse(Float64, matches_params[3].match)
        n = parse(Float64, matches_params[2].match)
        beta = parse(Float64,matches_params[4].match)
        npy_files = glob(subdirs*"/*.npy")
        mu_files = glob(subdirs*"/*mu.dat")
        iwn_files = glob(subdirs*"/*w.dat")
        mu = grasp_mu(mu_files...)
        iwn_list = load_iwn_list(iwn_files...)
        push!(list_of_npy_SE_files,npy_files[1])
        push!(list_of_npy_GE_files,npy_files[2])
        push!(list_mu,mu)
        push!(list_list_iwn,iwn_list)
        push!(list_params,(n,U,beta))
    end
end
mu_SE_arr = collect(zip(list_mu,list_of_npy_SE_files))
println(list_params)
function periodize_SC(arg::Array{Complex{Float64},2})
    nambu_periodized = zeros(Complex{Float64}, _Nc)
    g_up = arg[1:4,1:4]; f = arg[1:4,5:end]
    f_dag = arg[5:end,1:4]; g_down = arg[5:end,5:end]
    elperiodized = [g_up, f, f_dag, g_down]
    r_sites = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
    K_sites = pi*deepcopy(r_sites)

    for (l,el) in enumerate(elperiodized)
        for i in 1:_Nc
            for j in 1:_Nc
                for K in K_sites
                    nambu_periodized[l] += 1.0/_Nc * exp(1.0im*dot(K,r_sites[i]-r_sites[j])) * el[i,j] 
	 end
            end
        end
    end
    nambu_periodized = reshape(nambu_periodized,(2,2))
    return 0.25*nambu_periodized
end

list_Kintegral_stiff = Array{Tuple{Tuple{Float64,Float64,Float64},Array{Complex{Float64},2}},1}()
list_dop = Array{Float64,1}()
for l in 1:1:length(mu_SE_arr)
    SEvec = npzread(list_of_npy_SE_files[l])
    params_tuple = list_params[l]
    modelvec = PeriodizeSC.ModelVector(t,tp,tpp,list_mu[l],list_list_iwn[l],SEvec)
    cal = PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_SC)
    push!(list_Kintegral_stiff,(params_tuple,cal))
    push!(list_dop,params_tuple[1])
end 
println("list_dop: ", list_dop)
stiffness = Array{Float64,1}()
for j in 1:size(list_Kintegral_stiff)[1]
    stiffness_tmp = 0.0
    println(list_Kintegral_stiff[j][1])
    for ii in 1:size(list_Kintegral_stiff[j][2])[1]
        stiffness_tmp+=list_Kintegral_stiff[j][2][ii, 2]
    end
    beta = list_Kintegral_stiff[j][1][3]
    println("beta: ", beta)
    push!(stiffness,(1/beta)*real(stiffness_tmp))
    println((1/beta)*imag(stiffness_tmp))
end

@assertion(length(list_dop)==length(stiffness), "Length of list_dop differs from that of stiffness.")

dop_stiff = hcat(list_dop,stiffness)
println("Stiffness: ", stiffness)

writedlm("first_try_.dat", dop_stiff, "\t\t")

