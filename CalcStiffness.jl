using JSON
using NPZ
using Glob
using Memento
using StiffSquare2.PeriodizeSC

if isfile("./Julia.log")
    rm("./Julia.log")
end

print_with_color(:green, "Extracting info from params.json\n")

paramsfile = "params.json"
params = JSON.parsefile(paramsfile)
path_to_files = params["path_to_files"]
path_to_files_U12 = params["path_to_files_U12"]
pattern = params["pattern"]
beta = params["beta"]
print_mu_dop = params["Print_mu_dop"]
AFM_SC_NOCOEX = params["AFM_SC_NOCOEX"] ## Set to 1 if using AFM_SC function to evaluate NOCOEX data. Set to 0
w_discretization = params["w_discretization"]  # if not.
DOS = params["DOS"]
cumulant = params["cumulant"]
Periodization = params["Periodization"]
Other_U = params["Other_U"]

if DOS == 1
    fout_name = params["fout_name_DOS"]
elseif DOS == 0
    fout_name = params["fout_name"]
else 
    throw(ErrorException("DOS takes in either 0 or 1"))
end

print_with_color(:green, "Setting parameters for log file\n")

logger = getlogger(current_module())
push!(logger, DefaultHandler("Julia.log"))
Memento.config("debug";fmt="[{level}|{name}]:{msg}")

data_file_s = open(readdlm, params["data_loop"][1])
data_file = open(readdlm, params["data_loop"][2])
data_file_U12 = open(readdlm, params["data_loop_U12"])
data_file_header = data_file[1,:]
data_file_s_header = data_file_s[1,:]
data_file_U12_header = data_file_U12[1,:]
#ind=1
#index=0
#while getindex(data_file_header,ind)!="mu"
#    ind+=1
#    index=ind
#end
setlevel!(logger, "info")
info(logger, "Reading relevant data from averages*.dat or Loop*.dat files")
indmu_s = find(x->x=="mu",data_file_s_header)
inddop_s = find(x->x=="ave_mu",data_file_s_header)
indM_s = find(x->x=="ave_M",data_file_s_header)

indmu = find(x->x=="mu",data_file_header)
inddop = find(x->x=="ave_mu",data_file_header)
indM = find(x->x=="ave_M",data_file_header)

indU12 = find(x->x=="mu",data_file_U12_header)
inddopU12 = find(x->x=="ave_mu",data_file_U12_header)
indMU12 = find(x->x=="ave_M",data_file_U12_header)

data_file_s_mu = data_file_s[:,indmu_s]
data_file_s_dop = data_file_s[:,inddop_s]
data_file_s_M = data_file_s[:,indM_s]

data_file_mu = data_file[:,indmu]
data_file_dop = data_file[:,inddop]
data_file_M = data_file[:,indM]

data_file_mu_U12 = data_file_U12[:,indU12]
data_file_dop_U12 = data_file_U12[:,inddopU12]
data_file_M_U12 = data_file_U12[:,indMU12]

#data_mu = data_file[:,index]
data_file_s_mu = filter(x->x!="mu",data_file_s_mu)
data_file_s_dop = filter(x->x!="ave_mu",data_file_s_dop)
data_file_s_M = filter(x->x!="ave_M",data_file_s_M)
data_file_s_mu = convert(Array{Float64,1}, data_file_s_mu)

data_file_mu = filter(x->x!="mu",data_file_mu)
data_file_dop = filter(x->x!="ave_mu",data_file_dop)
data_file_M = filter(x->x!="ave_M",data_file_M)
data_file_mu = convert(Array{Float64,1}, data_file_mu)

data_file_mu_U12 = filter(x->x!="mu",data_file_mu_U12)
data_file_dop_U12 = filter(x->x!="ave_mu",data_file_dop_U12)
data_file_M_U12 = filter(x->x!="ave_M",data_file_M_U12)

datap_s = collect(zip(data_file_s_mu,data_file_s_dop))
datap = collect(zip(data_file_mu,data_file_dop))
#datap_new = [elem for elem in datap if isa(elem[2],Float64) && isa(elem[1],Float64)]

#list_dop=[]
#for elem in datap_new
#    push!(list_dop,elem[2])
#end

#list_mu=[]
#for elem in datap_new
#    push!(list_mu,elem[1])
#end

datap_s_h = hcat(data_file_s_mu,data_file_s_dop)
#datap_s_h_sorted = sort(datap_s_h,1; alg=MergeSort)
datap_other_U_h = hcat(data_file_mu_U12,data_file_dop_U12)
datap_other_U_h = convert(Array{Float64,2},datap_other_U_h)
println("datap_other_U_h = ", datap_other_U_h)
datap_h = hcat(data_file_mu,data_file_dop)
#datap_h_sorted = sort(datap_h,1; alg=MergeSort)
#println(datap_new_h_sorted)
super_datap = vcat(datap_h,datap_s_h)
if print_mu_dop == 1
    if Other_U == 0
        writedlm("stiffness.dat", super_datap, "\t\t")
        print_with_color(:red, "Printed stiffness.dat for later use (SEvec_fig.py)\n")
        exit(0)
    elseif Other_U == 1
        writedlm("stiffness.dat", datap_other_U_h, "\t\t")
        print_with_color(:red, "Printed stiffness.dat for later use (SEvec_fig.py)\n")
        exit(0)
    end
end
super_data_M = vcat(data_file_M,data_file_s_M)

setlevel!(logger, "info")
info(logger, "Reading SEvec*.dat.npy files obtained from qcm")
list_of_files = glob(string(path_to_files[2],pattern))
list_of_files_s = glob(string(path_to_files[1],pattern))
list_of_files_other_U = glob(string(path_to_files_U12,pattern))
#println("list_of_files_other_U = ",list_of_files_other_U)

list_num = []
for files in list_of_files
    m=match(r"\d+",files,length(files)-10)
    push!(list_num,m.match)
end

list_num_s = []
for files in list_of_files_s
    m=match(r"\d+",files,length(files)-10)
    push!(list_num_s,m.match)
end

m3 = match(r"^(.*?)(?=/)",params["data_loop_U12"])

if Other_U == 1 && m3.match == "NOCOEX"
    length_file = 6
else
    length_file = 10
end

list_num_other_U = []
for files in list_of_files_other_U
    m=match(r"\d+",files,length(files)-length_file)
    push!(list_num_other_U,m.match)
end

map!(x->parse(Int64,x),list_num,list_num)
map!(x->parse(Int64,x),list_num_s,list_num_s)
map!(x->parse(Int64,x),list_num_other_U,list_num_other_U)
zip_num_files = collect(zip(list_num,list_of_files))
zip_num_files_s = collect(zip(list_num_s,list_of_files_s))
println("list_num_other_U = ", list_num_other_U)
zip_num_files_other_U = collect(zip(list_num_other_U,list_of_files_other_U))
sort!(zip_num_files)
sort!(zip_num_files_s)
sort!(zip_num_files_other_U)

zvec = 1.0im*[(2*n+1)*pi/beta for n in 1:w_discretization]

list_SEvec_c = []
list_t_mu = []
for (num,files) in zip_num_files
    sEvec_c = npzread(files)
    push!(list_t_mu, (params["t"], params["tp"], params["tpp"], datap_h[num+1]))#datap[num+1][1]))
    push!(list_SEvec_c,sEvec_c)
end

list_SEvec_c_s = []
list_t_mu_s = []
for (num,files) in zip_num_files_s
    sEvec_c = npzread(files)
    push!(list_t_mu_s, (params["t"], params["tp"], params["tpp"], datap_s_h[num+1]))
    push!(list_SEvec_c_s,sEvec_c)
end

list_SEvec_c_other_U = []
list_t_mu_other_U = []
for (num,files) in zip_num_files_other_U
    println("num = ", num)
    println("files = ", files)
    sEvec_c = npzread(files)
    push!(list_t_mu_other_U,(params["t"],params["tp"],params["tpp"],datap_other_U_h[num+1]))
    push!(list_SEvec_c_other_U,sEvec_c)
end
println("list_t_mu_other_U = ", list_t_mu_other_U)
#println("list_SEvec_c_other_U = ", list_SEvec_c_other_U)

super_list_SEvec_c = vcat(list_SEvec_c, list_SEvec_c_s)
super_list_t_mu = vcat(list_t_mu,list_t_mu_s)

####################################################################################################################
##Verifying input validity##########################################################################################
####################################################################################################################
setlevel!(logger,"notice")
if cumulant == 1
    notice(logger, "Cumulant set to 1: only methods to compute superfluid stiffness exist")
elseif cumulant == 0
    notice(logger, "Cumulant set to 0")
else
    throw(ErrorException("Cumulant takes in 0 or 1 as arguments"))
end

if DOS == 1
    notice(logger, "DOS set to 1")
elseif DOS == 0
    notice(logger, "DOS set to 0")
else
    throw(ErrorException("DOS takes in 0 or 1 as arguments"))
end
####################################################################################################################
setlevel!(logger, "info")
info(logger, "Entering decision loop regarding if NOCOEX or COEX. Integration over BZ or RBZ")
m1 = match(r"^(.*?)(?=/)",params["data_loop"][1])
m2 = match(r"^(.*?)(?=/)",params["data_loop"][2])
#m3 = match(r"^(.*?)(?=/)",params["data_loop_U12"])
println("Printing m3 : ", m3)
println("Printing m1 and m2 : ", m1, m2)
M_mean_field_tol = 1e-5
list_kIntegral_stiff = []
list_kIntegral_DOS = []
if Other_U == 0
    #print_with_color(:red, "Other_U set to 0\n")
    setlevel!(logger, "notice")
    notice(logger, "Other_U option set to 0")
    if Periodization == 1
        #print_with_color(:red, "Periodization set to 1\n")
        setlevel!(logger, "notice")
        notice(logger, "Periodization option set to 1")
        if m1.match == "NOCOEX" && m2.match == "NOCOEX" && AFM_SC_NOCOEX == 1
            setlevel!(logger, "notice")
            notice(logger, "Entered loop with NOCOEX and AFM_SC_NOCOEX = 1")
            for l in 1:length(super_list_t_mu)
                modelvec = PeriodizeSC.ModelVector(super_list_t_mu[l][1], super_list_t_mu[l][2], super_list_t_mu[l][3], super_list_t_mu[l][4], zvec[1:w_discretization], super_list_SEvec_c[l][1, 1:w_discretization, :, :])
                if cumulant == 1
                    push!(list_kIntegral_stiff,0.5*PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_AFM_SC))
                elseif cumulant == 0
                    if DOS == 1
                        #push!(list_kIntegral_DOS,0.5*PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_kintegrand_DOS_k_coex))
                        push!(list_kIntegral_DOS,0.5*PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.DOS_k_coex,gridK=75))
                        println(size(list_kIntegral_DOS[l]))
                    elseif DOS == 0
                        #push!(list_kIntegral_stiff,0.5*PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_test)) ##Differences if using calcintegral_RBZ (worst) or calcintegral_BZ
                        push!(list_kIntegral_stiff,0.5*PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.stiffness_test,gridK=75))
                        println(size(list_kIntegral_stiff[l]))
                    end
                end
            end
        else
            #print_with_color(:red, "AFM_SC_NOCOEX option set to 0\n")
            setlevel!(logger, "notice")
            notice(logger, "AFM_SC_NOCOEX option set to 0")
            for l in 1:length(super_list_t_mu)
                modelvec = PeriodizeSC.ModelVector(super_list_t_mu[l][1], super_list_t_mu[l][2], super_list_t_mu[l][3], super_list_t_mu[l][4], zvec[1:w_discretization], super_list_SEvec_c[l][1, 1:w_discretization, :, :])
                if ismatch(r"N[\w]*X",params["data_loop"][1]) && ismatch(r"N[\w]*X",params["data_loop"][2])
                    if cumulant == 1
	         push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_SC))
                    elseif cumulant == 0
                        if DOS == 1
                            push!(list_kIntegral_DOS,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_kintegrand_DOS_k_nocoex))
	             #push!(list_kIntegral_DOS,2.0*PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.DOS_k_nocoex,gridK=75))
                            println(size(list_kIntegral_DOS[l]))
                        elseif DOS == 0
                            push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_SC))#Change kintegrand if COEX or NOCOEX
                            #push!(list_kIntegral_stiff,PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_SC,gridK=75))
	             println(size(list_kIntegral_stiff[l]))
	         end
                    end
                elseif ismatch(r"^C[\w]*X",params["data_loop"][1]) && ismatch(r"^C[\w]*X",params["data_loop"][2]) || error("Both averages*.dat must come whether from NOCOEX or COEX folders")
                    if cumulant == 1
	         if abs(super_data_M[l])>M_mean_field_tol
	             push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_AFM_SC))
	         else
	             push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_cum_SC))
	         end
                    elseif cumulant == 0
                        if DOS == 1
	             if abs(super_data_M[l])>M_mean_field_tol
	                 #push!(list_kIntegral_DOS,PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_kintegrand_DOS_k_coex))
	                 push!(list_kIntegral_DOS,PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.DOS_k_coex,gridK=60))
	             else
	                 #push!(list_kIntegral_DOS,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_kintegrand_DOS_k_nocoex))
	                 push!(list_kIntegral_DOS,PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.DOS_k_nocoex,gridK=60)) ##coex or nocoex
	             end
	             println(size(list_kIntegral_DOS[l]))
                        elseif DOS == 0
                            if abs(super_data_M[l])>M_mean_field_tol
                                push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_test))
	                 #push!(list_kIntegral_stiff,PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.stiffness_test,gridK=50))
                            else 
                                push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_test))
	                 #push!(list_kIntegral_stiff,2.0*PeriodizeSC.sum_RBZ(modelvec,PeriodizeSC.stiffness_SC,gridK=50))
	             end
	             println(size(list_kIntegral_stiff[l]))
	         end
                    end
                end        
            end
        end
    elseif Periodization == 0
        #print_with_color(:red, "Periodization set to 0\n")
        setlevel!(logger, "notice")
        notice(logger, "Periodization option set to 0")
        for l in 1:length(super_list_t_mu)
            modelvec = PeriodizeSC.ModelVector(super_list_t_mu[l][1], super_list_t_mu[l][2], super_list_t_mu[l][3], super_list_t_mu[l][4], zvec[1:w_discretization], super_list_SEvec_c[l][1, 1:w_discretization, :, :])
            if ismatch(r"N[\w]*X",params["data_loop"][1]) && ismatch(r"N[\w]*X",params["data_loop"][2])
                println("Function to come")
            elseif ismatch(r"^C[\w]*X",params["data_loop"][1]) && ismatch(r"^C[\w]*X",params["data_loop"][2])
                if abs(super_data_M[l])>M_mean_field_tol
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_cluster_G_kintegrand))
                else
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec, PeriodizeSC.make_stiffness_cluster_G_kintegrand))
                end
            end
        end
    else
        print_with_color(:red,"Periodization takes in as argument either 0 or 1")
    end
elseif Other_U == 1
    #print_with_color(:red, "Other_U set to 1\n")
    setlevel!(logger, "notice")
    notice(logger, "Other_U option set to 1")
    if Periodization == 1
        print_with_color(:red, "Periodization set to 1\n")
        setlevel!(logger, "notice")
        notice(logger, "Periodization option set to 1")
        for l in 1:length(list_t_mu_other_U)
            modelvec = PeriodizeSC.ModelVector(list_t_mu_other_U[l][1],list_t_mu_other_U[l][2],list_t_mu_other_U[l][3],list_t_mu_other_U[l][4],zvec[1:w_discretization],list_SEvec_c_other_U[l][1,1:w_discretization,:,:])
            if m3.match == "NOCOEX"
                push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_SC))
            elseif m3.match == "COEX"
                if abs(data_file_M_U12[l])>M_mean_field_tol
                    push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_test))
                else
	     push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_BZ(modelvec,PeriodizeSC.make_stiffness_kintegrand_SC))
                end
            end
        end
    elseif Periodization == 0
        println("This part is under development")
    end
end


setlevel!(logger, "info")
info(logger, "Summation in frequency beta = $(@sprintf("%.2f", beta))")
stiffness = []

if DOS == 1
    for j in 1:size(list_kIntegral_DOS)[1]-1
        stiffness_tmp = 0.0
        for ii in 1:size(list_kIntegral_DOS[j])[1]
            stiffness_tmp+=list_kIntegral_DOS[j][ii,2]
        end
        push!(stiffness,(1/beta)*real(stiffness_tmp))
        println((1/beta)*imag(stiffness_tmp))
    end
elseif DOS == 0     
    for j in 1:size(list_kIntegral_stiff)[1]-1
        stiffness_tmp = 0.0
        for ii in 1:size(list_kIntegral_stiff[j])[1]
            stiffness_tmp+=list_kIntegral_stiff[j][ii, 2]
            #push!(stiffness_tmp, 0.5 * (list_kIntegral_stiff[j][ii, 2] + list_kIntegral_stiff[j][ii+1, 2]) * (list_kIntegral_stiff[j][ii+1, 1] - list_kIntegral_stiff[j][ii, 1]) / (2.0*pi))
        end
        push!(stiffness,(1/beta)*real(stiffness_tmp))
        println((1/beta)*imag(stiffness_tmp))
    end
end
#stiffness*=2 #For the spin projection ##Eventually include this directly in calcintegral_BZ and calcintegral_RBZ
if Other_U == 0
    mu_dop_stiff = hcat(super_datap[1:length(stiffness),1],super_datap[1:length(stiffness),2],stiffness)
elseif Other_U == 1
    mu_dop_stiff = hcat(datap_other_U_h[1:length(stiffness),1],datap_other_U_h[1:length(stiffness),2],stiffness)
end
println("Stiffness = ", mu_dop_stiff)
writedlm(fout_name, mu_dop_stiff, "\t\t")

setlevel!(logger, "warn")
warn(logger, "Program exited")

#ll=[]
#for dop in collect(zip(list_mu,list_dop))
#    for i in 1:length(data_mu)
#        if dop[1]==data_mu[i]
#            push!(ll,dop)
#        end
#    end
#end
#
#u = unique(ll)

##still not enough: should get 101 elements instead of 93 for data_mu
#
#for el in u
#    if count(x->x[1]==el[1],u)!=1
#        print(el[2])
#    end
#end
#
#
