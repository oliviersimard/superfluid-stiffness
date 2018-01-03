using JSON
using NPZ
using Glob
using StiffSquare2.PeriodizeSC

paramsfile = "params.json"
params = JSON.parsefile(paramsfile)
path_to_files = params["path_to_files"]
pattern = params["pattern"]
beta = params["beta"]
fout_name = params["fout_name"]

data_file_s = open(readdlm, params["data_loop"][1])
data_file = open(readdlm, params["data_loop"][2])
data_file_header = data_file[1,:]
data_file_s_header = data_file_s[1,:]
#ind=1
#index=0
#while getindex(data_file_header,ind)!="mu"
#    ind+=1
#    index=ind
#end

indmu_s = find(x->x=="mu",data_file_s_header)
inddop_s = find(x->x=="ave_mu",data_file_s_header)

indmu = find(x->x=="mu",data_file_header)
inddop = find(x->x=="ave_mu",data_file_header)


data_file_s_mu = data_file_s[:,indmu_s]
data_file_s_dop = data_file_s[:,inddop_s]

data_file_mu = data_file[:,indmu]
data_file_dop = data_file[:,inddop]

#data_mu = data_file[:,index]
data_file_s_mu = filter(x->x!="mu",data_file_s_mu)
data_file_s_dop = filter(x->x!="ave_mu",data_file_s_dop)
data_file_s_mu = convert(Array{Float64,1}, data_file_s_mu)

data_file_mu = filter(x->x!="mu",data_file_mu)
data_file_dop = filter(x->x!="ave_mu",data_file_dop)
data_file_mu = convert(Array{Float64,1}, data_file_mu)

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

datap_h = hcat(data_file_mu,data_file_dop)
#datap_h_sorted = sort(datap_h,1; alg=MergeSort)
#println(datap_new_h_sorted)
super_datap = vcat(datap_h,datap_s_h)

list_of_files = glob(string(path_to_files[2],pattern))
list_of_files_s = glob(string(path_to_files[1],pattern))

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

map!(x->parse(Int64,x),list_num)
map!(x->parse(Int64,x),list_num_s)
zip_num_files = collect(zip(list_num,list_of_files))
zip_num_files_s = collect(zip(list_num_s,list_of_files_s))
sort!(zip_num_files)
sort!(zip_num_files_s)

zvec = 1.0im*[(2*n+1)*pi/beta for n in 1:400]

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

super_list_SEvec_c = vcat(list_SEvec_c, list_SEvec_c_s)
super_list_t_mu = vcat(list_t_mu,list_t_mu_s)

list_kIntegral_stiff = []
for l in 1:length(super_list_t_mu)
    modelvec = PeriodizeSC.ModelVector(super_list_t_mu[l][1], super_list_t_mu[l][2], super_list_t_mu[l][3], super_list_t_mu[l][4], zvec[1:400], super_list_SEvec_c[l][1, 1:400, :, :])
    push!(list_kIntegral_stiff,PeriodizeSC.calcintegral_RBZ(modelvec, PeriodizeSC.make_stiffness_kintegrand_test))#Change kintegrand if COEX or NOCOEX
    println(size(list_kIntegral_stiff[l]))
end

println("Summation in frequency")
stiffness = []
for j in 1:size(list_kIntegral_stiff)[1]-1
    stiffness_tmp = 0.0
    for ii in 1:size(list_kIntegral_stiff[j])[1]
        stiffness_tmp+=list_kIntegral_stiff[j][ii, 2]
        #push!(stiffness_tmp, 0.5 * (list_kIntegral_stiff[j][ii, 2] + list_kIntegral_stiff[j][ii+1, 2]) * (list_kIntegral_stiff[j][ii+1, 1] - list_kIntegral_stiff[j][ii, 1]) / (2.0*pi))
    end
    push!(stiffness,(1/beta)*real(stiffness_tmp))
end

stiffness*=2 #For the spin projections
mu_dop_stiff = hcat(super_datap[1:length(stiffness),1],super_datap[1:length(stiffness),2],stiffness)
println("Stiffness = ", mu_dop_stiff)
writedlm(fout_name, mu_dop_stiff, "\t\t")

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
