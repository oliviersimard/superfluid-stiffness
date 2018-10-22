import plot_utilities as pu # It is located in pylib directory (see .bash_profile)

## Main
################################################################

# Making dictionnary containing all the data of all the files. Sorted according to the relevant parameters.

if __name__ == '__main__':

    dict_vals_files = {}
    varied_parameter = ""

    list_tuple_vars = pu.get_tuple_vars()
    list_tuple_vars_sorted = sorted(list_tuple_vars, key=lambda x: (x[0],x[1],x[2],x[3]))
    tup_vals_gathered = []

    for tup_vals in list_tuple_vars_sorted:
        if tup_vals_gathered == []:
            tup_vals_gathered.append(tup_vals[0:4])
        elif tup_vals[0:4] not in tup_vals_gathered:
            tup_vals_gathered.append(tup_vals[0:4])

    for tup_vals in tup_vals_gathered: # While loop to collect all tuples sharing same parameters' values
        list_tmp = []
        #print("tup_vals: ", tup_vals)
        for tup_vals_file in list_tuple_vars_sorted:
            if tup_vals == tup_vals_file[0:4]:
                list_tmp.append(tup_vals_file[-1])
        dict_vals_files[tup_vals] = list_tmp

    ## Useful in order to convert quantities in Nambu-like formalism. Use with function convert_converged_green().

    dict_ones, dict_lasts = pu.get_last_ones(dict_vals_files, "greenR", 10) # Should be Green's function if order parameter to be plotted

    #print("dict_ones: ", dict_ones)

    pu.print_SC_order_parameter(dict_ones, "greenR")

    ## If you want to convert the Green's functions and self-energies into proper Nambu formalism, uncomment the following line. It is necessary to compute the superfluid stiffness. 
    ## Otherwise, comment.

    #pu.convert_converged_green(dict_ones, "selfR") #<------------------------------------------------------- To compute superfluid stiffness, uncomment!

    #pu.get_mean_plots(dict_vals_files, "greenR", 10, "g12_nambu_Re", "g12_nambu_Im") ## Can be a list. Input string is to specify which quantities to plot.

    #pu.get_mean_SC_order_parameter_plots(dict_ones, "greenR")

exit(0)
