import dill
import random
import numpy as np
import copy
import sys

import model
import AIMS_configure as config

print "Please select the correct grid type."
print "0: Main sequence"
print "1: Red giant"
raw_option = raw_input("Choice: ")
option = int(raw_option)

if option not in [0, 1]:
    print "The input you have selected does not correspond to one of the options"
    sys.exit(1)

raw_sample = raw_input("Select the number of test models to create: ")
samples = int(raw_sample)

random.seed()
np.random.seed()

# read in model grid for the grid currently in use from the AIMS_configure file
grid_file_name = config.binary_grid
grid_file = open(grid_file_name, "r")
data = dill.load(grid_file)
grid_file.close()

if config.user_params != data.user_params:
    print "User params in configuration file do not match user params in grid file"
    print "Grid: %s" % str(data.user_params)
    print "AIMS_configure: %s" % str(config.user_params)
    sys.exit(1)

# read in uncertainty grid for modes generated from LEGACY data
uncert_file = open("legacy_grid", "r")
uncert_dict = dill.load(uncert_file)
uncert_file.close()

# randomly select the track and model indices
track_index = random.randint(0, len(data.tracks))
model_index = random.randint(0, len(data.tracks[track_index].models))

# pick out the chosen track and model, extract its modes and reference frequency
trakk = data.tracks[track_index]
modelle = data.tracks[track_index].models[model_index]
mode_tuple = modelle.modes
mode_list = []
ref_freq = modelle.glb[model.ifreq_ref]

# generate a new list of modes with the frequency scaled back to uHz and the n, l columns swapped to be consistent
# with AIMS observable list file
for index in xrange(len(mode_tuple)):
    mode_list.append(list(mode_tuple[index]))
    mode_list[index][2] = mode_list[index][2] * ref_freq
    mode_list[index][0], mode_list[index][1] = mode_list[index][1], mode_list[index][0]

# dump true model parameters into a file for comparison with AIMS results
output_model_file = open("Model_params", "w")
for index in xrange(len(trakk.grid_params)):
    line = "%s = %f\n" % (trakk.grid_params[index], trakk.params[index])
    output_model_file.write(line)

line_age = "Age = %f\n" % modelle.glb[model.iage]
output_model_file.write(line_age)
output_model_file.close()

# extract constraint parameters
T_mean = modelle.glb[model.itemperature]
numax_mean = modelle.numax
FeH_mean = modelle.FeH

# generate 'n' model observable files chosen by user
for sample in xrange(samples):
    fname = "Model_%d" % sample
    output_file = open(fname, "w")

    if option == 0:  # generate file for main sequence star
        # print all allowed modes into the list file
        for mode in mode_list:
            if 0.5*numax_mean < mode[2] < 1.5*numax_mean:  # chosen mode cutoff to not overload AIMS
                mode_key = "%d:%d" % (mode[1], mode[0])

                if mode_key in uncert_dict:  # check if mode exists in LEGACY grid

                    # create a new copy of list in memory to allow for perturbation of modes without accidentally changing
                    # these values in the master mode list
                    new_mode = copy.deepcopy(mode)
                    mean = mode[2]
                    st_dev = uncert_dict[mode_key][0]

                    # copy over exact values of the frequencies if this is Model_0 else assign a random value sampled from
                    # a Gaussian distribution centred around the mean and with width given by LEGACY data
                    new_mode[2] = mean if sample == 0 else np.random.normal(mean, st_dev)
                    new_mode[3] = st_dev
                    mode_list_string = []

                    # print modes in a way that AIMS understands
                    for value in new_mode:
                        mode_list_string.append(str(value))
                    mode_string = " ".join(mode_list_string)
                    output_file.write(mode_string + "\n")
            else:
                print "%.2f ignored as it is not between %.2f and %.2f" % (mode[2], 0.5*numax_mean, 1.5*numax_mean)

        # uncertainties for constraints chosen as the modal average of the uncertainties in LEGACY data due to a 53 to 13
        # ratio of this value to all others
        T_width = 77.  # 53 to 13
        T_value = T_mean if sample == 0 else np.random.normal(T_mean, T_width)
        T_string = "Teff %f %f\n" % (T_value, T_width)
        output_file.write(T_string)

        # uncertainties for constraints chosen as the mean average of the uncertainties in LEGACY data due to no clear
        # relationship between uncertainties and the values themselves
        numax_width = 11.3605263158
        numax = numax_mean if sample == 0 else np.random.normal(numax_mean, numax_width)
        numax_string = "numax %f %f\n" % (numax, numax_width)
        output_file.write(numax_string)

        # uncertainties for constraints chosen as the modal average of the uncertainties in LEGACY data due to a 55 to 11
        # ratio of this value to all others
        FeH_width = 0.1  # 55 to 11
        FeH = FeH_mean if sample == 0 else np.random.normal(FeH_mean, FeH_width)
        FeH_string = "Fe_H %f %f\n" % (FeH, FeH_width)
        output_file.write(FeH_string)

        output_file.close()

    elif option == 1:  # generate file for red giant star
        # print all allowed modes into the list file
        for mode in mode_list:
            if 0.5 * numax_mean < mode[2] < 1.5 * numax_mean:  # chosen mode cutoff to not overload AIMS
                mode_key = "%d:%d" % (mode[1], mode[0])

                # create a new copy of list in memory to allow for perturbation of modes without accidentally changing
                # these values in the master mode list
                new_mode = copy.deepcopy(mode)
                mean = mode[2]
                st_dev = 0.0168  # Kepler56 mean radial mode uncertainty

                # copy over exact values of the frequencies if this is Model_0 else assign a random value sampled from
                # a Gaussian distribution centred around the mean and with width given by representative uncertainty
                # from Kepler56 data
                new_mode[2] = mean if sample == 0 else np.random.normal(mean, st_dev)
                new_mode[3] = st_dev
                mode_list_string = []

                # print modes in a way that AIMS understands
                for value in new_mode:
                    mode_list_string.append(str(value))
                mode_string = " ".join(mode_list_string)
                output_file.write(mode_string + "\n")
            else:
                print "%.2f ignored as it is not between %.2f and %.2f" % (mode[2], 0.5 * numax_mean, 1.5 * numax_mean)

        # uncertainties for constraints chosen from Kepler56 data
        T_width = 97.
        T_value = T_mean if sample == 0 else np.random.normal(T_mean, T_width)
        T_string = "Teff %f %f\n" % (T_value, T_width)
        output_file.write(T_string)

        numax_width = 1.4
        numax = numax_mean if sample == 0 else np.random.normal(numax_mean, numax_width)
        numax_string = "numax %f %f\n" % (numax, numax_width)
        output_file.write(numax_string)

        FeH_width = 0.16
        FeH = FeH_mean if sample == 0 else np.random.normal(FeH_mean, FeH_width)
        FeH_string = "Fe_H %f %f\n" % (FeH, FeH_width)
        output_file.write(FeH_string)

        output_file.close()

