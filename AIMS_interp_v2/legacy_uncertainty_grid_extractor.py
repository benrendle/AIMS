import numpy as np
import dill

legacy_file = open("LEGACY_mode_list", "r")

# initialise lists of IDs and mode dictionaries
data_list = []
mode_dict = {}

# append all unique IDs to data_list and append all unique modes to mode_dict
for line in legacy_file:
    column = line.split()
    KIC_id = int(column[0])  # KIC ID of star in catalogue
    n_val = int(column[1])
    l_val = int(column[2])
    freq = float(column[3])
    dfreq = (float(column[4]) + float(column[5])) / 2.  # mean of upper uncertainty and lower uncertainty

    mode_id = "%d:%d" % (n_val, l_val)  # generate a dictionary key in the form 'n:l' for easy lookup

    if KIC_id not in data_list:
        data_list.append(KIC_id)

    if mode_id not in mode_dict:
        mode_dict[mode_id] = []

    # each dictionary entry contains a list of frequency/uncertainty pairs for all occurrences of that mode
    mode_dict[mode_id].append([freq, dfreq])

legacy_file.close()  # finished reading from file here

# initialise a new dict that contains the same keys but each entry is a two-element list of the mean uncertainty of
# that mode and the standard deviation of the mean uncertainty
mode_uncert_dict = {}
for key in mode_dict:
    uncert_list = []

    for mode in mode_dict[key]:
        uncert_list.append(mode[1])

    mean = sum(uncert_list) / float(len(uncert_list))
    x_square = [(x - mean)**2 for x in uncert_list]
    var = sum(x_square) / float((len(uncert_list) if len(uncert_list) > 1 else 2.) - 1.)
    st_dev = np.sqrt(var)

    mode_uncert_dict[key] = [mean, st_dev]

uncert_file = open("legacy_grid", "w")
dill.dump(mode_uncert_dict, uncert_file)
uncert_file.close()
