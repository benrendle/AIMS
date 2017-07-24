import dill
import numpy as np
import sys

leg_tab = open("LEGACY_table", "r")
leg_numax = open("LEGACY_numax", "r")
leg_file = open("LEGACY_mode_list", "r")

star_dict = {}

for line in leg_tab:
    column = line.strip().split("\t")
    KIC = int(column[0])
    if KIC not in star_dict:  # check if ID is in catalogue and create a new dict and dict entry with KIC as key
        star_dict[KIC] = {}

    # isolate temperature and uncertainty and assign to correct dict entry
    temp_col = column[9]
    tStart = temp_col.find('-') + 2
    tEnd = temp_col.find('(')
    tEnd = tEnd if tEnd != -1 else None
    temp = float(temp_col[:4])
    dtemp = float(temp_col[tStart:tEnd])
    star_dict[KIC]['temp'], star_dict[KIC]['dtemp'] = temp, dtemp

    # isolate metallicity and uncertainty and assign to correct dict entry
    met_col = column[10]
    mStart = met_col.find('r') + 3
    mEnd = met_col.find('(')
    mEnd = mEnd if mEnd != -1 else None
    metal = float(met_col[:5] if met_col[0] == '-' else met_col[:4])
    dmetal = float(met_col[mStart:mEnd])
    star_dict[KIC]['metal'], star_dict[KIC]['dmetal'] = metal, dmetal

for line in leg_numax:
    column = line.strip().split("\t")
    KIC = int(column[0])
    if KIC not in star_dict:
        print "KIC not found"
        break

    nu_col = column[1]
    nuEnd = nu_col.find('$')
    loStart = nu_col.find('-') + 1
    loEnd = nu_col.find('^') - 1
    loNumax = float(nu_col[loStart:loEnd])
    hiStart = nu_col.find('+') + 1
    hiNumax = float(nu_col[hiStart:-2])
    numax = float(nu_col[:nuEnd])
    dnumax = (loNumax + hiNumax) / 2.

    frac_sol = 30. / 3090.
    sol_contrib = frac_sol * numax
    dnumax_better = np.sqrt(dnumax**2 + sol_contrib**2)

    star_dict[KIC]['numax'], star_dict[KIC]['dnumax'] = numax, dnumax_better

for line in leg_file:
    column = line.split()
    KIC = int(column[0])
    if KIC not in star_dict:
        print "KIC not found"
        break
    elif 'modes' not in star_dict[KIC]:
        star_dict[KIC]['modes'] = []

    nval = int(column[1])
    lval = int(column[2])
    freq = float(column[3])
    lo = float(column[4])
    hi = float(column[5])
    dfreq = (lo + hi) / 2.
    if lval < 3:
        star_dict[KIC]['modes'].append([lval, nval, freq, dfreq])

sort_list = [(key, star_dict[key]['numax']) for key in star_dict]
sort_list.sort(key=lambda x: x[1])
KIC_list = [pair[0] for pair in sort_list]

leg_tab.close()
leg_numax.close()
leg_file.close()

for index, KIC in enumerate(KIC_list):
    in_file = open("legacy_model_inputs/%d_%d" % (index, KIC), "w")
    for mode in star_dict[KIC]['modes']:
        mode_array = [str(val) for val in mode]
        mode_string = " ".join(mode_array) + "\n"
        in_file.write(mode_string)

    temp_string = "Teff %f %f\n" % (star_dict[KIC]['temp'], star_dict[KIC]['dtemp'])
    numax_string = "numax %f %f\n" % (star_dict[KIC]['numax'], star_dict[KIC]['dnumax'])
    FeH_string = "Fe_H %f %f" % (star_dict[KIC]['metal'], star_dict[KIC]['dmetal'])

    in_file.write(temp_string)
    in_file.write(numax_string)
    in_file.write(FeH_string)
    in_file.close()

temp_array =[star_dict[key]['temp'] for key in star_dict]
temp_array.sort()
print temp_array