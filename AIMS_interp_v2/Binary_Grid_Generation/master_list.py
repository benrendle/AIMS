# This is for making the list needed to make the binary grid for AIMS.
# Input is a list of structures and the paths to them.
# Can put the header line on manually

import numpy as np

########################################

def get_profiles(filename):
	names = []
	with open(filename) as f:
		for line in f:
			names.append(line.rstrip("\n"))

	return names



def get_MXZ(profile,data_list):
	# Inputs are name of profile file and list to put star's info into
	data_list[0] = float(profile[1:5]) * 1.9884e33 # mass in g
	data_list[4] = float(profile[7:12]) # hydrogen mass fraction
	data_list[3] = float(profile[14:20]) # metallicity



def get_others(profile,data_list):
	# Inputs are name of profile file and list to put star's info into
	# Construct name of track file using corresponding profile
	model = int(profile[-12:-8])
	prefix = profile[:20]
	suffix = "-sumN.txt"
	track = "/home/miglioa/GridCLESAIMS_D/" + prefix + "/" + prefix + suffix

	# Get data from track
	data = []
	with open(track) as f:
		for line in f:
			data.append(line.rstrip("\n"))
	data = data[3:] # Chop off headers
	for i in range(len(data)):
		data[i] = data[i].replace("D","e")
	data[0] = data[0].split()
	first_model = int(data[0][0])
	# Get line corresponding to model and get its components
	model_data = data[model-first_model].split()
	for i in range(len(model_data)):
		model_data[i] = float(model_data[i])

	data_list[1] = model_data[9] * 6.95508e10 # radius in cm
	data_list[2] = 10.0**(model_data[3]) * 3.8422e33 # luminosity in erg/s
	data_list[5] = model_data[1] * 1e-6 # age in Myr
	data_list[6] = 10.0**(model_data[2]) # effective temp in K
	data_list[7] = model_data[4] # central hydrogen fraction
	# data_list[8] = model_data[20] # period spacing (toggle depending upon grid)



########################################

input_file = "structure.lis"
output_file = "master"
profiles = get_profiles(input_file)


star_info = [None for i in range(8)]	#change to 9 to include DNl1

with open(output_file,"w") as f:
	for i in range(len(profiles)):
		f.write(profiles[i] + "\t")
		get_MXZ(profiles[i],star_info)
		get_others(profiles[i],star_info)
		for j in range(len(star_info)-1):
			f.write(str(star_info[j]) + "\t")
		f.write(str(star_info[7]) + "\n")	#Change to 8 to include DNl1
