# This is for making the list needed to make the binary grid for AIMS.
# Input is a list of structures and the paths to them.
# Can put the header line on manually

import numpy as np
import pandas as pd
import sys

########################################

def get_profiles(filename):
	names = []
	with open(filename) as f:
		for line in f:
			names.append(line.rstrip("\n"))

	return names

def get_others(profile,data_list):
	# Inputs are name of profile file and list to put star's info into
	# Construct name of track file using corresponding profile
	model = int(profile[-8:-4]) # CLES
	# if len(profile) == 59:
	# 	model = int(profile[42:45]) # MESA
	# elif len(profile) == 58:
	# 	model = int(profile[42:44]) # MESA
	# if len(profile) == 57:
	# 	model = int(profile[42]) # MESA
	prefix = profile[:20] # CLES
	# prefix = profile[:40] # MESA
	suffix = "-sum.txt" # CLES
	# suffix = ".track" # MESA
	track = "/home/bmr135/GridGiantsClesV0.3/models_grad_rad_under/tracks/" + prefix + suffix
	# track = "/home/miglioa/GridNGC6819_Fep0.25/LOGS/" + prefix + "/" + prefix + suffix
	# print track

	''' MESA index conversion to match model numbers input to those on the tracks '''
	# index = pd.read_csv('/home/bmr135/GridSunClesDV2/models/' + prefix + '/' + prefix + '.index',skiprows='1',names=['in','pri','out'],delimiter=r'\s+')
	# for i in range(len(index)):
	# 	if model == index['in'][i]:
	# 		model = index['out']

	# Get data from track
	data = []
	with open(track) as f:
		for line in f:
			data.append(line.rstrip("\n"))
	data = data[3:] # Chop off headers, CLES
	# data = data[6:] # Chop off headers, MESA
	for i in range(len(data)):
		data[i] = data[i].replace("D","e")
		data[i] = data[i].split()
		for j in range(len(data[i])):
			data[i][j] = float(data[i][j])
		if int(data[i][0]) == model:
			model_data = data[i]
			break
	# print model_data
	''' If models are in an ordered index, use the code below '''
	# first_model = data[0][0]
	# # Get line corresponding to model and get its components
	# model_data = data[model-first_model]

	data_list[0] = model_data[19] # mass in g
	data_list[1] = model_data[9] * 6.957e10 # radius in cm, CLES
	# data_list[1] = model_data[7] * 6.95508e10 # radius in cm, MESA
	data_list[2] = 10.0**(model_data[3]) * 3.828e33 # luminosity in erg/s, CLES
	# data_list[2] = model_data[6] * 3.8422e33 # luminosity in erg/s, MESA
	data_list[3] = model_data[17] # Z mass fraction, CLES
	data_list[4] = model_data[18] # hydrogen mass fraction, CLES
	data_list[5] = model_data[1] * 1e-6 # age in Myr, CLES
	# data_list[5] = model_data[1] * 1e-6 # age in Myr, MESA
	data_list[6] = 10.0**(model_data[2]) # effective temp in K, CLES
	# data_list[6] = model_data[3] # effective temp in K, MESA
	data_list[7] = model_data[22] * 1.988475415e33 # He-core mass, CLES
	data_list[8] = model_data[4] # central hydrogen fraction, CLES
	# data_list[8] = model_data[15] # central hydrogen fraction, MESA
	data_list[9] = model_data[20] # period spacing (toggle depending upon grid), CLES
	# data_list[8] = model_data[69] # period spacing (toggle depending upon grid), MESA




########################################

input_file = "structure.lis"
output_file = "master"
profiles = get_profiles(input_file)


star_info = [None for i in range(10)]	#change to 9 to include DNl1

with open(output_file,"w") as f:
	for i in range(len(profiles)):
		f.write(profiles[i] + "\t")
		get_others(profiles[i],star_info)
		for j in range(len(star_info)-1):
			f.write(str(star_info[j]) + "\t")
		f.write(str(star_info[9]) + "\n")	#Change to 8 to include DNl1
