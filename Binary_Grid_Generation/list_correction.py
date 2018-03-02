# Generates the full file locations for the final output table and saves this
# out to the desired location.

########################################

def get_lines(filename):
	lines = []
	with open(filename) as f:
		for line in f:
			lines.append(line)

	return lines

##########################################

input_file = "master"	# file containing the initial table
output_file = "/home/bmr135/AIMS_New/CLES_MS_DIFF2_list"	# output file
lines = get_lines(input_file)


names = []
directs = []
to_freqs = []
for line in lines:
	names.append(line[39:])
	directs.append(line[:25])
	to_freqs.append(line[:34])

commands = []
for i in range(len(names)):
	# stitch together sections of the file name to create the desired file location
	commands.append(directs[i] + "/AIMS/" + to_freqs[i] + ' ' + names[i]) # CLES
	# commands.append(directs[i] + "/" + names[i]) # MESA

with open(output_file,"w") as f:
	# write data to output file with extension to location of frequency files included
	f.write("/home/bmr135/GridSunClesDV2/models/	.freq \n") # CLES
	# f.write("/home/miglioa/GridNGC6819_Fep0.25/LOGS/	.sgyre_l0 \n") # MESA
	for command in commands:
		f.write(command)
