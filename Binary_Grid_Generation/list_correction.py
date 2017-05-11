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
output_file = "/home/miglioa/GridCLESAIMS_D/CLES_MS_DIFF"	# output file
lines = get_lines(input_file)


names = []
directs = []
for line in lines:
	names.append(line[:])
	directs.append(line[:20])

commands = []
for i in range(len(names)):
	# stitch together sections of the file name to create the desired file location
	commands.append(directs[i] + "/AIMS/" + names[i])

with open(output_file,"w") as f:
	# write data to output file with extension to location of frequency files included
	f.write("/home/miglioa/GridCLESAIMS_D/	.freq \n")
	for command in commands:
		f.write(command)
