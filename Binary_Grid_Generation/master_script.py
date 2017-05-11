# Makes script that extracts the tgz files

########################################

def get_lines(filename):
	lines = []
	with open(filename) as f:
		for line in f:
			lines.append(line.rstrip("\n"))

	return lines

##########################################

root = "/home/miglioa/GridCLESAIMS_D"	# location of .tgz folders
input_file = "files"
output_file = "extract"
lines = get_lines(input_file)

names = []
directs = []
for line in lines:
	names.append(line[2:])
	directs.append(root)

commands = []
for i in range(len(names)):
	# command line to create extraction code for each input file
	commands.append("cd " + directs[i] + "\ntar xzvf " + names[i] + "\n")

with open(output_file,"w") as f:
	f.write("#!/bin/bash\n\n")
	for command in commands:
		f.write(command)
