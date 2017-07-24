import os
import subprocess

file_list = os.listdir('.')
model_list = [fl for fl in file_list if 'Model_' in fl and fl != 'Model_params']

go = raw_input("Models found = %d, do you wish to continue? [y/n]: " % len(model_list))
if go == "y":
    for i in xrange(len(model_list)):
        model_name = "Model_%d" % i
        print "Calculating %s\n" % model_name
        p = subprocess.Popen(["python", "AIMS.py", model_name])
        p.wait()
        print "\nCompleted calculating %s\n\n\n\n" % model_name
