import os
import subprocess

folders = ['50', '75', '100']
commands = []

delay = []
jitter = []
variation = []
for f in folders:
	for i in range(10, 110, 10):
		path = "Washington-" + f + "/" + str(i)
		instance = os.listdir(path)
		for inst in instance:
			if (inst[0] == 'w' or inst[0] == '1'):
                                for i in range(1, 5, 1):
                                        commands.append("./MaxServiceRL " + str(i) + " " + path + "/" + inst + " " + path + "/param-" + inst + " results_" + str(i) + "_1" + "/result_" + inst + " 1 1 1 1 1.5 1000 50 1800")
                                        commands.append("./MaxServiceRL " + str(i) + " " + path + "/" + inst + " " + path + "/param-" + inst + " results_" + str(i) + "_0" + "/result_" + inst + " 1 1 1 0 1.5 1000 50 1800")

for c in commands:
   print(c)
   p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
   msg, err = p.communicate()
   if msg:
   	print(msg)
   print("OK!!")
