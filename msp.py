#python version of msp and xmsp
#cjm 20230618
# superuseful: https://realpython.com/python38-new-features/
# look at argparse

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

#print(f"Arguments count: {len(sys.argv)}")
argcount = len(sys.argv)
progname = sys.argv[0]

#print("argcount = ", argcount)
#print(f"Name of the script      : {sys.argv[0]=}")
#print(f"Arguments of the script : {sys.argv[1:]=}")

#filename = sys.argv[1]

if argcount < 2 :
	print("This python script plots output from the program 'rsfs'. You can plot multiple files" )
	print("Usage: ", progname, "filename, filename, filename, -n")
	print("files are assumed to have a 2 line header, with column name etc in line 1 and line 2 blank. Use the -n option for no header\n")
	exit()

mycolor = iter(cm.rainbow(np.linspace(0, 1, argcount)))

i=1
header = True
option = argcount - 1
flag = sys.argv[option]
if flag == "-n":
	argcount -= 1
	header = False 
	print("Using -n option")
	print("argcount =", argcount, "flag=", flag)

while(i < argcount):

	filename = sys.argv[i]
	print(filename, "argcount =", argcount, "i=", i)

	dispOrtime = []
	friction = []

	infile = open(filename,"r")
	lines = infile.readlines()

	print(lines[0])
	line = lines[0]

	if header == True:
		numbers = line.split(",")
		xaxis_label = numbers[0]
		yaxis_label = numbers[1]
		plot_label = numbers[2]
		
		for line in lines[2:]: #read header on line 1, data start at line 2
    			numbers = line.split("\t")
    			dispOrtime.append(float(numbers[0]))
    			friction.append(float(numbers[1]))
	else:
		xaxis_label = "x" 
		yaxis_label = "y" 

		for line in lines[0:]: #no header, data start at line 0
    			numbers = line.split("\t")
    			dispOrtime.append(float(numbers[0]))
    			friction.append(float(numbers[1]))

	infile.close()
 
	if(i == 1):
		fig = plt.figure()
		ax = plt.subplot(111)

	#	color = 'r'
	#if(i == 2):
	#	color = 'g'

	color = next(mycolor)	

	if header == True:
		ax.plot(dispOrtime,friction, c = color, label = plot_label)
		dist_above = 1.1 + i*0.05
		ax.legend(loc='upper center', bbox_to_anchor=(0.7, dist_above),
          	ncol=1, fancybox=False, shadow=True)
	else:
		ax.plot(dispOrtime,friction, c = color)

	plt.xlabel(xaxis_label, fontsize = 12)
	plt.ylabel(yaxis_label, fontsize = 12)

	i += 1 

plt.show(block=True)

exit()

