# RSF_solvers

some notes...


To compile the code, from the directory that has nrsfs.c and rsp.h:
gcc -o nrsfs nrsfs.c -lm

Then make sure the executable is located somewhere in your path (or from the directory where it lives just type ./nrsfs)


Workflow:
write the script, running things a few times to adjust the parameters if you don't know what they should be (my script is called runs.VS_ss)
	Note that my script provides a way to include the labels (via echo) in output because those normally just to go the terminal (STDerr)
make the script executable: e.g., chmod +x scriptName  (for my script this is chmod +x runs.VS_ss)
run it and redirect the output to a file for plotting (for my script this is: runs.VS_ss > runs.output)

Note that my code has an explicit parameter for the output step size. This makes the file size manageable but it means that max/min values are often missed. The code accounts for that by providing a table of important paramters with exact values.

Here's a little command line python plotter so that you can make quick plots from unix
msp.py
I suggest to set up an alias such as this: 
alias msp 'python ~/bin/msp.py'
then you can simply type:  msp filename.dis, etc. (for my filenames, e.g., msp runs1.dis or msp runs1.vel)
