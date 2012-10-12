#! /usr/bin/python

######  SUM STATES #######
# Python script for summing and ploting the data from the Density Of States
# files obtained from projwfc.x. It can sum also k-solved dos, and make a plot
# with mathplotlib (if not available, gnuplot, if not avaible, print to file)
# if there is not X11 forwarding, plots in terminal.
# It does something very similar to sumpdos.f90, but with
# some extra features (use "-h" option).
#
# it takes two different inputs, the first one is the pw.x output
# ("-o" option), which is used for parsing the Fermi energy for fitting
# the PDOS curve to the right energy. The other files are the pDOS files
# ("-s" option), that can be given with shell syntax, i.e.
# pdos_atm*Fe*wfc*d* for summing all the d orbitals of Fe.
# It can also handle k solved dos files.
#
# One of the most useful feature, compared to the sumpdos.x, is the
# fact that it also builds the picture directly, so it can be directly
# visualized and exported for inclusion in a document.
# It uses mathplotlib for plotting, but if no mathplotlib is found in
# the $PYTHONPATH, it tries to use gnuplot, if no gnuplot available,
# dumps the output data to a file.
# In the that no X11 forwarding is available (i.e. ssh to the cluster),
# it shows a rough graph in the terminal, so we get an idea of the shape
# of the results.
#
# Example of usage:
# cd ....../espresso-5.0/PP/examples/example02/results/
# ../../../src/sum_states.py -o ni.dos.out -s
# ni.pdos_atm#1\(Ni\)_wfc#2\(d\) -t "Example PP/02" -xr -6 2
#
#
# The procedure for obtaining the DOS files is explained 
# i.e. in (espresso-dir)/PP/examples/example02/
# 
# Author: Dr. Julen Larrucea
#         University of Bremen,
#         Bremen Centre for Computational Materials Science, HMI Group
#         julenl [at] gmail.com  or larrucea [at] hmi.uni-bremen.de
#
#  This file is distributed under the terms of the GNU General Public 
#  License. See the file `License'
#  in the root directory of the present distribution, 
#  or http://www.gnu.org/copyleft/gpl.txt .
#######################

import sys
import os
import fnmatch
import linecache

# Some default variables
version=0.2
pwout=""
selat="*"
graphtitle=""
min_x,max_x=-10,3
min_y,max_y="",""
output_file_name="sum_dos.out"
prt="no"

print " #### sum_states.py version "+str(version)+" #### "  


# Check if X11, mathplotlib and gnuplot are available
try:
 os.popen("gnuplot -V").read()
 prog_gnuplot="yes" # gnuplot is installed
except:
 prog_gnuplot="no"



# Parse command line options
if len(sys.argv)>1:
 for i in sys.argv:
  if i.startswith('-'):
   option=i.split('-')[1]
   if option=="o":
     pwout= sys.argv[sys.argv.index('-o')+1]
   if option=="s":
    selat= sys.argv[sys.argv.index('-s')+1]
   if option=="p":
    prt="yes"
    if len(sys.argv) > sys.argv.index('-p')+1: # if there is a name after "-p" take it as an output name
     if sys.argv[sys.argv.index('-p')+1] != "-": # otherwise default name sum_dos.out
        dos_out_name=sys.argv[sys.argv.index('-p')+1]
   if option=="t":
    graphtitle= sys.argv[sys.argv.index('-t')+1]
   if option=="xr":
    min_x,max_x= float(sys.argv[sys.argv.index('-xr')+1]),float(sys.argv[sys.argv.index('-xr')+2]) 
   if option=="yr":
    min_y,max_y= float(sys.argv[sys.argv.index('-yr')+1]),float(sys.argv[sys.argv.index('-yr')+2]) 
   if option=="v":
    print "sum_dos.py version: "+version
    sys.exit()
   if option=="h":
    print '''
    -o QE output file name (for grepping Fermi E)
    -s Selection of atoms for summing the DOSes. "*" for all, *1*Fe*d* for first Fe atom " (def. "*")
    -p Print output to a file and aditionaly provide an output name (def. no output and "sum_dos.out")
    -t set title in the head of the graph
    -xr set min and max x value for the axes in the graph
    -yr set min and max y value for the axes in the graph
    -h print this help
    -v print version
       Example: sum_states.py --s sys.pdos_atm#4\(Fe2\)_wfc#2\(d\) -t "Wustite LDA+U single Fe" -xr -9 4 
   '''
    sys.exit()


# Check for mathplotlib/gnuplot and import mpl if possible
if len(os.popen('echo $DISPLAY').read()) > 1:
  graphic_plot="yes"
  try:
   from pylab import *
   mplplot="yes"
   print "pylab imported"
  except:
   print "There is no mathplotlib installed. Using gnuplot."
   mplplot="no"
   prt="yes"
else:
  print "No X11. Trying to plot on terminal"
  graphic_plot="no"
  if prog_gnuplot=="no":
   prt="yes"


# if not specified, try to find the espresso output, in order to parse the Fermi energy
if pwout == "":
 for filen in filter(os.path.isfile, os.listdir('.')):
  if "Program PWSCF" in linecache.getline(filen, 2):
   print "Using " + filen + " as pw.x output. You can specify another one with the -o option."
   pwout=filen

# Parse Fermi energy from the pw.x output
if pwout!="":
 try:
  os.popen("grep -a 'the Fermi energy is' "+pwout ).read()
  fermi=float(os.popen("grep -a 'the Fermi energy is' "+pwout ).read().split()[4])
  print "Fermi energy = ", fermi, "a.u."
 except:
  print "WARNING: No Fermi energy found. Using 0 e.V. instead"
  fermi=0
else:
 print "WARNING: No pw.x output found. Using E Fermi = 0 e.V."
 fermi=0

# List of all DOS files to add 
dosfiles=[]
for dfile in os.listdir('.'):
   if fnmatch.fnmatch(dfile, selat):
     dosfiles.append(dfile) 
if len(dosfiles)==0:
 print "ERROR: Provide a (list of) valid DOS file(s)"
 sys.exit()

print "dosfiles list: ",
for dosfile in dosfiles:
  print dosfile,
print ""

# Check wetter we have k-solved DOS
if open(dosfiles[0],'r').readline().split()[1]=="E":
  ksolved="no"
  print "no ksolved"
elif open(dosfiles[0],'r').readline().split()[1]=="ik":
  ksolved="yes"
  print "ksolved"

# Sum over all k-points and files
mat=[]  # matrix with total sum of ldos
for i in range(len(dosfiles)):
 mati=[] # temporal matrix for each DOS file "i"
 k=0
 for line in open(dosfiles[i],'r'):
  if len(line) > 10 and line.split()[0] != "#":

   if ksolved=="no":
      mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])

   if ksolved=="yes":
      ik = int(line.split()[0])
      if ik > k:  #if it is a different k block
         k=int(line.split()[0])
         oldmat=[] # temporal matrix for each k-point
      if ik == 1:
         mati.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]) # append: energy, ldosup, ldosdw
      elif ik == k and k > 1:
         oldmat.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
      elif len(line) < 5 and k > 1:  #if blank line, sum k-frame to the total
       for j in range(len(oldmat)):  
         mati[j]=[mati[j][0],mati[j][1]+oldmat[j][1],mati[j][2]+oldmat[j][2]]

 if mat == []: # if it is the first dos file, copy total matrix (mat) = the first dos files's data
    mat=mati[:]
 else:
    for j in range(len(mati)): # if it is not the first file, sum values
        mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2]]  


print "...ploting..."


if prt=="yes":
 out=open(output_file_name,"w")
x,y1,y2=[],[],[]    
for i in mat:
 x.append(i[0]-fermi)
 y1.append(i[1])
 y2.append(-i[2])
 if prt=="yes":  # print to a file
  print>>out, i[0]-fermi, i[1], i[2]
if prt=="yes":
 out.close()


if graphic_plot=="yes":
  # if there is matplotlib, generate a plot with it
  if mplplot=="yes":
    plot(x,y1,linewidth=1.0)
    plot(x,y2,linewidth=1.0)
    print min(y2),max(y1)
    plt.title(graphtitle)
    plt.xlabel('E (eV)')
    plt.ylabel('States')
    plt.grid(True)
    plt.rcParams.update({'font.size': 22})
    plt.fill(x,y1,color='0.8')
    plt.fill(x,y2,color='0.9')
    if min_x and max_x:
     fromx,tox=min_x,max_x 
    plt.axis([fromx, tox, min(y2), max(y1)])
    show()   
  elif mplplot=="no" and prog_gnuplot=="yes":  # If no mathplotlib available, use gnuplot
     os.system("echo \"plot '"+ output_file_name + "' using ($1-"+str(fermi)+"):2 w l, '' u ($1"+str(fermi)+"):3 w l\" | gnuplot -persist")  
elif graphic_plot=="no":  # If no X forwarding available, show graph in terminal
  if prog_gnuplot=="yes":
     os.system("echo \"set terminal dumb; plot '"+ output_file_name + "' using ($1-"+str(fermi)+"):2 w l, '' u ($1-"+str(fermi)+"):3 w l\" | gnuplot -persist")




