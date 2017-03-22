#!/usr/bin/python
#
# Copyright (C) 2017 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# Written by Lorenzo Paulatto 2017
#
from os.path import realpath, basename
from sys import argv
myname = basename(argv[0])
#ibrav2cell_x = realpath(__file__)
#ibrav2cell_x = ibrav2cell_x.replace(myname,"ibrav2cell.x")
ibrav2cell_x = argv[0].replace(myname,"ibrav2cell.x")

def main() :
  from numpy import array
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('-i', action='append', type=int, metavar="ibrav", dest="ibrav_list", default=[], 
                       help="ibrav value to explore, can be repeated (default: all)")
  parser.add_argument('-c', type=float, nargs=9, dest="at",
                      metavar=("a1x","a1y","a1z","a2x","a2y","a2z", "a3x", "a3y", "a3z"), 
                      help="cell parameters (default: read from standard input)")
  units = parser.add_mutually_exclusive_group()
  units.add_argument('-A', default=False, action="store_true", dest="angst", help="cell is entered in Agstrom units")
  units.add_argument('-B', default=True, action="store_false", dest="angst", help="cell is entered in Bohr units (default)")
  units.add_argument('--alat', type=float, dest="alat", help="cell is entered in units of alat (in bohr)", metavar="ALAT")

  parser.add_argument('--celldm', type=float, default=[1, 1, 1, 0.5, 0.5, 0.5], nargs=6, dest="celldm",
                      metavar=("1","2","3","4","5","6"), help="initial values of celldm(1..6), you have to specify all 6 or none (default: 1 1 1 .5 .5 .5)")
  parser.add_argument('-t', type=float, default=1.e-3,  dest="mthr", help="match threshold", metavar="THR")
  parser.add_argument('-k', type=float, default=1.e-15, dest="kthr", help="convergence threshold", metavar="THR")
  parser.add_argument('-x', type=str, dest="ibrav2cell_x", help="full path to the ibrav2cell.x tool from PW/tools/ \
                      (default: look in the same directory as this program)", metavar="/path/to/ibrav2cell.x")

  args = parser.parse_args()

  ibrav_list = [1,2,3,4,5,-5,6,7,8,9,-9,10,11,12,-12,13,14]
  if args.ibrav_list == []:
    args.ibrav_list = ibrav_list

  global ibrav2cell_x
  if args.ibrav2cell_x:
    ibrav2cell_x = args.ibrav2cell_x
    
  check_ibrav2cell()

# test for ibrav = 5
#  cell = [    3.900896593574796,   -2.252183691403927,   18.043044401344225,
#              0.000000000000000,    4.504367382807855,   18.043044401344225,
#             -3.900896593574796,   -2.252183691403927,   18.043044401344225]
# test for ibrav = 14
#  cell = [ 1.200000000000000E+01,   0.000000000000000E+00,   0.000000000000000E+00,
#           1.224000000000000E+01,   7.585670702053973E+00,   0.000000000000000E+00,
#           7.979999999999999E+00,   1.474991525399383E+00,   2.168870673876157E+00]

  if args.alat and args.angst:
    print "You cannot specifiy both '-A' and '--alat'"
    exit(254)

  cell = args.at
  if not cell:
    cell = cell_stdin()

  if args.angst:
    cell = array(cell)*1.889725989
  elif args.alat:
    cell = array(cell)*args.alat

  bnds = ((0,None), (0,None), (0,None), (-1,1), (-1,1), (-1,1))
  
  print "Scanning..."
  for ibrav in args.ibrav_list:
      p=args.celldm 
      options = { "ibrav" : ibrav,
                  "cell"  : cell
                }
      #from scipy.optimize import basinhopping
      #r = basinhopping(recompute, p, minimizer_kwargs={"method":"L-BFGS-B", "bounds":bnds, "args":tuple([options])}, niter=10)
      from scipy.optimize import minimize
      r = minimize(recompute, p, args=tuple([options]), method="L-BFGS-B", tol=args.kthr, bounds=bnds)
      if r["fun"] < args.mthr:
        print "match found: ibrav =", ibrav 
        check_zero(r["x"], ibrav)
      elif r["fun"] < 100*args.mthr:
        print "possible noisy match found: ibrav =", ibrav 
        check_zero(r["x"], ibrav)

def check_ibrav2cell():
  from numpy import array
  from sys import argv
  if not is_exe(ibrav2cell_x):
    print " File '"+ibrav2cell_x+"' not found or not executable."
    print " Specify the position of ibrav2cell.x with the '-x' command line option." 
    print " Use '-h' for more help."
    exit(100)
  namelist = make_namelist(1,[1,0,0,0,0,0])
  at = compute_cell(namelist)
  if not at or len(at)!=9 or any(at != array([1,0,0,0,1,0,0,0,1])):
    print " File '"+ibrav2cell_x+"' not working as expected."
    print " Please verify that it comes from the same distribution as '"+argv[0]+"'"
    print " and that it is running correctly."
    exit(101)

def is_exe(fpath):
    import os
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)  

def run_command(executable_file, input_data):
    import subprocess
    proc = subprocess.Popen(executable_file,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            )
    output, error = proc.communicate(input_data)
    return (output, error)

def make_namelist(ibrav,p):
  replacements  = { "ibrav" : ibrav,
                   "celldm1" : p[0],
                   "celldm2" : p[1],
                   "celldm3" : p[2],
                   "celldm4" : p[3],
                   "celldm5" : p[4],
                   "celldm6" : p[5]
                   }
  systemnl = """
  &system
    ibrav = {ibrav}, 
    celldm(1) ={celldm1},
    celldm(2) ={celldm2},
    celldm(3) ={celldm3},
    celldm(4) ={celldm4},
    celldm(5) ={celldm5},
    celldm(6) ={celldm6}
  /
  """
  return systemnl.format(**replacements)

def constraint(p):
  result = all(p[0:3]>=0)
  result = result and all(abs(p[4:7])<=1)
  if result:
    return 0
  else:
    return 1

def cell_stdin():
  #import fileinput
  from sys import argv, stdin
  print "Type '{} -h' for help".format(argv[0])
  print
  print "Please enter cell parameters: "
  at = []
  while len(at)<9:
    line = stdin.readline()
    for w in line.split():
      at.append(float(w))
  print " ------------------------ ok"
  return at

def compute_cell(namelist):
  from numpy import isnan
  output,error = run_command(ibrav2cell_x,namelist)
  try:
    at = map(float, output.split())
  except ValueError:
    at = [0,0,0, 0,0,0, 0,0,0]
  if any(isnan(at)):
    at = [0,0,0, 0,0,0, 0,0,0]
  return at

def check_zero(p, ibrav):
  from numpy import array
  from numpy.linalg import norm
  namelist = make_namelist(ibrav,p)
  at0 = array(compute_cell(namelist))
  p_out=list(p)
  for i in range(len(p)):
    q = list(p)
    q[i] = 0
    namelist = make_namelist(ibrav,q)
    at1 = array(compute_cell(namelist))
    if norm(at0-at1)==0:
      p_out[i] = 0
    else:
      print "    celldm({:d}) = {:.6f}".format(i+1, p[i])
  return array(p_out)

def recompute(p, options):
  from numpy import array
  from numpy.linalg import norm
  namelist = make_namelist(options["ibrav"],p)
  at = compute_cell(namelist)
  at1 = array(at)
  at0 = array(options["cell"])
  return norm(at1-at0)

main()



