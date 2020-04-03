#!/usr/bin/python

import numpy as np
import argparse

def is_string(string):
   try:
      float(string)
   except ValueError:
      return True
   return False


def read_currents(filename='corrente_box', dic={}, compute_total=True):
   # open input/output files
   f = open( filename, 'r' )
   for line in f:
      values = line.split()
      if is_string(values[0]):
         if values[0] in dic:      # key already exists: append line
            dic[values[0]] = np.vstack([dic[values[0]], map(float, values[1:])])
         else:
            dic[values[0]] = np.array(map(float, values[1:]))
      else:
         pass     # not a string -> do nothing
   if compute_total:
      dic['total'] = compute_total_current( dic )
   return dic


def compute_total_current(dic):
   # compute total current
   if (len(dic['h&K-K'][1:4].shape) > 1):
      return dic['h&K-XC'] + dic['ionic:'] + dic['h&K-K'][:,1:4] + dic['h&K-H'] + dic['zero:']
   else:
      return dic['h&K-XC'] + dic['ionic:'] + dic['h&K-K'][1:4] + dic['h&K-H'] + dic['zero:']


def print_total_current( dic ):
   print("total:  {:.12E} {:.12E} {:.12E}".format(dic['total'][0], dic['total'][1], dic['total'][2]))
   return


def main():
   """This programs reads the heat current components, stores them in a dictionary and compute the total heat current."""
   parser = argparse.ArgumentParser()
   parser.add_argument('inputfile', help='heat current file to read')
   args = parser.parse_args()

   inputfile = args.inputfile

   # open input/output files
   cur = read_currents(inputfile)

   # print total current
   print_total_current(cur)

   return 0


if __name__ == "__main__":
   main()
