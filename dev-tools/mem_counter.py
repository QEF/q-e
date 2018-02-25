"""
This simple script adds a CALL statement (or, more precisely, whatever 
statement you want to add) after each call to ALLOCATE/DEALLOCATE.

Usage:
    python instrument.py file_name

this command will create a bkp file (if not present) and a .new file
which contains the print statements.

something like:
    files=`sift -i --filename -x f90 -w *allocate | cut -d ':' -f 1 | uniq`
    for f in $files
    do
      python instrument.py $f
      mv ${f}.new $f
    done

will add print statements to all the files in the directory

restore with

for file in *.bkp; do cp $file ${file%.bkp}; done

"""

import os, sys, re
from shutil import copyfile

fmt = """ CALL mem_counter ( SIZEOF({0}), {1}, "{0}" )
"""

if len(sys.argv) < 2:
  print("Specify file name!")
  sys.exit(-1)

fname = sys.argv[1]

if not os.path.exists(fname):
  print("File does not exists!")
  sys.exit(-1)
  

#if not os.path.exists(fname+".bkp"):
#  copyfile(fname, fname+".bkp")


def get_args_list(line):
  args = ""
  npar = 0
  for i, c in enumerate(line):
    if c == '(':
      npar +=1
      continue
    if c == ')':
      npar -=1
      continue
    if c == '&':
      continue
    if c == '\n':
      continue
    if c == '!':
      break
    if c == '%':
      # copy argument if needed
      if npar >= 2:
        continue
      if line[i-1] == ')':
        j = i
        while (j>0):
          j = j-1 # this can fail but should never happen
          if line[j] == '(':
            break
        
        args += line[j:i].replace(',','!') # temp replace , with ! to 
                                           #avoid confusion in arg splitting.
                                           # btw, we go back in characters so
                                           # there shouldn't be any '!' 
           
    if npar == 1:
      args += c
    
  return [x.replace('!',',') for x in args.split(',')]


def write_counters(fout, sign, args, tag=[], prepend=""):
  for a, arg in enumerate(args):
    if 'stat' in arg.lower():
      continue
    if (a == 0) and (tag != []):
      fout.write(tag[0] + ' ')
    fout.write(prepend + fmt.format(*[arg, sign]))

with open(fname,'r') as fin:
  with open(fname+".new",'w') as fout:
    lines = fin.readlines()
    i = 0
    
    while (i<len(lines)):
      if len(lines[i].strip()) > 0:
        
        # write out comments straight away
        if ('!' in lines[i].strip()[0]):
          fout.write(lines[i])
          i += 1
          continue
        
        # possibly get rid of comments
        buff = lines[i].split('!')[0]
        
        # add end of line if needed
        if buff[-1] != '\n':
          buff += '\n'
        
        # merge continuation in single string. Never forget \n
        j = 0
        while '&' in lines[i+j].strip()[1:]: # line continuation may begin with '&'
          j += 1
          buff += lines[i+j].split('!')[0]
          if buff[-1] != '\n':
            buff += '\n'

        # find (de)allocate in strings
        str_allocate = re.findall("[\"\'].*?(allocate).*?[\"\']", buff.lower())
        
        # match all allocate
        allocate = re.findall(r'\bdeallocate\b', buff.lower()) + re.findall(r'\ballocate\b', buff.lower())
        
        # if no allocate or deallocate present, write all previous lines as they were originally
        if ((allocate == []) or (len(str_allocate) >= len(allocate))):
          for _ in range(j+1):
            fout.write(lines[i])
            i += 1
            
          continue

        
        if re.findall(r'\bdeallocate\b', buff.lower()) and re.findall(r'\bif\b', buff.lower()):
          
          # split clause with 'deallocate' as delimiter
          ifcls = re.split(r'\bdeallocate\b',buff.lower())
          
          args = get_args_list(ifcls[1])

          if (len(args) > 1):
            print("WARNING! Multiple deallocate in if clause")
          
          tag = re.findall(r'^\d+', buff.lower().strip())
          
          write_counters(fout, -1, args, tag, ifcls[0] + ' ')
            
          #fout.write(ifcls[0] + ' ' + fmt.format(*[args[0], -1]))
          if tag:
            fout.write(re.sub("^[ X]*\d+","",buff))
          else:
            fout.write(buff)
          
        elif re.findall(r'\ballocate\b', buff.lower()) and re.findall(r'\bif\b', buff.lower()):
          ifcls = re.split(r'\ballocate\b',buff.lower())
          args = get_args_list(ifcls[1])
          if (len(args) > 1):
            print("WARNING! Multiple allocate in if clause")

          # is there a tag?
          tag = re.findall(r'^\d+', buff.lower().strip())
          if tag:
            print("WARNING! deallocate with tag! Check (and change?) file " + fname + ".new (original line " + str(i) + ")")
          
          
          fout.write(ifcls[0] + ' THEN \n ALLOCATE ' + ifcls[1] + '\n ')
          write_counters(fout, +1, args)
          fout.write( '\n END IF \n')
          
        elif re.findall(r'\bdeallocate\b', buff.lower()):
          
          args = get_args_list(buff)

          # is there a tag?
          tag = re.findall(r'^\d+', buff.lower().strip())
          if tag:
            print("WARNING! deallocate with tag! Check (and change?) file " + fname + ".new (original line " + str(i) + ")")
          
          
          write_counters(fout, -1, args, tag)
          #for a, arg in enumerate(args):
          #  if 'stat' in arg.lower():
          #    continue
          #  if (a == 0) and (tag != []):
          #    fout.write(tag[0] + ' ')
          #  fout.write(fmt.format(*[arg, -1]))
            
          if tag:
            fout.write(re.sub("^[ X]*\d+","",buff))
          else:
            fout.write(buff)
          
        elif re.findall(r'\ballocate\b', buff.lower()):
          fout.write(buff)
          args = get_args_list(buff)

          # is there a tag?
          tag = re.findall(r'^\d+', buff.lower().strip())
          if tag:
            print("WARNING! allocate with tag! Check (and change!) file " + fname + ".new (original line " + str(i) + ")")

          write_counters(fout, +1, args)
          #for arg in args:
          #  if 'stat' in arg.lower():
          #    continue
          #
          #  fout.write(fmt.format(*[arg, 1]))
        else:
          fout.write(buff)
        
        i += 1 + j
        continue
      
      fout.write(lines[i])
      i += 1
  
    
  
