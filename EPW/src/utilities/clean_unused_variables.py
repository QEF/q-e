#!/usr/bin/env python3

# Remove unused variables
# 1. Edit QE/make.inc to  "FFLAGS = -O0 -g  -Wall -fbounds-check -frange-check" (gfortran)
# 2. In EPW/src, run make clean; make epw > out 2>&1.
# 3. Run ./utilities/clean_unused_variables.py out
# 4. Repeat 2&3 several times.

# Removes unused variable and the following FORD comment (starting with !!).

# Limitations
# 1. If there is comma in variable declaration (multiple variables declared, array with dim > 1)
#    the line will not be deleted.


import os
import sys

with open(sys.argv[1], "r") as f:
    lines = f.readlines()


unused_variables = []
for i, line in enumerate(lines):
    if "Unused variable" in line:
        variable = line.strip().split()[3][1:-1]
        filename = lines[i-4].split(":")[0]
        iline = int(lines[i-2].strip().split()[0])
        unused_variables += [(filename, variable, iline)]

# Sort to have lines in decreasing order
def myKey(x):
    filename, variable, iline = x
    return (filename, -iline)
unused_variables.sort(key = myKey)


for filename, variable, iline_to_delete in unused_variables:

    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename + ".tmp", "w") as f:
        remove_this_line = False
        removed_variable = False

        for iline, line in enumerate(lines):
            if variable in line.lower() and "::" in line and iline+1 == iline_to_delete:
                if "," in line:
                    # Multiplie variables declared in one line
                    print(f"### File {filename} line {iline+1} contains unused variable {variable} but not removed.")
                else:
                    # Declaration of unused variable
                    remove_this_line = True
                    removed_variable = True
                    print(f"REMOVE {filename}:{iline+1} : {line.strip()}")

            elif removed_variable and len(line.strip()) >= 2 and line.strip()[0:2] == "!!":
                # Comment following declaration of unused variable
                remove_this_line = True
                removed_variable = True
                print(f"REMOVE {filename}:{iline+1} : {line.strip()}")

            else:
                # Do nothing
                remove_this_line = False
                removed_variable = False
            
            if not remove_this_line:
                f.write(line)

    os.rename(filename + ".tmp", filename)


for filename, variable, iline_to_delete in unused_variables:
    print(f"{filename:20s}:{iline_to_delete:4d} {variable}")