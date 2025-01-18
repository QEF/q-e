#!/usr/bin/env python3

# Remove unused module USE statements
# 1. Edit QE/make.inc to  "FFLAGS = -O0 -g  -Wall -fbounds-check -frange-check" (gfortran)
# 2. In EPW/src, run make clean; make epw > out 2>&1.
# 3. Run ./utilities/clean_unused_variables.py out
# 4. Repeat 2&3 several times.

# Limitations
# 1. If there are multiple USE of the same module, some of them will not be deleted.
# 2. Whitespace needs to be manually fixed

import os
import sys

with open(sys.argv[1], "r") as f:
    lines = f.readlines()


unused_variables = []
for i, line in enumerate(lines):
    if "Unused module variable" in line:
        variable = line.strip().split()[4][1:-1]
        filename = lines[i-4].split(":")[0]
        iline = int(lines[i-2].strip().split()[0])
        unused_variables += [(filename, variable, iline)]

# Sort to have lines in decreasing order
def myKey(x):
    filename, variable, iline = x
    return (filename, -iline)
unused_variables.sort(key = myKey)

for filename, variable, iline_to_delete in unused_variables:

    if filename == "qdabs.f90":
        continue

    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename + ".tmp", "w") as f:
        remove_this_line = False
        removed_variable = False
        iline_continued = -1

        for iline, line in enumerate(lines):

            # If this line is a continuation of a previous line, do not update iline_continued.
            # Otherwise, update iline_continued.
            if iline > 0 and len(lines[iline-1].strip()) > 0 and lines[iline-1].strip()[-1] == "&":
                pass
            else:
                iline_continued = iline

            if variable in line.lower() and iline_continued + 1 == iline_to_delete:
                # print(line, iline, iline_continued, iline_to_delete)

                if "," in line.split(":")[-1]:
                    # Multiplie variables declared in one line, remove only "variable"

                    # Add everything before ":"
                    if ":" in line:
                        part1, part2 = line.split(":")
                        line_new = part1 + ":"
                    else:
                        part2 = line
                        line_new = ""

                    count = 0
                    for item in part2.split(","):
                        if item.strip().lower() == variable:
                            # Variable to be removed
                            continue

                        else:
                            # Variable to be kept
                            if count == 0:
                                line_new += item.strip('\n')
                            else:
                                line_new += "," + item.strip('\n')
                            count += 1

                    line_new += "\n"
                    print(f"MODIFY {filename}:{iline+1} : {line_new.strip()}")
                    f.write(line_new)
                    remove_this_line = True

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


print()
print("=" * 50)
print("Variables to be manually deleted")
for filename, variable, iline_to_delete in unused_variables:
    print(f"{filename:20s}:{iline_to_delete:4d} {variable}")