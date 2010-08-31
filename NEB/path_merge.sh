#!/bin/bash
################################################################################
##                    Copyright (C) 2007 Guido Fratesi.                       ##
##                 This file is distributed under the terms                   ##
##                 of the GNU General Public License.                         ##
##                 See http://www.gnu.org/copyleft/gpl.txt .                  ##
##                                                                            ##
##    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         ##
##    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      ##
##    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                   ##
##    NONINFRINGEMENT.  IN NO EVENT SHALL GUIDO FRATESI BE LIABLE FOR ANY     ##
##    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    ##
##    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       ##
##    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  ##
################################################################################
#                                                                              #
# Suitable for PWSCF v. 4.0, might require revisions for future releases.      #
# Execute with no arguments for help.                                          #
#                                                                              #
################################################################################

function showhelp () {
cat <<EOF > /dev/stderr

  USAGE:  path_merge.bash  file0.path n0 m0  file1.path n1 m1

    merges the two path files file0.path and file1.path by taking the
    images n0 to m0 from file0.path and images n1 to m1 from file file1.path.
    Merged file is written to standard output.

EOF

if [ "$err" ]; then cat << EOF > /dev/stderr

* ERROR:  $err

EOF
fi
exit
}

awk=`which awk`

#
# ... Parse and check input parameters
#
inpf[0]=$1   # first path file (file0.path)
imlw[0]=$2   # from this image index (n0)
imup[0]=$3   # to this image index (m0)
inpf[1]=$4   # second path file (file1.path)
imlw[1]=$5   # from this image index (n1)
imup[1]=$6   # to this image index (m1)
#
if (($#==0)); then err=""; showhelp; fi
if (($#<6)); then err="Too few arguments!"; showhelp; fi
if (($#>6)); then err="Too many arguments!"; showhelp; fi
for i in 0 1; do
    if [ ! -f ${inpf[$i]} ]; then
	err="File n.$i \"${inpf[$i]}\" not found!"
	showhelp
    fi
    #
    # ... Extract number of atoms
    #
    nat[$i]=`$awk '
      ($1=="Image:")&&($2==1) {n=NR};
      ($1=="Image:")&&($2==2) {printf NR-n-2; exit};
    ' ${inpf[$i]}`
    #
    # ... Read number of images
    #
    nim[$i]=`$awk '/NUMBER OF IMAGES/ {getline; printf $1; exit}' ${inpf[$i]}`
    #
    if  ((${imlw[$i]}<1)) || ((${imlw[$i]}>${nim[$i]})) ||
	((${imup[$i]}<1)) || ((${imup[$i]}>${nim[$i]})) ||
	((${imup[$i]}<${imlw[$i]})) ; then
	err="Check the images requested from file n.$((i+1))"
	showhelp
    fi
    #
done
#
if ((${nat[0]}!=${nat[1]})); then
    err="The number of atoms in the two path files is not the same!"
    showhelp
fi


#
# ... Write the header
#
cat <<EOF
RESTART INFORMATION
       0
       0
       0
NUMBER OF IMAGES
`printf %4i $(( ${imup[0]}-${imlw[0]}+1 + ${imup[1]}-${imlw[1]}+1 ))`
ENERGIES, POSITIONS AND GRADIENTS
EOF


#
# ... Write the images sections
#
iim=1
for i in 0 1; do
    #
    awk -v iim=$iim -v nat=${nat[$i]} -v imlw=${imlw[$i]} -v imup=${imup[$i]} '
      #
      BEGIN {HUGE=99999; nfix=-HUGE; n=-HUGE};
      #
      # ... Get indexes fixing coordinates
      #
      /ENERGIES, POSITIONS AND GRADIENTS/ {nfix=NR+2};
      (NR>=nfix+1)&&(NR<=nfix+nat) {fix[NR-nfix]=sprintf("%3i%3i%3i",$7,$8,$9)};
      #
      # ... Write images in the range requested
      #
      ($1=="Image:")&&($2>=imlw)&&($2<=imup) { n=NR+1;
                                               printf("Image:%5i\n",iim); }
      (NR==n)
      (NR>=n+1)&&(NR<=n+nat) { for(i=1;i<=6;i++) printf("%20.12f",$i);
                               if(iim==1) printf("%9s",fix[NR-n]);
                               printf("\n"); };
      (NR==n+nat)            { iim++; };
    ' ${inpf[$i]}
    #
    iim=$(($iim+${imup[$i]}-${imlw[$i]}+1))
    #
done

