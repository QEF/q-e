#!/bin/bash --noprofile
#
################################################################################
##                    Copyright (C) 2004 Carlo Sbraccia.                      ##
##                 This file is distributed under the terms                   ##
##                 of the GNU General Public License.                         ##
##                 See http://www.gnu.org/copyleft/gpl.txt .                  ##
##                                                                            ##
##    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         ##
##    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      ##
##    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                   ##
##    NONINFRINGEMENT.  IN NO EVENT SHALL CARLO SBRACCIA BE LIABLE FOR ANY    ##
##    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    ##
##    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       ##
##    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  ##
################################################################################
#
################################################################################
##                Set these variables according to your needs                 ##
################################################################################
#
# ... root directory of the PWscf package
#
ROOT_DIR=""
#
# ... old and new restart file name 
# ... ( usually old_prefix.neb and new_prefix.neb )
#
old_restart_file=""
new_restart_file=""
#
# ... number of images of the old and of the new path
#
old_num_of_images=
new_num_of_images=
#
# ... interpolation is performed between the first_image and the last_image
# ... of the old path
#
first_image=
last_image=  
#
# ... lattice parameter ( celldm(1) in the pw input file )
#
alat=1.D0
#
# ... interpolation is possible only in the ibrav=0 case
# ... this card is the same used in the pw input file
#
cat > CELL_PARAMETERS << EOF
CELL_PARAMETERS
  12.00000  0.00000  0.00000
   0.00000  5.00000  0.00000
   0.00000  0.00000  5.00000
EOF
#
# ... needed to gereate visualization files ( xyz and axsf formats ) 
#
# ... put a symbol for each atomic specie that compose the system
#
list_of_atoms=""
#
# ... number of atoms of each specie
#
N[1]=
#
################################################################################
################################################################################
######             DO NOT MODIFY THE SCRIPT UNDER THESE LINES             ######
################################################################################
################################################################################
#
GAWK=$( which gawk )
#
################################################################################
##                    the total number of atoms is computed                   ##
################################################################################
i=1
Ntot=0
while [ ${i} -le ${#N[@]} ]; do
  Ntot=$(( ${Ntot} + ${N[${i}]} ))
  i=$(( ${i} + 1 ))
done
#
################################################################################
##                     lattice vectors for gawk scripts                       ##
################################################################################
a11=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==2) {printf "%12.8f", $1} }' )
a12=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==2) {printf "%12.8f", $2} }' )
a13=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==2) {printf "%12.8f", $3} }' )
a21=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==3) {printf "%12.8f", $1} }' )
a22=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==3) {printf "%12.8f", $2} }' )
a23=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==3) {printf "%12.8f", $3} }' )
a31=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==4) {printf "%12.8f", $1} }' )
a32=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==4) {printf "%12.8f", $2} }' )
a33=$( cat CELL_PARAMETERS | $GAWK '{ if (NR==4) {printf "%12.8f", $3} }' )
#
################################################################################
##            the input file for the iterpolator code is generated            ##
################################################################################
#
if [ ! -f ${old_restart_file} ]; then
  echo "Error: file ${old_restart_file} not fount"; exit
fi
#
cat > input << EOF
${Ntot}
${old_num_of_images}
${new_num_of_images}
${first_image}
${last_image}
${old_restart_file}
${new_restart_file}
${alat}
EOF
#
cat CELL_PARAMETERS | $GAWK '{ if ( NR == 1 ) { print }; if ( NR > 1 ) \
{ printf "  %12.8f  %12.8f  %12.8f\n", $1, $2, $3} }' >> input
#
################################################################################
##       dynamical generation of the "from_restart_to_xmol.gawk" script       ##
################################################################################
file="from_restart_to_xmol.gawk"
cat > ${file} << EOF
BEGIN { a_zero = 0.529177 }
{
if ( \$1 == "Image:" ) {
  count = -1 ;
  printf "${Ntot}" ;
  printf "     \n" ;
  printf "     \n" ;
  }
else {
  count++;
EOF
#
ref1=0
ref2=0
#
index=0
#
for atom in ${list_of_atoms}; do
  #
  index=$(( ${index} + 1 ))
  #
  ref1=$(( ${ref2} + 1 ))
  ref2=$(( ${ref2} + ${N[${index}]} ))
  #
  echo "  if ( count >= ${ref1} && count <= ${ref2} ) { " >> ${file}
  echo "    printf \"${atom}  \";                       " >> ${file}
  echo "    printf \" %8.4f   \", \$1 * a_zero ;        " >> ${file}
  echo "    printf \" %8.4f   \", \$2 * a_zero ;        " >> ${file}
  echo "    printf \" %8.4f   \", \$3 * a_zero ;        " >> ${file}
  echo "    printf \" \\n\";                            " >> ${file}
  echo "    }                                           " >> ${file}
  #
done
#
echo "  }                                               " >> ${file}
echo "}                                                 " >> ${file}
#
################################################################################
##       dynamical generation of the "from_restart_to_axfs.gawk" script       ##
################################################################################
file="from_restart_to_axsf.gawk"
cat > ${file} << EOF
BEGIN{ 
  a_0  = 0.529177 ;
  printf " ANIMSTEPS %3i \n", ${new_num_of_images} ;
  printf " CRYSTAL       \n" ;
  printf " PRIMVEC       \n" ;
  printf "  %16.10f  %16.10f  %16.10f \n", ${a11}*a_0, ${a12}*a_0, ${a13}*a_0;
  printf "  %16.10f  %16.10f  %16.10f \n", ${a21}*a_0, ${a22}*a_0, ${a23}*a_0;
  printf "  %16.10f  %16.10f  %16.10f \n", ${a31}*a_0, ${a32}*a_0, ${a33}*a_0;
}
{
  if ( \$1 == "Image:" ) {
    count = -1 ;
    printf " PRIMCOORD %3i \n", \$2 ;
    printf "%4i  1 \n", ${Ntot} ;
    }
  else {
    count++;
EOF
#
ref1=0
ref2=0
#
index=0
#
for atom in ${list_of_atoms}; do
  #
  index=$(( ${index} + 1 ))
  #
  ref1=$(( ${ref2} + 1 ))
  ref2=$(( ${ref2} + ${N[${index}]} ))
  #
  echo "  if ( count >= ${ref1} && count <= ${ref2} ) { " >> ${file}
  echo "    printf \"${atom}  \";                       " >> ${file}
  echo "    printf \" %16.10f   \", \$1 * a_0 ;         " >> ${file}
  echo "    printf \" %16.10f   \", \$2 * a_0 ;         " >> ${file}
  echo "    printf \" %16.10f   \", \$3 * a_0 ;         " >> ${file}
  echo "    printf \" %16.10f   \", \$4 / a_0 ;         " >> ${file}
  echo "    printf \" %16.10f   \", \$5 / a_0 ;         " >> ${file}
  echo "    printf \" %16.10f   \", \$6 / a_0 ;         " >> ${file}  
  echo "    printf \" \\n\";                            " >> ${file}
  echo "    }                                           " >> ${file}
  #
done
#
echo "  }                                               " >> ${file}
echo "}                                                 " >> ${file}
#
################################################################################
#
$ROOT_DIR/bin/path_int.x < input
#
$GAWK -f from_restart_to_xmol.gawk ${new_restart_file} > \
${new_restart_file}.xmol
#
$GAWK -f from_restart_to_axsf.gawk ${new_restart_file} > \
${new_restart_file}.axsf
#
rm -f input 
rm -f CELL_PARAMETERS
rm -f from_restart_to_xmol.gawk \
      from_restart_to_axsf.gawk
