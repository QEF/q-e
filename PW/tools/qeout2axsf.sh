#!/bin/bash --noprofile
################################################################################
##                    Copyright (C) 2006 Carlo Sbraccia.                      ##
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
if [ "$1" == ""  ]; then
  echo "input file missing"; exit
fi
#
filename=$1
#
alat=$( grep 'alat\|a_0' ${filename} | \
        head -n 1 | awk -F = '{print $2}' | awk '{print $1}' )
#
if grep "PWSCF" ${filename} &> /dev/null ; then code="PW"; fi
if grep "CP:"   ${filename} &> /dev/null ; then code="CP"; fi
#
string=$( cat ${filename} | awk -v code=${code} -v alat=${alat} ' \
BEGIN{ \
  iter = 0 ; \
  done = 0 ; \
  b2a  = 0.529177 ; \
  if ( code == "PW" ) { runt = b2a*alat } \
  if ( code == "CP" ) { runt = b2a } \
} \
{ \
  if ( $1 == "ATOMIC_POSITIONS" ){ \
    if ( match( toupper( $0 ), "ANGSTROM" ) ){ runt = 1.0 } ; \
    if ( match( toupper( $0 ), "BOHR"     ) ){ runt = b2a } ; \
    iter++ ; \
    if ( done == 0 ){ \
      nat = 0 ; \
      getline ; \
      while ( NF == 4 || NF == 7 ){ ++nat ; getline } \
      done = 1 ; \
    } \
  } \
} \
END{ printf "%d %d %10.8f", nat, iter, runt } ' )
#
nat=$(   echo ${string} | awk '{print $1}' )
niter=$( echo ${string} | awk '{print $2}' )
runt=$(  echo ${string} | awk '{print $3}' )
#
cat ${filename} | awk -v runt=${runt} -v nat=${nat} \
                       -v alat=${alat} -v niter=${niter} ' \
BEGIN{ \
  iter = 0 ; \
  b2a  = 0.529177 ; \
  printf "ANIMSTEPS %5d\n", niter ; \
  printf "CRYSTAL\n" ; \
  printf "PRIMVEC\n" ; \
  printf "%14.10f %14.10f %14.10f\n", alat*b2a, 0.0, 0.0 ; \
  printf "%14.10f %14.10f %14.10f\n", 0.0, alat*b2a, 0.0 ; \
  printf "%14.10f %14.10f %14.10f\n", 0.0, 0.0, alat*b2a ; \
} \
{ \
  if ( $1 == "ATOMIC_POSITIONS" ){ \
    printf "PRIMCOORD %5d\n", ++iter ; \
    printf "%5d 1\n", nat ; \
    for ( i = 1; i <= nat; ++i ){ \
      getline ; \
      printf "%3s   %14.9f%14.9f%14.9f\n", $1, $2*runt, $3*runt, $4*runt ; \
    } \
  } \
} ' > ${filename}.axsf
#
printf "\nalat = %12.8f Bohr\n" ${alat}
printf "\npositions in alat coordinates :\n\n"
#
tail -n ${nat} ${filename}.axsf | awk -v alat=${alat} ' \
BEGIN{ \
  angstrom2alat = 1.0 / 0.529177 / alat ; \
} \
{ \
  printf "%3s   %14.9f%14.9f%14.9f\n", $1, $2*angstrom2alat,  \
                                           $3*angstrom2alat,  \
                                           $4*angstrom2alat ; \
} '
