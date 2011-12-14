#!/bin/bash --noprofile
################################################################################
##                    Copyright (C) 2004 Guido Fratesi.                       ##
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
#
# script modified by Carlo Sbraccia
#
if [ -z $1 ] || [ ! -e $1 ]; then
   #
   echo "which file? usage:";
   echo "md_analyzer <file> [step_low step_high]"; exit;
   #
fi
#
user=$( whoami )
#
file=$1
fileout=/tmp/$file.dat
#
if [ -n $3 ]; then lra=$2; hra=$3; else lra="1"; hra="*"; fi
#
# while [ 1 -le 2 ]; do
   #
   awk 'BEGIN{ start = 1 } \
        ( $1 == "!" )    { if ( start == 1 ) { epot0 = $5 } \
                           epot = $5 - epot0 } \
        /kinetic energy/ { if ( start == 1 ) { ekin0 = $5 } \
                           ekin = $5 - ekin0; getline ; \
                           temp = $3;         getline ; \
                           if ( start == 1 ) { etot0 = $6 } \
                           etot = $6 - etot0 ; \
                           if ( start == 1 ) { start = 0 } \
        printf "%3i  %16.10f  %16.10f  %16.10f  %8.3f\n", \
                           ++it, etot, ekin, epot, temp }' $file > $fileout
   #
   kill $( ps -u $user | grep gnuplot_x11 | awk '{print $1}' ) &> /dev/null
   #
   #for term in X11 post; do
      #
cat  << EOF | gnuplot -persist
#set term $term
#set out '/tmp/$file.ps'
set da s l
set xra [$lra:$hra]
set lmargin 10
set origin 0,0; set size 1,1
set multiplot
set origin 0,0; set size 1,0.3
plot '$fileout' u 1:2 t "Etot"
unset xlabel
set origin 0,0.3; set size 1,0.4
plot '$fileout' u 1:2 t "Etot", '$fileout' u 1:3 t "Ekin", \
     '$fileout' u 1:4 t "Epot"
set origin 0,0.7; set size 1,0.3
set title 'MD - $file'
plot '$fileout' u 1:5 t "T"
set nomultiplot
#reread
EOF
      #
   #done
   #
   # sleep 5
   #
# done
