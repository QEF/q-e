#!/bin/bash
###############################################
#This script reads the most intense peak position from a post processing output
#and prepares a file for a second stage run giving projections
#OBM 2010
################################################


################################################
#Prelimineries
################################################
 if [ $# -ne 2 ] ; then
  echo "Usage:"
  echo "create_projection.sh <tddfpt first stage input file> <tddfpt_pp first stage output file>"
  exit
 fi
 if [ ! -f $1 ] ; then
  echo "File $1 not found"
  exit
 fi
 if [ ! -f $2 ] ; then
  echo "File $2 not found"
  exit
 fi
 fprefix=`basename $1`
 fprefix=${fprefix%.*}
 fprefix2=`basename $2`
 fprefix2=${fprefix2%.*}
 if [ $fprefix != $fprefix2 ] ; then
  echo "WARNING: a possible mismatch of files!"
 fi

  echo "Processing $1 and $2"
################################################
#Reading the input file
################################################
  prefix=`grep -i prefix $1`
  prefix=${prefix#*=}
  if [ -z $prefix ] ; then 
   echo "Error reading prefix"
   exit
  fi
 
  outdir=`grep -i outdir $1`
  outdir=${outdir#*=}
  if [ -z $outdir ] ; then 
   echo "Error reading outdir"
   exit
  fi
 
  restart_step=`grep -i restart_step $1`
  restart_step=${restart_step#*=}
  if [ -z $restart_step ] ; then 
   echo "Error reading restart_step"
   exit
  fi
 
  itermax=`grep -i itermax $1`
  itermax=${itermax#*=}
  if [ -z $itermax ] ; then 
   echo "Error reading itermax"
   exit
  fi
 
  ipol=`grep -i ipol $1`
  ipol=${ipol#*=}
  if [ -z $ipol ] ; then 
   echo "Error reading ipol"
   exit
  fi
  
  if [ $ipol -eq 4 ] ; then
   npol=3
  else
   npol=$ipol
  fi  


###############################################
#Taking a copy of the Liovillian
###############################################
  outdir_t=${outdir//\"/}
  outdir_t=${outdir_t//,/}
  outdir_t=${outdir_t//\'/}
  prefix_t=${prefix//\"/}
  prefix_t=${prefix_t//,/}
  prefix_t=${prefix_t//\'/}
  for pol in `jot $npol 1`
  do
   cp $outdir_t/$prefix_t.beta_gamma_z.$pol $outdir_t/$prefix_t-orig.beta_gamma_z.$pol
  done
  echo "$npol beta_gamma_z files are copied from $outdir_t/$prefix_t to $outdir_t/$prefix_t-orig"
################################################
#Reading the output file
################################################
  
  peak_pos=`grep "Possible peak at" $2 |gawk 'BEGIN {intensity=-999999;value=0} \$7 > intensity {intensity=\$7;value=\$4} END {print value}'`
  echo "Maximum intensity peak is at $peak_pos Ry"
 
  epsil=`grep -i Broadening $2`
  epsil=${epsil#*=}
  epsil=${epsil%% Ry}
  if [ -z $epsil ] ; then 
   echo "Error reading epsil"
   exit
  fi
 
###############################################
#Writing the file
###############################################



  cat > $fprefix.tddfpt-in-s2 << EOF
&lr_input
    prefix=$prefix
    outdir=$outdir
    restart_step=$restart_step
    restart=.false.
/
&lr_control
    itermax=$itermax
    ipol=$ipol
    charge_response=2
    project=.true.
/
&lr_post
   beta_gamma_z_prefix='$prefix_t-orig'
   omeg=$peak_pos
   epsil=$epsil
   w_T_npol=$npol
   plot_type=3
/

EOF

