#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwo2xsf.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#------------------------------------------------------------------------
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# pwo2xsf: PW-output--to--XSF conversion; created for the purpose 
#          of PW spring college !!!
#
# Usage: pwo2xsf [options] pw-output-file
#
#------------------------------------------------------------------------

#############
BOHR=0.529177
#############

cat > pwo2xsfUsage.$$ <<EOF

 Usage: pwo2xsf.sh [options] [pw-output-file]

 Options are:

               pwo2xsf.sh --initcoor|-ic [pw-output-file]
                             Will extract the initial (i.e. input) 
                             ionic coordinates.

               pwo2xsf.sh --latestcoor|-lc [pw-output-file]
                             Will extract latest estimation of ionic 
                             coordinates from pw-output file. The coordinates 
                             can be either taken from "Search of equilibrium 
                             positions" record or from "Final estimate of 
                             positions" record.

               pwo2xsf.sh --optcoor|-oc [-A|-B|-F] [pw-output-file]
			     Similar to "--latestcoor", but extract just the 
			     optimized coordinates from "Final estimate of 
                             positions" record.

	       pwo2xsf.sh --animxsf|-a [pw-output-file1] ...
			     Similar to "--latestcoor", but extract the
			     coordinates from all ionic steps and make
			     a AXSF file for animation.

               pwo2xsf.sh -r option [pw-output-file]
			     one must specify i.e. ityp->nat conversion, 
			     and the corresponding data are written to file 
			     nuclei.charges. The -r flag deletes this file.
                             Here "option" is one of the above options
			     (i.e -lc|-oc|-a).
EOF

pwoExit() {
    rm -f pwo2xsfUsage.$$
    exit 1
}
pwoUsage() {
    if [ $1 ]; then
	echo "
Usage: $2
" 1>&2
	pwoExit
    fi
}
pwError() {
    echo "
 ========================================================================
    $1
 ========================================================================
" 1>&2
    if [ "$2" -ge 0 ]; then
	exit $2
    fi
}

# --------------------------------------------------------------------------
# FUNCTION:   pwNucleiCharges --
#
# Purpose:    ityp->nat conversion data 
#
# Usage:      pwNucleiCharges pw_input|pw_output outfile
#
# Side efect: creates nuclei.charges file
# --------------------------------------------------------------------------

pwNucleiCharges() {
    #
    # if file nuclei.charges does not exists prompt for ityp->nat conversion !!
    #

    if [ \( "$1" = "" \) -o \( "$2" = "" \) ]; then
	pwError "Usage:  pwNucleiCharges  pw_input|pw_output  outfile" 1
    fi
    
    # do we have PW-INPUT or PW-OUTPUT file ???
    
    if [ "`cat $1 | egrep '&input|&system'`" != "" ]; then
	# it is PW-INPUT
	ntyp=`cat "$1" | awk '{gsub(",","\n"); print}' | grep ntyp \
		| awk '{split($0,a,"=|,"); print a[2];}'`
    else
	# PW-OUTPUT
	ntyp=`cat "$1" | grep 'number of atomic types' | \
	    head -1 | awk '{print $NF}'`
	#echo 'NTYP=$ntyp'
	if [ "$ntyp" = "" ]; then
	    # some older PWSCF versions didn't have "number of atomic
	    # types" printout -> user will have to make nuclei.charges
	    # file by himself/herself !!!
	    pwError "This is either non PW-output file or is a PW-output file 
    produced with some old PWSCF version" -1
	    echo -n "How many ityp->nat replacements ? " 1>&2
	    read ntyp
	fi
    fi		     

    if [ ! -f nuclei.charges ]; then
	echo -n "Please enter $ntyp ityp->nat replacements !!! " 1>&2
	echo $ntyp > nuclei.charges
	i=0
	while [ $i -lt "$ntyp" ]
	do
	    i=`echo "$i + 1"|bc`
	    echo "" 1>&2
	    echo "Replacement #${i}: ityp->nat" 1>&2
	    echo -n "ityp[$i]=$i; nat[$i]=" 1>&2; read nat
	    echo "$i $nat" >> nuclei.charges
	done
    fi
	
    cat nuclei.charges > "$2"
}

# ---------------------------------------------------------------------------
# read PW-output file and print the XSF file according to specified flags
pwoPrintCoor() {
    #set -x
    pwoUsage "$# -lt 1" \
	"$0 --latestcoor|-lc [pw-output-file]   or   $0 --optcoor|-oc [pw-output-file]"
   
    case $1 in
	--latestcoor|-lc) type=lc; shift;;
	--optcoor|-oc)    type=oc; shift;;
    esac

    if [ $# -eq 0 ]; then
	pwNucleiCharges -  pw.$$
	cat - >> pw.$$
    else
	pwNucleiCharges $1  pw.$$
	for i in `seq 1 $#`
	do
	    cat $1 >> pw.$$
    	    shift	    
	done
    fi
    inp=pw.$$

    cat "$inp" | awk -v bohr=$BOHR -v t=$type -- '
function isInt(__num) {
   if ( __num ~ /^[0-9]+$/ ) return 1;
   else return 0;
} 
function PrintItyp(__ityp) {
   if ( __ityp ~ /^[0-9]+$/ ) return atn[ __ityp ];
   else return __ityp;
}
BEGIN { 
    line=0; optc=0; 
    for (i=0; i<100; i++) atn[i]=i;
    getline;
    ntyp=$1;
    for (i=0; i<ntyp; i++) {
	getline;
	atn[$1] = $2;
    }
}

/bravais-lattice index/       { ibrav=$NF; }
/celldm\(1\)=/                { a0=$2*bohr; }
/crystal axes:/               {
    for (i=0; i<3; i++) {
	getline;
	for (j=4; j<7; j++) {
	    vec[i,j-4] = $j * a0;
	}
    }
    printf " DIM-GROUP\n 3 1\n PRIMVEC\n";
    printf "%15.10f %15.10f %15.10f\n", vec[0,0], vec[0,1], vec[0,2];
    printf "%15.10f %15.10f %15.10f\n", vec[1,0], vec[1,1], vec[1,2];
    printf "%15.10f %15.10f %15.10f\n", vec[2,0], vec[2,1], vec[2,2];
}
/number of atoms/             { 
    nat=$NF; 
}
/Final estimate of positions/ { 
    optc=1; line=1; 	    
    printf " PRIMCOORD\n %d 1\n", nat;
    next; 
}
/Cartesian coordinates/ {
    if (optc==1 && line==1) {
	next;
    }
}
/Search of equilibrium positions/ { if ( t == "lc" ) {
    getline; getline; 
    # maybe current line is: "Cartesian coordinates"
    if ( NF != 4 ) getline;
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=a0*$1; y[i]=a0*$2; z[i]=a0*$3; ityp[i]=$4;
        } else { 
           x[i]=a0*$2; y[i]=a0*$3; z[i]=a0*$4; ityp[i]=$1;
        }
	getline;
    }
  }
}
/Entering Dynamics;/ {
    getline; getline; getline; nstep++;
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=a0*$1; y[i]=a0*$2; z[i]=a0*$3; ityp[i]=$4;
        } else { 
           x[i]=a0*$2; y[i]=a0*$3; z[i]=a0*$4; ityp[i]=$1;
        }
	getline;
    }
}

/a*/ { if (line == 1) {
	    if (NF != 4) exit 0;	    	    
            if (isInt($4)) printf "% 3d   % 15.10f  % 15.10f  % 15.10f\n",
		                  atn[$4], a0*$1, a0*$2, a0*$3; 
            else printf "% 3s   % 15.10f  % 15.10f  % 15.10f\n",
		        $1, a0*$2, a0*$3, a0*$4;
       } 
}
END { if ( t == "lc" && line == 0) {
    printf " PRIMCOORD\n %d 1\n", nat;
    for(i=1; i<=nat; i++)
       printf "% 3s   % 15.10f  % 15.10f  % 15.10f\n", 
	    PrintItyp(ityp[i]), x[i], y[i], z[i];
       }
    }' | tee pwo2xsf.xsf_out
}

# ---------------------------------------------------------------------------
# read PW-output file and print the animated-XSF file
pwoAnimXSF() {
    #set -x
    pwoUsage "$# -lt 1" "$0 --animxsf|-a [pw-output-file1] ..."
    
    only_init=0
    case $1 in
	--inicoor|-ic) only_init=1;;
    esac

    if [ $# -eq 1 ]; then
	pwNucleiCharges -  pw.$$    
	cat - >> pw.$$
    else
	pwNucleiCharges $2  pw.$$
	# if "pw-output-file" is glob expresion -> must merge all outputs
	for i in `seq 2 $#`
	do
	    cat $2 >> pw.$$
	    shift	    
	done
    fi
    inp=pw.$$

    nstep=`egrep "Final estimate of positions|Search of equilibrium positions|Entering Dynamics;" $inp | wc | awk '{print $1}'`
    # add also initial coordinates
    nstep=`expr $nstep + 1`

    cat "$inp" | awk -v bohr=$BOHR -v astep=$nstep -v onlyinit=$only_init -- '
function isInt(__num) {
   if ( __num ~ /^[0-9]+$/ ) return 1;
   else return 0;
}
function PrintItyp(__ityp) {
   if ( __ityp ~ /^[0-9]+$/ ) return atn[ __ityp ];
   else return __ityp;
}

function PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz) {
    if (onlyinit) {
       print " PRIMCOORD";
    } else {
       print " PRIMCOORD", nstep;
    }
    print nat, 1;
    for(i=1; i<=nat; i++) {
	printf "% 3s ", PrintItyp(ityp[i]); 
        printf "% 15.10f  % 15.10f  % 15.10f   % 15.10f  % 15.10f  % 15.10f\n", x[i], y[i], z[i], fx[i], fy[i], fz[i];
    }
}
function GetInitCoor(nstep, nat, a0, ityp, x, y, z) {
    for(i=1; i<=nat; i++) {
        ityp[i]=$2;
	split($0,rec,"("); split(rec[3],coor," ");
	x[i]= a0*coor[1]; y[i]=a0*coor[2]; z[i]=a0*coor[3];
	getline;
    }
}
function GetPrimCoor(nstep, nat, a0, ityp, x, y, z) {
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=a0*$1; y[i]=a0*$2; z[i]=a0*$3; ityp[i]=$4;
        } else { 
           x[i]=a0*$2; y[i]=a0*$3; z[i]=a0*$4; ityp[i]=$1;
        }
	getline;
    }
}
function GetForces(nstep, nat, ityp, fx, fy, fz) {
    for(i=1; i<=nat; i++) {
	if (nstep==1) ityp[i]=$4;
	fx[i]=$7; fy[i]=$8; fz[i]=$9;
	getline;
    }
}

BEGIN { 
    nstep=1; print_celldm=1; 
    for (i=0; i<100; i++) atn[i]=i;
    getline;
    ntyp=$1;
    for (i=0; i<ntyp; i++) {
	getline;
	atn[$1] = $2;
    }
}

/bravais-lattice index/       { ibrav=$NF; }

/celldm\(1\)=/                { a0=$2*bohr; }

/crystal axes:/               {
    if (print_celldm) {
	for (i=0; i<3; i++) {
	    getline;
	    for (j=4; j<7; j++) {
		vec[i,j-4] = $j * a0;
	    }
	}
        #---
        # this is now donw posteriori (see below)
        #if (!onlyinit) printf " ANIMSTEPS %d\n", astep
        #---
	printf " DIM-GROUP\n 3 1\n PRIMVEC\n";
	printf "%15.10f %15.10f %15.10f\n", vec[0,0], vec[0,1], vec[0,2];
	printf "%15.10f %15.10f %15.10f\n", vec[1,0], vec[1,1], vec[1,2];
	printf "%15.10f %15.10f %15.10f\n", vec[2,0], vec[2,1], vec[2,2];
	print_celldm=0;
    }
}

/number of atoms/  { nat=$NF; }

/Cartesian axes|Cartesian axes/  { 
    getline; getline; getline;
    if (nstep == 1) GetInitCoor(nstep, nat, a0, ityp, x, y, z);
    if (onlyinit == 1) {
       PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz);
       exit;
    }
}
/Final estimate of positions/     { 
    getline; nstep++;
    # maybe current line is: "Cartesian coordinates"
    if ( NF != 4 ) getline;
    GetPrimCoor(nstep, nat, a0, ityp, x, y, z);
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz); 
}

/Search of equilibrium positions/ { 
    getline; getline; nstep++;
    # maybe current line is: "Cartesian coordinates"
    if ( NF != 4 ) getline;
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    GetPrimCoor(nstep, nat, a0, ityp, x, y, z);
}
/Entering Dynamics;/ {
    getline; getline; getline; nstep++;
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    GetPrimCoor(nstep, nat, a0, ityp, x, y, z);
}
/Forces acting on atoms/          { 
    getline; getline; 
    GetForces(nstep, nat, ityp, fx, fy, fz); 
    PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz); 
}' > pwo2xsf.xsf_out

    if [ $only_init -eq 0 ]; then
    # Assign the number of ANIMSTEPS here. The reason is that the
    # output file (queue runs) is the result of several job runs, then
    # some of them might be terminated on the "wrong" place, and the
    # initial ANIMSTEPS might be wrong. The most secure way is to extract the 
    # sequential digit from the last "PRIMCOORD id" record.
	nsteps=`grep PRIMCOORD pwo2xsf.xsf_out | tail -1 | awk '{print $2}'`    
	echo " ANIMSTEPS $nsteps" > pwo2xsf.xsf_out1
	cp pwo2xsf.xsf_out pwo2xsf.xsf_out2
	cat pwo2xsf.xsf_out1 pwo2xsf.xsf_out2 | tee pwo2xsf.xsf_out
    else
	cat pwo2xsf.xsf_out
    fi
}


#######################################################################
####                              MAIN                              ###
#######################################################################
if [ $# -eq 0 ]; then
    cat pwo2xsfUsage.$$
    pwoExit
fi

r=0
if [ "$1" = "-r" ]; then
    r=1
    shift
fi

case $1 in
    --inicoor|-ic)    pwoAnimXSF $@;;
    --latestcoor|-lc) pwoPrintCoor $@;;
    --optcoor|-oc)    pwoPrintCoor $@;;
    --animxsf|-a)     pwoAnimXSF $@;;    
    *) cat pwo2xsfUsage.$$; pwoExit;;
esac

rm -f  pwo2xsfUsage.$$ pw.$$
if [ $r -eq 1 ]; then
    rm nuclei.charges
fi

exit 0
