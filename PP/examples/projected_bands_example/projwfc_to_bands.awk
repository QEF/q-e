# Analyzes the output of projwfc.x to extract the bands Epsilon(k)
# and possibly the weights onto specific atomic states
#
# execution:
# awk -v ef=[EF] -v firststate=[FIRSTSTATE] -v laststate=[LASTSTATE] -f projwfc_to_bands.awk [PROJWFCOUTPUT]
#
# variables to be initialized when calling (all are optional)
# ef : Fermi energy                         (if unspecified/0: don't translate)
# firststate : first atomic wfc to consider (if unspecified/0: write bands only)
# laststate  :  last atomic wfc to consider (if unspecified/0: write bands only)
#
# Guido Fratesi 2017-01-23
# Tested on output of QuantumESPRESSO 6.0

# set the atomic wavefunctions to be included
BEGIN {
    natwfc=0;
    firststate=strtonum(firststate);
    laststate=strtonum(laststate);
    if ((firststate)&&(laststate)) {
	for (i=firststate;i<=laststate;i++) {
	    latwfc[natwfc]=sprintf("[#%4s]", i);
	    natwfc++;
	}
    }
}

# as an alternative, the atomic wavefunctions can be specified manually here:
#BEGIN {
#    natwfc=0;        #1234567#
#    latwfc[natwfc++]="[#   1]";
#    latwfc[natwfc++]="[#   2]";
#    latwfc[natwfc++]="[#   3]";
#    latwfc[natwfc++]="[#   4]";
#}

# initialize length of the path in reciprocal space
BEGIN {
    klen=0;
    ik=0;
}

# specify which states have been included
/ state #/ {
    iii=index($0,":");
    jstate=strtonum(substr($0,iii-4,4));
    for (iatwfc=0;iatwfc<natwfc;iatwfc++) {
	istate=strtonum(substr(latwfc[iatwfc],3,4));
	if (istate==jstate) printf "# %s\n", $0;
    }
}

# write file header
/ k = / && (ik==0) {
    printf "# energies were translated by aligning to the reference value EF= %11.5f\n", ef;
    printf "%3s %11s", "#iK", "K-length";
    printf " %11s", "E-EF(eV)";
    if (natwfc) {
	printf " %7s", "[TOTAL]";
	for (iatwfc=0;iatwfc<natwfc;iatwfc++) {
	    dum=latwfc[iatwfc];
	    gsub(" ","_",dum);
	    printf " %7s", dum;
	}
    }
    printf "\n";
}

# measure length of the path in reciprocal space
/ k = / {
    kx=$3; ky=$4; kz=$5;
    if (ik) {
	dk=sqrt((kx-kxo)^2+(ky-kyo)^2+(kz-kzo)^2);
	klen+=dk;
    }
    ik++;
    kxo=kx; kyo=ky; kzo=kz;
}

# new wavefunction: get energy and set weights to zero
/==== e/ {
    e=$(NF-2);
    if (natwfc) {
	for (iatwfc=0;iatwfc<natwfc;iatwfc++) wfcweight[iatwfc]=0;
	totweight=0;
    }
}

# scan for weights on selected atomic wavefunctions
{
    for (iatwfc=0;iatwfc<natwfc;iatwfc++) {
	iii=index($0,latwfc[iatwfc]);
	if (iii) {
	    wfcweight[iatwfc]=substr($0,iii-6,5);
	    totweight+=strtonum(wfcweight[iatwfc]);
#DEBUG#	    print $0;
#DEBUG#	    print wfcweight[iatwfc];
	}
    }
}

# end of information about the current wavefunction: printout
/\|psi\|/ {
    printf "%3i %11.5f", ik, klen;
    printf " %11.5f", e-ef;
    if (natwfc) {
	printf " %7.3f", totweight;
	for (iatwfc=0;iatwfc<natwfc;iatwfc++) printf " %7.3f", wfcweight[iatwfc];
    }
    printf "\n";
}
