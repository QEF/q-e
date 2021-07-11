plot "JDOS_Gaus_equil.dat" u 1:(sqrt($2/19562791)) w l,\
     "JDOS_Gaus_333_300K.dat" u 1:(sqrt($2/25088078)) w l,\
     "JDOS_Gaus_333_300K_nosym.dat" u 1:(sqrt($2/25044853)) w l,\
     "JDOS_Gaus_333_300K_qA.dat" u 1:(sqrt($2/25220884)) w l
