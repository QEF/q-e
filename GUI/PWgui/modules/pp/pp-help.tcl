
help  prefix -vartype character -helpfmt txt2html -helptext {
     prefix of files saved by program pwscf
}

help  outdir -vartype character -helpfmt txt2html -helptext {
     temporary directory where pwscf files resides
}
  
help  filplot  -vartype character -helpfmt txt2html -helptext {
     punch file, contains the quantity selected by plot_num
}

help  plot_num -vartype integer -helpfmt txt2html -helptext {
   Selects what is saved in filplot:
                 0 = charge
                 1 = total potential V_bare+V_H + V_xc
                 2 = local ionic potential
                 3 = local density of states at e_fermi
                 4 = local density of electronic entropy
                 5 = STM images
                 6 = spin polarization (rho(up)-rho(down))
                 7 = |psi|^2
                 8 = electron localization function (ELF)
                 9 = planar average of all |psi|^2
                10 = integrated local density of states from
                     emin to emax (emin, emax in eV)
                     if emax is not specified, emax=E_fermi
                11 = the V_bare + V_H potential
                12 = the electric field potential
                13 = the noncolinear magnetization.
}
  
help  spin_component -vartype integer -helpfmt txt2html -helptext {
   For charge and potential:
                 0 = total charge/potential 
                 1 = spin up charge/potential,
                 2 = spin down charge/potential.

   For noncolinear magnetization:
                 0 = absolute value (default value)
                 1 = x component of the magnetization
                 2 = y component of the magnetization
                 3 = z component of the magnetization

<p> ( default = 0 )
}
  
help  sample_bias  -vartype real -helpfmt txt2html -helptext {
     the bias of the sample (Ryd) in STM images
}

help  stm_wfc_matching  -vartype logical -helpfmt txt2html -helptext {
     if .t. match the wavefunctions to an exponentially vanishing function
     if .t. specify also (in celldm(1) units):
                  z    height of matching
                 dz    distance of next stm image calculation
<p> ( default = .false. )
}

help z -vartype real -helpfmt txt2html -helptext {
     height of matching
}
help dz -vartype real -helpfmt txt2html -helptext {
     distance of next stm image calculation
}
  
help  kpoint  -vartype integer -helpfmt txt2html -helptext {
     for which k-point to compute property
}

help  kband   -vartype integer -helpfmt txt2html -helptext {
     for which band to compute property
}

help  lsign   -vartype logical -helpfmt txt2html -helptext {
     if .true. and k point is Gamma, save |psi|^2 sign(psi)
}

set _e  {
    emin        lower energy boundary (in eV)
    emax        upper energy boundary (in eV), i.e. compute ILDOS from emin to emax
}
foreach var {emin emax} {
    help $var -vartype real -helpfmt txt2html -helptext $_e
}