
help atom -vartype character -helpfmt txt2html -helptext {
            atomic symbol: atom='H', 'He', 'Be', etc. (REQUIRED)
}

help config -vartype character -helpfmt html  -helptext {
           a string containing the electronic configuration (REQUIRED). 
           Example:<br>   '[Ar] 3d10 4s2 4p2.5' <p>
           In spin-polarised calculations (LSDA), spin-up and spin-down 
           state may appear twice with the respective occupancy: <br>
           3p4 3p2 = 4 up, 2 down. Otherwise, the Hund's rule is assumed.<p>
           In fully relativistic calculations, states with jj=l-1/2 are
           filled first. If a state appears twice, the first one has jj=l-1/2,
           the second one jj=l+1/2 (except S states)<p>
           Negative occupancies can be used to flag unbound states;
           they are not actually used<p>
}

help lsd -vartype integer -helpfmt txt2html -helptext {
           LSDA spin-polarized calculation are not allowed in
           pseudopotential generation or in full relativistic calculation
}

help dft_ -vartype character  -helpfmt txt2html -helptext {
           See module "which_dft" in the Modules/ directory
           for a complete list of acceptable values
}

help title -vartype character -helpfmt txt2html -helptext {
            a string describing the job (optional)
}

help vdw -vartype logicalr -helpfmt txt2html -helptext {
             If .true., the frequency dependent polarizability and van der
             Waals coefficient C6 will be computed in Thomas-Fermi and 
             von Weizsaecker approximation(only for closed-shell ions).
             Gradient-corrected DFT not yet implemented.
}

help prefix -vartype character -helpfmt html -helptext {
            prefix for output file names  containing the orbitals,
            logarithmic derivatives, tests. Optional, default: 'ld1'.<p>
            For the N-th electronic configuration, the following files
            are written:<pre>
            "prefix"N.wfc     all-electron orbitals
            "prefix"Nps.wfc   pseudo-orbitals</pre><br>
            if parameters for logarithmic derivatives are specified:<pre>
            "prefix"Nps.dlog  all-electron logarithmic derivatives
            "prefix"Nps.dlog  pseudo logarithmic derivatives<pre><br>
If pseudopotential testing is performed:<pre>
            "prefix".test    results of transferability test</pre>
            "N" is not present if there is just one electronic configuration
}

help xmin -vartype real -helpfmt txt2html -helptext {
            radial grid parameter (optional, default: -7.0)
}
help dx -vartype real -helpfmt txt2html -helptext {
            radial grid parameter (optional, default: 0.0125)
}
help rmax -vartype real -helpfmt txt2html -helptext {
            outermost grid point (optional, default: 100.0 a.u.)
}

help beta -vartype real -helpfmt txt2html -helptext {
            parameter for potential mixing (optional, default: 0.2)
}
help tr2 -vartype real -helpfmt txt2html -helptext {
            convergence threshold for self-consistency (optional, default: 1e-14)
}

help file_pseudopw -vartype character -helpfmt txt2html -helptext {
            file where the generated PP is written (REQUIRED).

            If the file name ends with "upf" or "UPF",
            or in any case for spin-orbit PP (rel=2),
            the file is written in UPF format

            if the file name ends with 'psp' it is written
            in native CPMD format (currently an experimental feature)

            otherwise it is written in the old "NC" format if it is 
            a norm-conserving PP with one channel per angular momentum
            (single-projector); in the old RRKJ format otherwise
            (these formats are obsolescent)
}

help file_pseudo -vartype character -helpfmt txt2html -helptext {
            file containing the PP to be read (REQUIRED).
            IMPORTANT: for norm-conserving PP with one channel per
            angular momentum, all calculations are done using the
            SEMILOCAL form, not the separable nonlocal form
}

help file_recon  -vartype character -helpfmt txt2html -helptext {
            file containing data needed for PAW reconstruction 
            of all-electron wavefunctions from PP results.
            If you want to use additional states to perform the
            reconstruction, add them at the end of the list
            of all-electron states (optional, default: ' ')
}

help lloc -vartype integer -helpfmt txt2html -helptext {
           local channel: specify either the all-electron potential, 
           or the angular momentum of a channel. In the latter case
           the local channel must be the last in the list of
           states to be pseudized (after the namelist &inputp)
}

help rcloc -vartype real -helpfmt txt2html -helptext {
           matching radius for the local channel: must be specified 
           only if  the all-electron potential is pseudized, 
           otherwise the corresponding value of the matching radius is used
}
help rcore -vartype real -helpfmt txt2html -helptext {
           matching radius for the smoothing of the core charge,
           If not specified, the matching radius is determined 
           by the condition rho_core(rcore) = 2*rho_valence(rcore)
}

help rho0 -vartype real -helpfmt txt2html -helptext {
           charge at the origin: when the Rabe-Rappe-Kaxiras-Joannopoulos
           method with 3 Bessel functions fails, specifying a nonzero charge
           at the origin may allow to override the problem (using 4 Bessel 
           functions) .Typical values are in the order of 0.01-0.02
}

help lpaw -vartype logical -helpfmt txt2html -helptext {
           implemented only for pseudotype=3
}

help zval -vartype real -helpfmt html -helptext {
          Zval is automatically calculated from available data. 
          If one value of Zval is provided in input, it will be
          checked versus the calculated value. The only case in
          which you need to explicitly provide the value of Zval
          is for nointeger Zval (i.e. half core-hole pseudopotentials)
}

help PP_nwfs -vartype integer -helpfmt txt2html -helptext {
           number of wavefunctions to be pseudized
}

help PP_wfs -helpfmt html -helptext {
<I>Label:</I> wavefunction label, the same as in the reference
              all-electron configuration (example: 4s, 4p, ...)<p>
<I>N:</I>  principal quantum number (referred to the PSEUDOPOTENTIAL case!
           1 for lowest s, 2 for lowest p, and so on)<br>
<I>L:</I>  angular momentum quantum number<br>
<I>Occupation:</I>  the same as in the reference all-electron configuration<br>
<I>Energy:</I>  energy used to pseudize the corresponding state. If set to 0
              or not set: use the one-electron energy of the all-electron 
              state. Specify a nonzero energy for unbound states.<br>
<I>Rcut:</I>    matching radius for norm conserving pseudization<br>
<I>Rcutus:</I>  matching radius for ultrasoft pseudization (if any)<p>
      NB: if a state with a given L is used as local channel, that state
      must be the last in the list of states to be pseudized
}

help pseudotype -vartype integer -helpfmt txt2html -helptext {
            IMPORTANT: for norm-conserving PP with one channel per
            angular momentum, all calculations are done using the
            SEMILOCAL form, not the separable nonlocal form
}

help configts -vartype character -helpfmt html  -helptext {
           testing electronic configuration(s). Only the VALENCE part,
           i.e. states that have been pseudized, must be specified.
           The CORE states are the same as specified in the all-electron
           configuration. Example: '4s2 4p2.5'. <p>
           For PP generation you do not need to specify a test configuration,
           UNLESS:<p>
           1) you want to use a different configuration for unscreening 
              with respect to the one used to generate the PP. This is
              useful for PP with semicore states: use semicore states ONLY
              to produce the PP, use semicore AND valence states (if occupied)
              to make the unscreening<p>
           2) you want to specify additional states for PAW reconstruction
              of all-electron orbitals from pseudo-orbitals
}



help ecutmax -vartype real -helpfmt html  -helptext {
              ecutmin, ecutmax, decut: parameters (Ry) used for test with a
              basis set of spherical Bessel functions j<sub>l</sub>(qr) . The
              hamiltonian at fixed scf potential is diagonalized for various
              values of ecut:
              ecutmin, ecutmin+decut, ecutmin+2*decut ... up to ecutmax.
              This yields an indication of convergence with the corresponding
              plane-wave cutoff in solids, and shows in an unambigouous way if
              there are "ghost" states <p>
              Default: decut=5.0 Ry, ecutmin=ecutmax=0Ry; specify ecutmin and 
              ecutmax if you want to perform this test
}

help rm -vartype real -helpfmt html  -helptext {
                    Radius (a.u.) of the box used with spherical Bessel
                    functions (default: 30 a.u.)
}
