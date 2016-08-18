
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help title -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>title</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Status: </em> OPTIONAL
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> A string describing the job.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help zed -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zed</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>See: </em> atom
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The nuclear charge (1 &lt; zed &lt; 100).

<b>IMPORTANT:</b>
Specify either "zed" OR "atom", not both!
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help atom -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>atom</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>See: </em> zed
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Atomic symbol: atom='H', 'He', 'Be', etc.

<b>IMPORTANT:</b>
Specify either "atom" OR "zed", not both!
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help xmin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>xmin</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
-7.0 if "iswitch"&gt;1 or "rel"=0,
-8.0 otherwise
            </li>
<br><li> <em>See: </em> dx
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Radial grid parameter.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dx -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dx</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
0.0125 if "iswitch"&gt;1,
0.008 otherwise
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Radial grid parameter.

The radial grid is: r(i+1) = exp(xmin+i*dx)/zed  a.u.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rmax -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rmax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 100.0 a.u.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Outermost grid point.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help beta -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>beta</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> parameter for potential mixing
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tr2 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tr2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1e-14
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> convergence threshold for scf
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help iswitch -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>iswitch</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
1    all-electron calculation
2    PP test calculation
3    PP generation
4    LDA-1/2 correction, needs a previously generated PP file
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nld -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nld</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> the number of logarithmic derivatives to be calculated
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rlderiv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rlderiv</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> radius (a.u.) at which logarithmic derivatives are calculated
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {eminld emaxld} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>eminld, emaxld</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Energy range (min, max energy, in Ry) at which
logarithmic derivatives are calculated.
            </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help deld -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>deld</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Delta e (Ry) of energy for logarithmic derivatives.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rpwe -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rpwe</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> rlderiv
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> radius (a.u.) at which partial wave expansions are calculated
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rel -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rel</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em>
0 for Z &lt;= 18;
1 for Z &gt;  18
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 ... non relativistic calculation
1 ... scalar relativistic calculation
2 ... full relativistic calculation with spin-orbit
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsmall -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsmall</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. writes on files the small component
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help max_out_wfc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>max_out_wfc</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 7
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of atomic wavefunctions written in the output
file.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help noscf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>noscf</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. the charge is not computed. The occupations are
not used and the eigenvalues and eigenfunctions are those
of a hydrogen-like atom.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsd -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsd</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 ... non spin polarized calculation
1 ... spin-polarized calculation

<b>BEWARE:</b>
not allowed if iswitch=3 (PP generation) or with full
relativistic calculation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dft -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dft</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Exchange-correlation functional.

Examples:
<b>'PZ'</b>    Perdew and Zunger formula for LDA
<b>'PW91'</b>  Perdew and Wang GGA
<b>'BP'</b>    Becke and Perdew GGA
<b>'PBE'</b>   Perdew, Becke and Ernzerhof GGA
<b>'BLYP'</b>  ...

For the complete list, see module "functionals" in ../Modules/
The default is <b>'PZ'</b> for all-electron calculations,
it is read from the PP file in a PP calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help latt -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>latt</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 ... no Latter correction
1 ... apply Latter correction
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help isic -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>isic</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Status: </em>
only for all-electron calculation
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 ... no Self-interaction correction
1 ... apply Self-interaction correction
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rytoev_fact -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rytoev_fact</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> as specified in file Modules/constants.f90
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Factor used to convert Ry into eV.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help cau_fact -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>cau_fact</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> as specified in file Modules/constants.f90
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Speed of light in a.u..

(Be careful the default value is always used in the
 relativistic exchange.)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help vdw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>vdw</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Status: </em>
Gradient-corrected DFT not yet implemented.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true., the frequency dependent polarizability and van der
Waals coefficient C6 will be computed in Thomas-Fermi and
von Weizsaecker approximation(only for closed-shell ions).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help prefix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>prefix</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'ld1'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Prefix for file names - only for output file names
containing the orbitals, logarithmic derivatives, tests
See below for file names and the content of the file.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help verbosity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>verbosity</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'low'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
<b>'low'</b> or <b>'high'</b>

if <b>'high'</b> with iswitch=2,3 prints separately core and
valence contributions to the energies. Print the
frozen-core energy.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_charge -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_charge</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Name of the file where the code writes the all-electron
total charge. No charge is written if file_charge=' '.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help config -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>config</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
A string with the electronic configuration.

Example:
  '[Ar] 3d10 4s2 4p2.5'

* If "lsd"=1, spin-up and spin-down state may appear twice
  with the respective occupancy: 3p4 3p2 = 4 up,
  2 down. Otherwise, the Hund's rule is assumed.

* If "rel"=2, states with jj=l-1/2 are filled first.
  If a state appears twice, the first one has jj=l-1/2,
  the second one jj=l+1/2 (except S states)
  (Use rel_dist if you want to average the electrons
  over all available states.)

* If config='default' the code uses "zed" to set the ground
  state electronic configuration for the atom.

Negative occupancies are used to flag unbound states;
they are not actually used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help relpert -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>relpert</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the relativistic corrections to the non-relativistic
Kohn-Sham energy levels ("rel"=0 .and. "lsd"=0) are computed using
first-order perturbation theory in all-electron calculations.
The corrections consist of the following terms:
   E_vel: velocity (p^4) correction
   E_Dar: Darwin term
   E_S-O: spin-orbit coupling
The spin-orbit term vanishes for s-electron states and gives
rise to a splitting of (2*l+1)*E_S-O for the other states.
The separate contributions are printed only if verbosity='high'.

Formulas and notation are based on the Herman-Skillman book:
F. Herman and S. Skillman, "Atomic Structure Calculations",
Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1963
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rel_dist -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rel_dist</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'energy'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
<b>'energy'</b> or <b>'average'</b>

* if <b>'energy'</b> the relativistic l-1/2 states are filled first.

* if <b>'average'</b> the electrons are uniformly distributed
  among all the states with the given l.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help write_coulomb -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>write_coulomb</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true., a fake pseudo-potential file with name X.UPF,
where X is the atomic symbol, is written. It contains
the radial grid and the wavefunctions as specified in input,
plus the info needed to build the Coulomb potential
for an all-electron calculation - for testing only.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nwf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nwf</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of wavefunctions
                     </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help AE_wfs -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>nl</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> wavefunction label (e.g. 1s, 2s, etc.)
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>n</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> principal quantum number
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>l</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> angular quantum number
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>oc</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> occupation number
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>isw</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> the spin index (1-2) used only in the lsda case
                        </pre></blockquote>
</ul>   
    

    <ul>
<li> <em>Variable: </em><big><b>jj</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The total angular momentum (0.0 is allowed for complete
shells: the codes fills 2l states with jj=l-1/2,
2l+2 with jj=l+1/2).
                        </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help zval -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zval</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> (calculated)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Valence charge.

zval is automatically calculated from available data.
If the value of zval is provided in input, it will be
checked versus the calculated value. The only case in
which you need to explicitly provide the value of zval
for noninteger zval (i.e. half core-hole pseudo-potentials).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help pseudotype -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>pseudotype</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
1 ... norm-conserving, single-projector PP
      IMPORTANT: if pseudotype=1 all calculations are done
      using the SEMILOCAL form, not the separable nonlocal form

2 ... norm-conserving PP in separable form (obsolescent)
      All calculations are done using SEPARABLE non-local form
      IMPORTANT: multiple projectors allowed but not properly
      implemented, use only if you know what you are doing

3 ... ultrasoft PP or PAW
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_pseudopw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_pseudopw</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
File where the generated PP is written.

* if the file name ends with "upf" or "UPF",
or in any case for spin-orbit PP (rel=2),
the file is written in UPF format;

* if the file name ends with 'psp' it is
written in native CPMD format (this is currently
an experimental feature); otherwise it is written
in the old "NC" format if pseudotype=1, or
in the old RRKJ format if pseudotype=2 or 3
(no default, must be specified).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_recon -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_recon</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
File containing data needed for GIPAW reconstruction
of all-electron wavefunctions from PP results.
If you want to use additional states to perform the
reconstruction, add them at the end of the list
of all-electron states.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lloc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lloc</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> -1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Angular momentum of the local channel.

* lloc=-1 or lloc=-2 pseudizes the all-electron potential
  if lloc=-2 the original recipe of Troullier-Martins
  is used (zero first and second derivatives at r=0)
* lloc&gt;-1 uses the corresponding channel as local PP

NB: if lloc&gt;-1, the corresponding channel must be the last in the
list of wavefunctions appearing after the namelist &amp;inputp
In the relativistic case, if lloc &gt; 0 both the j=lloc-1/2 and
the j=lloc+1/2 wavefunctions must be at the end of the list.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rcloc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rcloc</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em>
Must be specified only if "lloc"=-1, otherwise the
corresponding value of "rcut" is used.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Matching radius (a.u.) for local pseudo-potential (no default).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nlcc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nlcc</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. produce a PP with the nonlinear core
correction of Louie, Froyen, and Cohen
["PRB 26, 1738 (1982)"].
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help new_core_ps -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>new_core_ps</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Status: </em> requires nlcc=.true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. pseudizes the core charge with bessel functions.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rcore -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rcore</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Matching radius (a.u.) for the smoothing of the core charge.
If not specified, the matching radius is determined
by the condition:  rho_core(rcore) = 2*rho_valence(rcore)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tm -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tm</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
* .true. for Troullier-Martins pseudization ["PRB 43, 1993 (1991)"]

* .false. for Rappe-Rabe-Kaxiras-Joannopoulos pseudization
  ["PRB 41, 1227 (1990)", erratum "PRB 44, 13175 (1991)"]
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rho0 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rho0</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Charge at the origin: when the Rappe-Rabe-Kaxiras-Joannopoulos
method with 3 Bessel functions fails, specifying rho0 &gt; 0
may allow to override the problem (using 4 Bessel functions).
Typical values are in the order of 0.01-0.02
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lpaw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lpaw</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. produce a PAW dataset, experimental feature
only for "pseudotype"=3
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help which_augfun -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>which_augfun</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em>
<tt>'AE'</tt> for Vanderbilt-Ultrasoft pseudo-potentials and 'BESSEL' for PAW datasets.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If different from <b>'AE'</b> the augmentation functions are pseudized
before "rmatch_augfun". The pseudization options are:

* <b>'PSQ'</b>        Use Bessel functions to pseudize Q
               from the origin to rmatch_augfun.

These features are available only for PAW:

* <b>'BESSEL'</b>     Use Bessel functions to pseudize the Q.
* <b>'GAUSS'</b>      Use 2 Gaussian functions to pseudize the Q.
* <b>'BG'</b>         Use original Bloechl's recipe with a single gaussian.

Note: if lpaw is true and which_augfun is set to AE real all-
electron charge will be used, which will produce extremely
hard augmentation.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rmatch_augfun -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rmatch_augfun</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5 a.u.
            </li>
<br><li> <em>Status: </em> Used only if which_augfun is different from 'AE'.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Pseudization radius for the augmentation functions. Presently
it has the same value for all L.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rmatch_augfun_nc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rmatch_augfun_nc</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> .false.
            </li>
<br><li> <em>Status: </em> Used only if which_augfun is 'PSQ'.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the augmentation functions are pseudized
from the origin to min(rcut(ns),rcut(ns1)) where ns
and ns1 are the two channels for that Q. In this case
rmatch_augfun is not used.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsave_wfc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsave_wfc</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false. if .not. lpaw, otherwise .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set it to .true. to save all-electron and pseudo wavefunctions
used in the pseudo-potential generation in the UPF file. Only
works for UPFv2 format.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lgipaw_reconstruction -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lgipaw_reconstruction</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set it to .true. to generate pseudo-potentials containing the
additional info required for reconstruction of all-electron
orbitals, used by GIPAW. You will typically need to specify
additional projectors beyond those used in the generation of
pseudo-potentials. You should also specify "file_recon".

All projectors used in the reconstruction must be listed BOTH
in the test configuration after namelist &amp;test AND in the
all-electron configuration (variable 'config', namelist &amp;inputp,
Use negative occupancies for projectors on unbound states). The
core radii in the test configuration should be the same as in
the pseudo-potential generation section and will be used as the
radius of reconstruction. Projectors not used to generate the
pseudo-potential should have zero occupation number.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help use_paw_as_gipaw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>use_paw_as_gipaw</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
When generating a PAW dataset, setting this option to .true. will
save the core all-electron wavefunctions to the UPF file.
The GIPAW reconstruction to be performed using the PAW data and
projectors for the valence wavefunctions.

In the default case, the GIPAW valence wavefunction and projectors
are independent from the PAW ones and must be then specified as
explained above in lgipaw_reconstruction.

Setting this to .true. always implies "lgipaw_reconstruction" = .true.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help author -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>author</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'anonymous'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Name of the author.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_chi -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_chi</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file containing output PP chi functions
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_beta -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_beta</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file containing output PP beta functions
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_qvan -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_qvan</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file containing output PP qvan functions
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_screen -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_screen</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file containing output screening potential
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_core -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_core</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file containing output total and core charge
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_wfcaegen -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_wfcaegen</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file with the all-electron wfc for generation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_wfcncgen -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_wfcncgen</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file with the norm-conserving wfc for generation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_wfcusgen -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_wfcusgen</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> file with the ultra-soft wfc for generation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nwfs -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nwfs</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> number of wavefunctions to be pseudized
                     </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help PP_wfs -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>nls</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Wavefunction label (same as in the all-electron configuration).
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>nns</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Principal quantum number (referred to the PSEUDOPOTENTIAL case;
nns=1 for lowest s, nns=2 for lowest p, and so on).
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>lls</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Angular momentum quantum number.
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>ocs</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Occupation number  (same as in the all-electron configuration).
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>ener</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Energy (Ry) used to pseudize the corresponding state.
If 0.d0, use the one-electron energy of the all-electron state.
Do not use 0.d0 for unbound states!
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>rcut</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Matching radius (a.u.) for norm conserving PP.
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>rcutus</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Matching radius (a.u.) for ultrasoft PP - only for pseudotype=3.
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>jjs</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> The total angular momentum (0.0 is allowed for complete shells).
                        </pre></blockquote>
</ul>   
    

       
    
}


# ------------------------------------------------------------------------
help nconf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nconf</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> the number of configurations to be tested. For iswitch = 4 nconf=2
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help file_pseudo -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>file_pseudo</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Status: </em> ignored if iswitch=3
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
File containing the PP.

* If the file name contains  ".upf" or ".UPF",
the file is assumed to be in UPF format;

* else if the file name contains ".rrkj3" or ".RRKJ3",
the old RRKJ format is first tried;

* otherwise, the old NC format is read.

IMPORTANT: in the latter case, all calculations are done
using the SEMILOCAL form, not the separable nonlocal form.
Use the UPF format if you want to test the separable form!
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {ecutmin ecutmax decut} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>ecutmin, ecutmax, decut</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
decut=5.0 Ry; ecutmin=ecutmax=0Ry
         </li>
<br><li> <em>Status: </em> specify "ecutmin" and "ecutmax" if you want to perform this test
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Parameters (Ry) used for test with a basis set of spherical
Bessel functions j_l(qr) . The hamiltonian at fixed scf
potential is diagonalized for various values of ecut:
"ecutmin", "ecutmin"+"decut", "ecutmin"+2*"decut" ... up to "ecutmax".
This yields an indication of convergence with the
corresponding plane-wave cutoff in solids, and shows
in an unambiguous way if there are "ghost" states
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help rm -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rm</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 30 a.u.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Radius of the box used with spherical Bessel functions.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help configts -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>configts(i), i=1,nconf</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
A string containing the test valence electronic
configuration nc, nc=1,nconf. Same syntax as for "config".
If configts(i) is not set, the electron configuration
is read from the cards following the namelist.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsdts -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>lsdts(i), i=1,nconf</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>See: </em> lsd
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 or 1. It is the value of lsd used in the i-th test.
Allows to make simultaneously spin-polarized and
spin-unpolarized tests.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help frozen_core -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>frozen_core</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. only the core wavefunctions of the first
configuration are calculated. The eigenvalues, orbitals
and energies of the other configurations are calculated
with the core of the first configuration.
The first configuration must be spin-unpolarized.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rcutv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rcutv</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -1.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Cutoff distance (CUT) for the inclusion of LDA-1/2 potential.
                  Needed (mandatory) only if iswitch = 4
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nwfts_1 nwfts_2 nwfts_3 nwfts_4 nwfts_5 nwfts_6 nwfts_7 nwfts_8 nwfts_9 nwfts_10 nwfts_11 nwfts_12 nwfts_13 nwfts_14 nwfts_15 nwfts_16 nwfts_17 nwfts_18 nwfts_20 nwfts_19} -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nwfts</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> number of wavefunctions
                     </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {test_wfs_10 test_wfs_11 test_wfs_12 test_wfs_13 test_wfs_14 test_wfs_15 test_wfs_16 test_wfs_17 test_wfs_18 test_wfs_20 test_wfs_19 test_wfs_1 test_wfs_2 test_wfs_3 test_wfs_4 test_wfs_5 test_wfs_6 test_wfs_7 test_wfs_8 test_wfs_9} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>elts</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>See: </em> nls
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>nnts</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>See: </em> nns
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>llts</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>See: </em> lls
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>octs</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>See: </em> ocs
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>enerts</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> not used
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>rcutts</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> not used
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>rcutusts</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> not used
                        </li>
<br>
</ul><ul>
<li> <em>Variable: </em><big><b>iswts</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> spin index (1 or 2, used in lsda case)
                        </pre></blockquote>
</ul>   
    

    <ul>
<li> <em>Variable: </em><big><b>jjts</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> total angular momentum of the state
                        </pre></blockquote>
</ul>   
    

       
    
}

