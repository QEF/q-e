
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help prefix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>prefix</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
prefix of files saved by program pw.x
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help outdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>outdir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em>
value of the <tt>ESPRESSO_TMPDIR</tt> environment variable if set;
current directory ('./') otherwise
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
directory containing the input data, i.e. the same as in pw.x
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filplot -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filplot</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
file "filplot" contains the quantity selected by plot_num
(can be saved for further processing)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help plot_num -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>plot_num</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Selects what to save in filplot:

   0  = electron (pseudo-)charge density

   1  = total potential V_bare + V_H + V_xc

   2  = local ionic potential V_bare

   3  = local density of states at specific energy or grid of energies
        (number of states per volume, in bohr^3, per energy unit, in Ry)

   4  = local density of electronic entropy

   5  = STM images
        Tersoff and Hamann, "PRB 31, 805 (1985)"

   6  = spin polarization (rho(up)-rho(down))

   7  = contribution of a selected wavefunction to the
        (pseudo-)charge density. For norm-conserving PPs,
        |psi|^2 (psi=selected wavefunction). Noncollinear case:
        contribution of the given state to the charge or
        to the magnetization along the direction indicated
        by spin_component (0 = charge, 1 = x, 2 = y, 3 = z )

   8  = electron localization function (ELF)

   9  = charge density minus superposition of atomic densities

   10 = integrated local density of states (ILDOS)
        from "emin" to "emax" (emin, emax in eV)
        if "emax" is not specified, "emax"=E_fermi

   11 = the V_bare + V_H potential

   12 = the sawtooth electric field potential (if present)

   13 = the noncollinear magnetization.

   17 = all-electron valence charge density
        can be performed for PAW calculations only
        requires a very dense real-space grid!

   18 = The exchange and correlation magnetic field in the noncollinear case

   19 = Reduced density gradient
        ( J. Chem. Theory Comput. 7, 625 (2011), "doi:10.1021/ct100641a" )
        Set the isosurface between 0.3 and 0.6 to plot the
        non-covalent interactions (see also plot_num = 20)

   20 = Product of the electron density (charge) and the second
        eigenvalue of the electron-density Hessian matrix;
        used to colorize the RDG plot (plot_num = 19)

   21 = all-electron charge density (valence+core).
        For PAW calculations only; requires a very dense real-space grid.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help spin_component -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = total charge (default value),
1 = spin up charge,
2 = spin down charge.
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = spin averaged potential (default value),
1 = spin up potential,
2 = spin down potential.
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variables: </em><big><b>spin_component(i), i=1,2</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Status: </em> OPTIONAL
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
<b>Noncollinear case only:</b>
plot the contribution of the given state to the charge
or to the magnetization along the direction indicated
by spin_component:
        0 = charge (default),
        1 = x,
        2 = y,
        3 = z.

Ignored in unpolarized or LSDA case

When both spin_component(1) and spin_component(2) are specified,
all components in the range [spin_component(1),spin_component(2)]
are plotted.
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
for LSDA case only: plot the contribution to ILDOS of
0 = spin-up + spin-down (default)
1 = spin-up   only
2 = spin-down only
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = absolute value (default value)
1 = x component of the magnetization
2 = y component of the magnetization
3 = z component of the magnetization
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = total charge (default value),
1 = spin up charge,
2 = spin down charge.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help emin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>emin</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> e_fermi
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
lower boundary of energy grid (in eV).

Defaults to Fermi energy.
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>emin</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
lower energy boundary (in eV)
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help emax -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>emax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> OPTIONAL
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
upper boundary of energy grid (in eV).

If not specified, LDOS is computed just for energy "emin"
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variable: </em><big><b>emax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
upper energy boundary (in eV),
i.e. compute ILDOS from "emin" to "emax"
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help delta_e -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>delta_e</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1
               </li>
<br><li> <em>Status: </em> OPTIONAL
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
spacing of energy grid (in eV).
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help degauss_ldos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>degauss_ldos</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> degauss (converted to eV)
               </li>
<br><li> <em>Status: </em> OPTIONAL
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
broadening of energy levels for LDOS (in eV).

Defaults to broadening degauss specified for electronic smearing
in pw.x calculation.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sample_bias -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sample_bias</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
the bias of the sample (Ry) in stm images
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help kpoint -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>kpoint(i), i=1,2</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Unpolarized and noncollinear case:
        k-point to be plotted
LSDA:
        k-point and spin polarization to be plotted
        (spin-up and spin-down correspond to different k-points!)

When both kpoint(1) and kpoint(2) are specified, all
kpoints in the range [kpoint(1),kpoint(2)] are plotted.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help kband -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>kband(i), i=1,2</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Band to be plotted.

When both kband(1) and kband(2) are specified, all
bands in the range [kband(1),kband(2)] are plotted.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsign -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsign</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if true and k point is Gamma, plot |psi|^2 sign(psi)
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nfile -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nfile</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Status: </em> OPTIONAL
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
the number of data files to read
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filepp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>filepp(i), i=1,nfile</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> filepp(1)=filplot
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
nfile = 1 : file containing the quantity to be plotted
nfile &gt; 1 : see "weight"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help weight -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>weight(i), i=1,nfile</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> weight(1)=1.0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
weighing factors: assuming that rho(i) is the quantity
read from filepp(i), the quantity that will be plotted is:

weight(1)*rho(1) + weight(2)*rho(2) + weight(3)*rho(3) + ...
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help iflag -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>iflag</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = 1D plot of the spherical average
1 = 1D plot
2 = 2D plot
3 = 3D plot
4 = 2D polar plot on a sphere
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help output_format -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>output_format</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
(ignored on 1D plot)

0  = format suitable for gnuplot   (1D)

1  = format suitable for contour.x (2D)

2  = format suitable for plotrho   (2D)

3  = format suitable for XCRYSDEN  (2D or user-supplied 3D region)

4  = format suitable for gOpenMol  (3D)
     (formatted: convert to unformatted *.plt)

5  = format suitable for XCRYSDEN  (3D, using entire FFT grid)

6  = format as gaussian cube file  (3D)
     (can be read by many programs)

7  = format suitable for gnuplot   (2D) x, y, f(x,y)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fileout -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fileout</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> standard output
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
name of the file to which the plot is written
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help interpolation -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>interpolation</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'fourier'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Type of interpolation:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fourier'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bspline'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> (EXPERIMENTAL)
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help e1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>e1(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vector which determines the plotting line (in alat units)
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help x0 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>x0(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vector, origin of the line (in alat units)
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variables: </em><big><b>x0(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vector, origin of the plane (in alat units)
               </pre></blockquote>
</ul>      
      

      <ul>
<li> <em>Variables: </em><big><b>x0(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vector, origin of the parallelepiped

"x0" is in alat units !
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nx -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nx</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of points in the line:

rho(i) = rho( x0 + e1 * (i-1)/(nx-1) ), i=1, nx
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {e1 e2} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>e1(i), 
		e2(i), 
		 
		  i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vectors which determine the plotting plane (in alat units)

BEWARE: <b>e1</b> and <b>e2</b> must be orthogonal
               </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {nx ny} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nx, ny</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of points in the plane:

rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)
                + e2 * (j-1)/(ny-1) ), i=1,nx ; j=1,ny
               </pre></blockquote>
</ul>
    

    <ul>
<li> <em>Variables: </em><big><b>nx, ny</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of points in the polar plane:

phi(i)   = 2 pi * (i - 1)/(nx-1), i=1, nx
theta(j) =   pi * (j - 1)/(ny-1), j=1, ny
               </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {e1 e2 e3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>e1(i), 
		e2(i), 
		e3(i), 
		 
		  i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
3D vectors which determine the plotting parallelepiped
(if present, must be orthogonal)

"e1", "e2", and "e3" are in alat units !
               </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {nx ny nz} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nx, ny, nz</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of points in the parallelepiped:

rho(i,j,k) = rho( x0 + e1 * (i-1)/nx
                  + e2 * (j-1)/ny
                  + e3 * (k-1)/nz ),
             i = 1, nx ; j = 1, ny ; k = 1, nz

- If "output_format" = 3 (XCRYSDEN), the above variables
  are used to determine the grid to plot.

- If "output_format" = 5 (XCRYSDEN), the above variables
  are ignored, the entire FFT grid is written in the
  XCRYSDEN format - works for any crystal axis (VERY FAST)

- If "e1", "e2", "e3", "x0" are present, and "e1", "e2", "e3" are parallel
  to xyz and parallel to crystal axis, a subset of the
  FFT grid that approximately covers the parallelepiped
  defined by "e1", "e2", "e3", "x0", is written (presently only
  if "output_format" = 4, i.e. gopenmol format) - works only
  if the crystal axis are parallel to xyz

- Otherwise, the required 3D grid is generated from the
  Fourier components (may be VERY slow)
               </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help radius -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>radius</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Radius of the sphere (alat units), centered at (0,0,0)
               </pre></blockquote>
</ul>      
      
}

