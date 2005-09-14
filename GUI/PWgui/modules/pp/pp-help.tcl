
help  prefix -vartype character -helpfmt txt2html -helptext {
     prefix of files saved by program pwscf
}

help  outdir -vartype character -helpfmt txt2html -helptext {
     temporary directory where pwscf files resides
}
  
help  filplot  -vartype character -helpfmt txt2html -helptext {
    file "filplot" contains the quantity selected by plot_num
    (can be saved for further processing)
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

help nfile -vartype integer -helpfmt txt2html -helptext {
            the number of data files (OPTIONAL)
<p> ( default = 1 )
}

help filepp -vartype character -helpfmt txt2html -helptext {
    nfile=1: file containing the quantity to be plotted
    nfile>1: see "weight"
    (default: filepp(1)=filplot)
}
help weight -vartype real -helpfmt txt2html -helptext {
            Weighing factors: assuming that rho(i) is the quantity read
            from filepp(i), the quantity that will be plotted is:
            weight(1)*rho(1) + weight(2)*rho(2) + weight(3)*rho(3)+...
            (OPTIONAL: default weight(1)=1.0)

            BEWARE: atomic coordinates are read from the first file;
            if their number is different for different files,
            the first file must have the largest number of atoms
<p> (  default : weight(1) = 1.0)
}


help iflag -vartype integer -helpfmt txt2html -helptext {
          1 1D plot of spherical average
          1 if a 1D plot is required (DEFAULT)
          2 if a 2D plot is required
          3 if a 3D plot is required
          4 if a 2D polar plot on a sphere is required
<p> ( default = 1 )
}

help output_format -vartype integer -helpfmt txt2html -helptext {
         (ignored on 1D plot)

          0  format suitable for gnuplot   (1D) (DEFAULT)
          1  format suitable for contour.x (2D)
          2  format suitable for plotrho   (2D)
          3  format suitable for XCRYSDEN  (1D, 2D, 3D)
          4  format suitable for gOpenMol  (3D)
             (formatted: convert to unformatted *.plt)
          5  format suitable for XCRYSDEN  (3D)
          6  Gaussian cube-file format (3D)
             (can be read by many programs)
<p> ( default = 0 )
}

set _volume {
  IF iflag = 1
                    REQUIRED:
                e1  3D vector which determines the plotting line
                x0  3D vector, origin of the line
                nx  number of points in the line:
                    rho(i) = rho( x0 + e1 * (i-1)/(nx-1) ), i=1, nx
  
  ELSEIF iflag = 2
                    REQUIRED:
            e1, e2  3D vectors which determine the plotting plane
                    (must be orthogonal)
                x0  3D vector, origin of the plane
            nx, ny  number of points in the plane:
                    rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)
                                       + e2 * (j-1)/(ny-1) ), i=1,nx ; j=1,ny
  
  ELSEIF iflag = 3
                    OPTIONAL:
        e1, e2, e3  3D vectors which determine the plotting parallelepiped
                    (if present, must be orthogonal)
                x0  3D vector, origin of the parallelepiped
          nx,ny,nz  number of points in the parallelepiped:
                    rho(i,j,k) = rho( x0 + e1 * (i-1)/(nx-1)
                                         + e2 * (j-1)/(ny-1)
                                         + e3 * (k-1)/(nz-1) ),
                                 i = 1, nx ; j = 1, ny ; k = 1, nz
  
                  - If output_format = 3 (XCRYSDEN), the above variables
                    are used to determine the grid to plot. 
                  - If output_format = 5 (XCRYSDEN), the above variables
                    are ignored, the entire FFT grid is written in the
                    XCRYSDEN format - works for any crystal axis (VERY FAST)
                  - If e1, e2, e3, x0 are present, e1 e2 e3 are parallel
                    to xyz and parallel to crystal axis, a subset of the
                    FFT grid that approximately covers the parallelepiped
                    defined by e1, e2, e3, x0, is written (presently only
                    in gopenmol "formatted" file format) - works only
                    if the crystal axis are parallel to xyz
                  - Otherwise, the required 3D grid is generated from the
                    Fourier components (may be VERY slow)
}

set _e {
     3D vectors which determine the plotting line (e1), plane (e1,e2), or 
     parallelepiped (e1, e2, e3)
}
foreach e {e1 e2 e3} {
    help $e -vartype real -helpfmt txt2html -helptext "$_e\n\n$_volume"
}

set _x {
    3D vector, origin of either line (1D plot), plane (2D plot), 
    or parallelepiped (3D plot)
}
help x0 -vartype real -helpfmt txt2html -helptext "$_x\n\n$_volume"

set _n {
    number of points for line (nx), plane (nx,ny), or parallelepiped (nx,ny,nz)
}
foreach n {nx ny nz} {
    help $n -vartype real -helpfmt txt2html -helptext "$_n\n\n$_volume"
}

help radius -vartype real -helpfmt txt2html -helptext {
    radius    Radius of the sphere (alat units), centered at (0,0,0)
    nx, ny    number of points in the polar plane:
                      phi(i)  = 2 pi * (i - 1)/(nx-1), i=1, nx
                      theta(j)=   pi * (j - 1)/(ny-1), j=1, ny
}

help fileout -vartype character -helpfmt txt2html -helptext {  
    name of the file to which the plot is written
<p> ( default = standard output )
}

