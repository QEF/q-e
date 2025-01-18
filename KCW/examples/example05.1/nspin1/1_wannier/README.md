# Initial DFT and Wannier90 calculations, and interface with the KCW code.

* DFT calculation 
 1) Self-consistent field: 
  
    pw.x -in scf.pwi | tee scf.pwo

 2) Non-self-consisten field on a regular k-point mesh (needed by Wannier90)

    pw.x -in nscf.pwi | tee nscf.pwo

* Separate wannierization for occupied and empty orbitals
 1) mv to the occ folder
 
    cd occ

 2) Wannier90 pre-processing

    wannier90.x -pp wann

 3) interface between PW and Wannier90
 
    pw2wannier90.x -in pw2wan.p2wi | tee pw2wan.p2wo

 4) Wannierization of the occupied manifold

    wannier90.x wann

 5) mv to the emp folder

    cd ../emp

 6) Wannier90 pre-processing

    wannier90.x -pp wann

 7) interface between PW and Wannier90

    pw2wannier90.x -in pw2wan.p2wi | tee pw2wan.p2wo

 8) Wannierization of the (lower part of the) empty manifold

    wannier90.x wann

* Interface between PWscf,Wannier90 and KCW

 1) go back to the working dir

    cd ../

 2) create a symbolic link to the rotational matrices produced by Wannier90 

    ./link_wann.sh

 3) Run kcw.x 

    kcw.x -in kc.w2ki | tee kc.w2ko
