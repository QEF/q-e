  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !--------------------------------------------------------------------------
  subroutine rotate_eigenm ( iq_first, nqc, isym, s, invs, irt, &
      rtau, xq, cz1, cz2) 
  !--------------------------------------------------------------------------
  !!
  !!  Here:
  !! 
  !!  1). determine the unitary matrix which gives the transformed eigenmodes
  !!     from the first q in the star to the present q
  !! 
  !!  2). rotate eigenmodes from the first q in the star (cz1) to the current 
  !!      q (cz2)
  !!  
  !!  In rotate_epmat.f90:
  !! 
  !!  3). rotate the electron-phonon matrix from the cartesian representation
  !!     of the first qpoint of the star to the eigenmode representation 
  !!     (using cz1).
  !!  
  !!  4). rotate the electron-phonon matrix from the eigenmode representation
  !!     to the cartesian representation of the qpoint iq in the star (with cz2).
  !! 
  !!     Step 3 requires using the rotated eigenmodes to set the phases 
  !!     (the diagonalizers of the rotated dyn mat will not
  !!     work because of the gages in deg. subspaces and the phases).
  !! 
  !!  W.r.t. the standard method of q2qstar_ph.f90, here I construct the
  !!  unitary matrix Gamma defined in Maradudin & Vosko, RMP 1968.
  !!  In this way all rotations are just zgemm operations and we throw away  
  !!  all nasty indices.
  !!
  !
  !--------------------------------------------------------------------------
  USE kinds,         ONLY : DP
  use io_global,     ONLY : stdout
  use elph2,         ONLY : dynq
  use phcom,         ONLY : nmodes
  USE constants_epw, ONLY : cone, czero, twopi, rydcm1 
  use control_flags, ONLY : iverbosity
  use cell_base,     ONLY : at, bg
  use ions_base,     ONLY : nat, amass, nat, ityp
  implicit none
  !
  integer :: iq_first, nqc, isym
  !  originating q point of the star
  !  total number of qpoints
  !  fractional translation
  !  the symmetry under consideration
  !  ...
  !  degeneracy of the star of q
  !  number of q in the star
  !  number of global (crystal) symmetries
  !  index of this q in the star
  real(kind=DP) :: xq(3), rtau(3,48,nat)
  !  the rotated q vector
  !  the relative position of the rotated atom to the original one
  integer :: s(3,3,48), irt(48,nat), invs(48)
  !  the symmetry operation for the eigenmodes
  !  the rotated of each atom
  !  the folding G-vector
  !  the index of the inverse operation
  complex(kind=DP) :: gamma(nmodes, nmodes), cz1( nmodes, nmodes), &
       cz2(nmodes, nmodes)
  !  the Gamma matrix for the symmetry operation on the dyn mat
  !  the eigenvectors for the first q in the star
  !  the rotated eigenvectors, for the current q in the star
  !
  ! variables for lapack zhpevx
  !
  complex(kind=DP) :: cwork( 2*nmodes ), dynp( nmodes*(nmodes+1)/2 )
  integer :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
  real(kind=DP) :: w1( nmodes ), w2(nmodes), rwork( 7*nmodes )
  !  dynp: complex dynmat packed (upper triangular part for zhpevx)
  !  after hermitian-ization
  !
  real(kind=DP), parameter :: eps = 1.d-10
  !
  ! work variables 
  !
  integer :: imode, jmode, nu, mu, sna, ism1, na, nb
  real(kind=DP) :: wtmp( nmodes ), arg, massfac, scart(3,3)
  complex(kind=DP) :: cfac, dyn1(nmodes, nmodes), dyn2( nmodes, nmodes)
  !
  ! ------------------------------------------------------------------
  ! diagonalize dynq(iq_first) --> w1, cz1
  ! ------------------------------------------------------------------
  !
  !  first divide by the square root of masses 
  !
  DO na = 1, nat
   DO nb = 1, nat
     massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
     dyn1(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
     dynq(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, iq_first) * massfac
   END DO
  END DO
  !
  DO jmode = 1, nmodes
   DO imode = 1, jmode
      dynp (imode + (jmode - 1) * jmode/2 ) = &
      ( dyn1 ( imode, jmode) + conjg ( dyn1 ( jmode, imode) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nmodes, dynp , 0.0, 0.0, &
               0, 0,-1.0, neig, w1, cz1, nmodes, cwork, &
               rwork, iwork, ifail, info)
  !
  IF (iverbosity.eq.1) then
    !
    !  check the frequencies
    !
    DO nu = 1, nmodes
      IF ( w1 (nu) .gt. 0.d0 ) then
         wtmp(nu) =  sqrt(abs( w1 (nu) ))
      ELSE
         wtmp(nu) = -sqrt(abs( w1 (nu) ))
      ENDIF
    ENDDO
    WRITE ( stdout, '(5x,"Frequencies of the matrix for the first q in the star")') 
    WRITE ( stdout, '(6(2x,f10.5))' ) (wtmp(nu)*rydcm1, nu = 1, nmodes)
    !
  ENDIF
  !
  !  here dyn1 is dynq(:,:,iq_first) after dividing by the masses
  !  (the true dyn matrix D_q)
  !
  !  -----------------------------------------------------------------------
  !  the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)   
  !  -----------------------------------------------------------------------
  !
  !  I have built the matrix by following step-by-step q2star_ph.f90 and
  !  rotate_and_add_dyn.f90
  !
  ism1 = invs (isym)
  !
  !  the symmetry matrix in cartesian coordinates 
  !  (so that we avoid going back and forth with the dynmat)  
  !  note the presence of both at and bg in the transform!
  !
  scart = dble ( s ( :, :, ism1) )
  scart = matmul ( matmul ( bg, scart), transpose (at) )
  ! 
  gamma = czero
  DO na = 1, nat
    !
    ! the corresponding of na in {S|v}
    sna = irt (isym, na)
    !
    ! cfac = exp[iSq*(tau_K - {S|v} tau_k)]   (Maradudin&Vosko RMP Eq. 2.33)
    ! [v can be ignored since it cancels out, see endnotes. xq is really Sq]
    ! rtau(:,isym,na) = s(:,:,invs(isym)) * tau(:, na) - tau(:,irt(isym,na))) (cartesian)
    !
    arg = twopi * dot_product (xq, rtau (:, isym, na)) 
    cfac = dcmplx (cos(arg),-sin(arg))
    !
    !  the submatrix (sna,na) contains the rotation scart 
    !
    gamma ( 3*(sna-1)+1:3*sna, 3*(na-1)+1:3*na ) = cfac * scart 
    !
  ENDDO
  !
  !  possibly run some consistency checks
  !
  IF ( iverbosity .eq. 1 ) then
    !
    !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5) 
    ! 
    CALL zgemm ('n', 'n', nmodes, nmodes, nmodes , cone  , gamma  , &
           nmodes, dyn1 , nmodes, czero , dyn2, nmodes)
    CALL zgemm ('n', 'c', nmodes, nmodes, nmodes , cone  , dyn2, &
           nmodes, gamma, nmodes, czero , dyn1   , nmodes)
    !
    DO jmode = 1, nmodes
     DO imode = 1, jmode
       dynp (imode + (jmode - 1) * jmode/2 ) = &
        ( dyn1 ( imode, jmode) + conjg ( dyn1 ( jmode, imode) ) ) / 2.d0
     ENDDO
    ENDDO
    !
    CALL zhpevx ('V', 'A', 'U', nmodes, dynp , 0.0, 0.0, &
                 0, 0,-1.0, neig, w2, cz2, nmodes, cwork, &
                 rwork, iwork, ifail, info)
    !
    !  check the frequencies
    !
    DO nu = 1, nmodes
      IF ( w2 (nu) .gt. 0.d0 ) then
         wtmp(nu) =  sqrt(abs( w2 (nu) ))
      ELSE
         wtmp(nu) = -sqrt(abs( w2 (nu) ))
      ENDIF
    ENDDO
    WRITE ( stdout, '(6(2x,f10.5) )' ) (wtmp(nu)*rydcm1*0.12398d0, nu = 1, nmodes)
    !
  ENDIF
  !
  !
  !  -----------------------------------------------------------------------
  !  transformation of the eigenvectors: e_{Sq} = gamma * e_q  (MV eq. 2.36)
  !  -----------------------------------------------------------------------
  !
  CALL zgemm ('n', 'n', nmodes, nmodes, nmodes , cone  , gamma  , &
       nmodes, cz1 , nmodes, czero , cz2, nmodes)
  !
  !   now check that cz2 diagonalizes dyn2 (after dividing by the masses)
  !
  DO na = 1, nat
   DO nb = 1, nat
     massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
     dyn2(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
       dynq(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, nqc) * massfac
   END DO
  END DO
  !
  ! the rotated matrix and the one read from file
  IF (iverbosity.eq.1) write (6,'(2f15.10)') dyn2-dyn1 
  !
  ! here I have checked that the matrix rotated with gamma 
  ! is perfectly equal to the one read from file for this q in the star
  !
  CALL zgemm ('c', 'n', nmodes, nmodes, nmodes, cone, cz2, &
      nmodes, dyn2, nmodes, czero, dyn1, nmodes)
  CALL zgemm ('n', 'n', nmodes, nmodes, nmodes, cone, dyn1, &
     nmodes, cz2, nmodes, czero, dyn2, nmodes)
  !
  DO nu = 1, nmodes
    w2(nu) = abs(dyn2(nu,nu))
    DO mu = 1, nmodes
    IF ( mu.ne.nu .and. abs(dyn2(mu,nu)).gt.eps ) call errore &
      ('rotate_eigenm','problem with rotated eigenmodes',0)
    ENDDO
  ENDDO
  !
  IF (iverbosity.eq.1) then
    !
    !  a simple check on the frequencies
    !
    DO nu = 1, nmodes
      IF ( w2 (nu) .gt. 0.d0 ) then
         wtmp(nu) =  sqrt(abs( w2 (nu) ))
      ELSE
         wtmp(nu) = -sqrt(abs( w2 (nu) ))
      ENDIF
    ENDDO
    WRITE ( stdout, '(5x,"Frequencies of the matrix for the current q in the star")') 
    WRITE ( stdout, '(6(2x,f10.5))' ) (wtmp(nu)*rydcm1, nu = 1, nmodes)
    !
  ENDIF
  !
  end subroutine rotate_eigenm
