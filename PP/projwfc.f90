! Copyright (C) 2001-2003 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!----------------------------------------------------------------------- 
program projwfc 
  !----------------------------------------------------------------------- 
  ! 
  ! projects wavefunctions onto orthogonalized atomic wavefunctions, 
  ! calculates Lowdin charges, spilling parameter, projected DOS 
  ! (separated into up and down components for lSDA)
  !
  ! Input (namelist &inputpp ... &end):                   Default value
  !
  !    prefix        prefix of input file produced by pw.x    'pwscf'
  !                    (wavefunctions are needed)
  !    outdir        directory containing the input file       ./
  !    ngauss        type of gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            = 99  Fermi-Dirac function
  !    degauss       gaussian broadening, Ry (not eV!)          0.0
  !    Emin, Emax    min, max energy (eV) for DOS plot          band extrema
  !    DeltaE        energy grid step (eV)                      none
  !    filpdos       prefix for output files containing PDOS(E) prefix
  !
  ! Output:
  !
  !   Symmetrized projections are written to standard output
  !   The total DOS and the sum of projected DOS are written to file 
  !   "filpdos".pdos_tot.
  !    The format for the spin-unpolarized case is
  !        E DOS(E) PDOS(E)
  !    The format for the spin-polarized case is
  !        E DOSup(E) DOSdw(E)  PDOSup(E) PDOSdw(E) 
  !   
  !   Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   (one file per atomic wavefunction found in the pseudopotential file) 
  !    The format for the spin-unpolarized case is
  !        E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
  !    where LDOS = \sum m=1,2l+1 PDOS_m(E)
  !    and PDOS_m(E) = projected DOS on atomic wfc with component m
  !    The format for the spin-polarized case is as above with two 
  !    components for both  LDOS(E) and PDOS_m(E)
  !
  !   All DOS(E) are in states/eV plotted vs E in eV
  !
  ! Important notice:
  !
  !    The tetrahedron method is presently not implemented.
  !    Gaussian broadening is used in all cases:
  !    - if degauss is set to some value in namelist &inputpp, that value
  !      (and the optional value for ngauss) is used
  !    - if degauss is NOT set to any value in namelist &inputpp, the 
  !      value of degauss and of ngauss are read from the input data
  !      file (they will be the same used in the pw.x calculations)
  !    - if degauss is NOT set to any value in namelist &inputpp, AND
  !      there is no value of degauss and of ngauss in the input data
  !      file, degauss=DeltaE (in Ry) and ngauss=0 will be used
  ! Obsolete variables, ignored:
  !   io_choice
  !   smoothing
  ! 
#include "f_defs.h" 
  USE io_global,  ONLY : stdout 
  USE constants,  ONLY : rytoev 
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : degauss, ngauss, lgauss
  use io_files,   only : nd_nmbr, prefix, tmp_dir 
#ifdef __PARA 
  use para, only : me, mypool
  use mp,   only : mp_bcast      
#endif 
  implicit none 
  character (len=80) :: filpdos, io_choice 
  character(len=256) :: outdir 
  real (kind=DP)     :: Emin, Emax, DeltaE, degauss1, smoothing
  integer :: ngauss1, ios, ionode_id = 0  
  namelist / inputpp / outdir, prefix, ngauss, degauss,&
             Emin, Emax, DeltaE, io_choice, smoothing

  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr, ilen
  INTEGER, EXTERNAL   :: iargc

  ! 
  call start_postproc (nd_nmbr) 
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf' 
  outdir = './' 
  filpdos= ' '
  Emin   =-1000000. 
  Emax   =+1000000. 
  DeltaE = 0.01 
  ngauss = 0
  degauss=0.d0
#ifdef __PARA 
  if (me == 1)  then 
#endif 
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                   & ' not found' , ierr )
        !
     END IF
     !
  END DO

  read (5, inputpp, err = 200, iostat = ios) 
200 call errore ('projwfc', 'reading inputpp namelist', abs (ios) ) 
  ! 
  tmp_dir = trim(outdir) 
  ! save the value of degauss and ngauss: they are read from file
  ngauss1 = ngauss
  degauss1=degauss
  ! 
#ifdef __PARA 
  end if 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix,  ionode_id ) 
  CALL mp_bcast( ngauss1, ionode_id ) 
  CALL mp_bcast( degauss1,ionode_id ) 
  CALL mp_bcast( DeltaE,  ionode_id ) 
  CALL mp_bcast( Emin, ionode_id ) 
  CALL mp_bcast( Emax, ionode_id ) 
#endif 
  ! 
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  call read_file 
  call openfil_pp
  ! 
  !   decide Gaussian broadening
  !
  if (degauss1.ne.0.d0) then
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  else if (lgauss) then
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  else
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  end if
  !
  call projwave ( ) 
  !
#ifdef __PARA 
  if (me == 1 .AND. mypool == 1) then 
#endif 
  if ( filpdos == ' ') filpdos = prefix
  call partialdos (Emin, Emax, DeltaE, filpdos)
#ifdef __PARA 
  end if
#endif 
  !
  call stop_pp 
  ! 
end program projwfc 

module projections
  USE kinds, only : DP 

  type wfc_label 
     integer na, n, l, m 
  end type wfc_label 
  type(wfc_label), allocatable :: nlmchi(:) 

  real (kind=DP), allocatable :: proj (:,:,:)

end module projections

!----------------------------------------------------------------------- 
subroutine projwave  ( )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout 
  USE atom 
  USE ions_base, ONLY : nat, ityp, atm, ntyp => nsp
  USE basis,     ONLY : natomwfc
  use cell_base 
  use constants, only: rytoev 
  use gvect 
  use klist, only: xk, nks, nkstot, nelec
  use ldaU 
  use lsda_mod, only: nspin, isk, current_spin
  use symme, only: nsym, irt 
  use wvfct 
  use uspp, only: nkb, vkb
  use becmod,   only: becp, rbecp
  use io_files, only: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc 
  use wavefunctions_module, only: evc 
#ifdef __PARA 
  use para, only: me, mypool
#endif
  !
  use projections
  implicit none 
  !
  integer :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, lm, nwfc,& 
       nwfc1, lmax_wfc, is 
  real(kind=DP), allocatable :: e (:)
  complex(kind=DP), allocatable :: wfcatom (:,:)
  complex(kind=DP), allocatable :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ... 
  real   (kind=DP), allocatable ::roverlap(:,:),          rwork1(:),rproj0(:,:)
  ! ... or for gamma-point. 
  real(kind=DP), allocatable :: charges(:,:,:), proj1 (:)
  real(kind=DP) :: psum, totcharge(nspinx)
  integer, allocatable :: index(:) 
  ! 
  ! 
  ! 
  WRITE( stdout, '(/5x,"Calling projwave .... ")') 
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")') 
  END IF
  ! 
  if (lmax_wfc > 3) call errore ('projwave', 'l > 3 not yet implemented', 1) 
  ! 
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1 
  ! 
  call d_matrix (d1, d2, d3)   
  ! 
  ! fill structure nlmchi 
  ! 
  allocate (nlmchi(natomwfc)) 
  nwfc=0 
  lmax_wfc = 0 
  do na = 1, nat 
     nt = ityp (na) 
     do n = 1, nchi (nt) 
        if (oc (n, nt) >= 0.d0) then 
           l = lchi (n, nt) 
           lmax_wfc = max (lmax_wfc, l ) 
           do m = 1, 2 * l + 1 
              nwfc=nwfc+1 
              nlmchi(nwfc)%na = na 
              nlmchi(nwfc)%n  =  n 
              nlmchi(nwfc)%l  =  l 
              nlmchi(nwfc)%m  =  m 
           enddo 
        endif 
     enddo 
  enddo 
  if (nwfc /= natomwfc) call errore ('projwave', 'wrong # of atomic wfcs?', 1) 
  ! 

! 
  allocate(proj (natomwfc, nbnd, nkstot) ) 
  proj   = 0.d0 
  if (.not. lda_plus_u) allocate(swfcatom (npwx , natomwfc ) ) 
  allocate(wfcatom (npwx, natomwfc) ) 
  allocate(overlap (natomwfc, natomwfc) ) 
  overlap= (0.d0,0.d0) 
  if ( gamma_only ) then 
     allocate(roverlap (natomwfc, natomwfc) ) 
     roverlap= 0.d0 
     allocate (rbecp (nkb,natomwfc)) 
  else
     allocate ( becp (nkb,natomwfc)) 
  end if 
  allocate(e (natomwfc) ) 
  ! 
  !    loop on k points 
  ! 
  call init_us_1 
  call init_at_1 
  ! 
  do ik = 1, nks 
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) 
     call davcio (evc, nwordwfc, iunwfc, ik, - 1) 
 
     call atomic_wfc (ik, wfcatom) 
 
     call init_us_2 (npw, igk, xk (1, ik), vkb) 
 
     if ( gamma_only ) then 
        call pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, wfcatom, npwx, rbecp, nkb)   
     else 
        call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom) 
     end if 
 
     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j> 
     ! 
     if ( gamma_only ) then 
        call pw_gemm ('Y', natomwfc, natomwfc, npw, wfcatom, npwx, swfcatom, &
             npwx, roverlap, natomwfc) 
        overlap(:,:)=cmplx(roverlap(:,:),0.d0)
        ! TEMP: diagonalization routine for real matrix should be used instead 
     else 
        call ccalbec (natomwfc, npwx, npw, natomwfc, overlap, wfcatom, swfcatom) 
     end if 
 
     ! 
     ! calculate O^{-1/2} 
     ! 
     allocate(work (natomwfc, natomwfc) ) 
     call cdiagh (natomwfc, overlap, natomwfc, e, work) 
     do i = 1, natomwfc 
        e (i) = 1.d0 / dsqrt (e (i) ) 
     enddo 
     do i = 1, natomwfc 
        do j = i, natomwfc 
           overlap (i, j) = (0.d0, 0.d0) 
           do k = 1, natomwfc 
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * conjg (work (i, k) ) 
           enddo 
           if (j /= i) overlap (j, i) = conjg (overlap (i, j)) 
        enddo 
     enddo 
     deallocate (work) 
     ! 
     ! calculate wfcatom = O^{-1/2} \hat S | phi> 
     ! 
     if ( gamma_only ) then 
        roverlap(:,:)=real(overlap(:,:),DP)
        ! TEMP: diagonalization routine for real matrix should be used instead 
        call DGEMM ('n', 't', 2*npw, natomwfc, natomwfc, 1.d0 , & 
             swfcatom, 2*npwx,  roverlap, natomwfc, 0.d0, wfcatom, 2*npwx) 
     else 
        call ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , & 
             swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx) 
     end if 
 
     ! 
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j> 
     ! 
     if ( gamma_only ) then 
        allocate(rproj0(natomwfc,nbnd), rwork1 (nbnd) ) 
        call pw_gemm ('Y', natomwfc, nbnd, npw, wfcatom, npwx, evc, npwx, rproj0, natomwfc) 
     else 
        allocate(proj0(natomwfc,nbnd), work1 (nbnd) ) 
        call ccalbec (natomwfc, npwx, npw, nbnd, proj0, wfcatom, evc) 
     end if 
     ! 
     ! symmetrize the projections 
     ! 
     do nwfc = 1, natomwfc 
        ! 
        !  atomic wavefunction nwfc is on atom na 
        ! 
        na= nlmchi(nwfc)%na 
        n = nlmchi(nwfc)%n 
        l = nlmchi(nwfc)%l 
        m = nlmchi(nwfc)%m 
        ! 
        do isym = 1, nsym 
           nb = irt (isym, na) 
           do nwfc1 =1, natomwfc 
              if (nlmchi(nwfc1)%na == nb             .and. & 
                   nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. & 
                   nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. & 
                   nlmchi(nwfc1)%m == 1 ) go to 10 
           end do 
           call errore('projwave','cannot symmetrize',1) 
10         nwfc1=nwfc1-1 
           ! 
           !  nwfc1 is the first rotated atomic wfc corresponding to nwfc 
           ! 
           if ( gamma_only ) then 
              if (l == 0) then 
                 rwork1(:) = rproj0 (nwfc1 + 1,:) 
              else if (l == 1) then  
                 rwork1(:) = 0.d0   
                 do m1 = 1, 3   
                    rwork1(:) = rwork1(:) + d1 (m1, m, isym) * rproj0 (nwfc1 + m1,:) 
                 enddo 
              else if (l == 2) then  
                 rwork1(:) = 0.d0   
                 do m1 = 1, 5   
                    rwork1(:) = rwork1(:) + d2 (m1, m, isym) * rproj0 (nwfc1 + m1,:) 
                 enddo 
              else if (l == 3) then  
                 rwork1(:) = 0.d0   
                 do m1 = 1, 7   
                    rwork1(:) = rwork1(:) + d3 (m1, m, isym) * rproj0 (nwfc1 + m1,:) 
                 enddo 
              endif 
              do ibnd = 1, nbnd 
                 proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                      rwork1(ibnd) * rwork1(ibnd) / nsym 
              enddo 
           else 
              if (l == 0) then 
                 work1(:) = proj0 (nwfc1 + 1,:) 
              else if (l == 1) then  
                 work1(:) = 0.d0   
                 do m1 = 1, 3   
                    work1(:) = work1(:) + d1 (m1, m, isym) * proj0 (nwfc1 + m1,:) 
                 enddo 
              else if (l == 2) then  
                 work1(:) = 0.d0   
                 do m1 = 1, 5   
                    work1(:) = work1(:) + d2 (m1, m, isym) * proj0 (nwfc1 + m1,:) 
                 enddo 
              else if (l == 3) then  
                 work1(:) = 0.d0   
                 do m1 = 1, 7   
                    work1(:) = work1(:) + d3 (m1, m, isym) * proj0 (nwfc1 + m1,:) 
                 enddo 
              endif 
              do ibnd = 1, nbnd 
                 proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                      work1(ibnd) * conjg (work1(ibnd)) / nsym 
              enddo 
           end if 
        enddo 
     enddo 
     if ( gamma_only ) then 
        deallocate (rwork1) 
        deallocate (rproj0) 
     else 
        deallocate (work1) 
        deallocate (proj0) 
     end if 
     ! on k-points 
  enddo 
  !
  deallocate (e)
  if ( gamma_only ) then 
     deallocate (roverlap) 
     deallocate (rbecp)
  else
     deallocate ( becp)
  end if 
  deallocate (overlap) 
  deallocate (wfcatom) 
  if (.not. lda_plus_u) deallocate (swfcatom) 
  ! 
#ifdef __PARA 
  ! 
  !   vectors et and proj are distributed across the pools 
  !   collect data for all k-points to the first pool
  ! 
  call poolrecover (et, nbnd, nkstot, nks) 
  call poolrecover (proj, nbnd * natomwfc, nkstot, nks) 
  ! 
  if (me == 1 .AND. mypool == 1) then 
#endif 
     ! 
     ! write on the standard output file 
     ! 
     WRITE( stdout,'(/"Projection on atomic states:"/)') 
     do nwfc = 1, natomwfc 
        WRITE(stdout,1000) & 
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m 
     end do
1000 FORMAT (5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")") 
     ! 
     allocate(index (natomwfc), proj1 (natomwfc) ) 
     do ik = 1, nkstot 
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3) 
        do ibnd = 1, nbnd 
           WRITE( stdout, '(5x,"e = ",f14.10," eV")') et (ibnd, ik) * rytoev 
           ! 
           ! sort projections by magnitude, in decreasing order 
           ! 
           do nwfc = 1, natomwfc 
              index (nwfc) = 0 
              proj1 (nwfc) = - proj (nwfc, ibnd, ik) 
           end do
           call hpsort (natomwfc, proj1, index) 
           ! 
           !  only projections that are larger than 0.001 are written 
           ! 
           do nwfc = 1, natomwfc 
              proj1 (nwfc) = - proj1(nwfc) 
              if ( abs (proj1(nwfc)) < 0.001 ) go to 20 
           end do
           nwfc = natomwfc + 1 
20         nwfc = nwfc -1 
           ! 
           ! fancy (?!?) formatting 
           ! 
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i3,"]+"))') & 
                (proj1 (i), index(i), i = 1, min(5,nwfc)) 
           do j = 1, (nwfc-1)/5 
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i3,"]+"))') & 
                   (proj1 (i), index(i), i = 5*j+1, min(5*(j+1),nwfc)) 
           end do
           psum = 0.d0 
           do nwfc = 1, natomwfc 
              psum = psum + proj (nwfc, ibnd, ik) 
           end do
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum 
           ! 
        enddo
     enddo
     deallocate (index, proj1) 
     ! 
     ! estimate partial charges (Loewdin) on each atom 
     ! 
     allocate ( charges (nat, 0:lmax_wfc, nspin ) ) 
     charges = 0.0 
     do ik = 1, nkstot 
        if ( nspin == 1 ) then
           current_spin = 1
        else if ( nspin == 2 ) then
           current_spin = isk ( ik )
        else
           call errore ('projwfc',' non collinear case not implemented ',1)
        end if
        do ibnd = 1, nbnd 
           do nwfc = 1, natomwfc 
              na= nlmchi(nwfc)%na 
              l = nlmchi(nwfc)%l 
              charges(na,l,current_spin) = charges(na,l,current_spin) + &
                   wg (ibnd,ik) * proj (nwfc, ibnd, ik) 
           enddo 
        end do 
     end do 
     ! 
     WRITE( stdout, '(/"Lowdin Charges: "/)') 
     ! 
     do na = 1, nat 
        do is = 1, nspin
           totcharge(is) = SUM(charges(na,0:lmax_wfc,is))
        end do
        if ( nspin == 1) then
           WRITE( stdout, 2000) na, totcharge(1), &
                ( charges(na,l,is), l= 0,lmax_wfc)
        else if ( nspin == 2) then 
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( charges(na,l,1), l= 0,lmax_wfc) 
           WRITE( stdout, 2002) totcharge(2), &
                ( charges(na,l,2), l= 0,lmax_wfc) 
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                 ( charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        end if
     end do 
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4 ,&
          & ", s, p, d, f = ",4f8.4) 
2001 FORMAT (15x,"  spin up      = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
2002 FORMAT (15x,"  spin down    = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
2003 FORMAT (15x,"  polarization = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
     !
     psum = SUM(charges(:,:,:)) / nelec 
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0 - psum 
     ! 
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995). 
     ! The spilling parameter measures the ability of the basis provided by 
     ! the pseudo-atomic wfc to represent the PW eigenstates, 
     ! by measuring how much of the subspace of the Hamiltonian 
     ! eigenstates falls outside the subspace spanned by the atomic basis 
     ! 
     deallocate (charges) 
     !
#ifdef __PARA 
  endif
#endif 
  return
end subroutine projwave

!----------------------------------------------------------------------- 
subroutine  partialdos (Emin, Emax, DeltaE, filpdos)
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout 
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev 
  !
  USE projections
  !
  implicit none 
  character (len=80) :: filpdos
  real(kind=DP) :: Emin, Emax, DeltaE
  !
  character (len=33) :: filextension 
  character (len=256):: fileout
  character (len=1)  :: l_label(0:3)=(/'s','p','d','f'/) 
  ! 
  integer :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, lm, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is 
  real(kind=DP) :: etev, delta, Elw, Eup 
  real(kind=DP), allocatable :: dostot(:,:), pdos(:,:,:), pdostot(:,:), &
       ldos(:,:)
  real(kind=DP), external :: w0gauss 
  ! 
  ! 
  ! find band extrema 
  ! 
  Elw = et (1, 1) 
  Eup = et (nbnd, 1) 
  do ik = 2, nkstot 
     Elw = min (Elw, et (1, ik) ) 
     Eup = max (Eup, et (nbnd, ik) ) 
  enddo
  if (degauss.ne.0.d0) then
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  endif
  Emin = max (Emin/rytoev, Elw) 
  Emax = min (Emax/rytoev, Eup) 
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001)
  !
  allocate (pdos(0:ne,natomwfc,nspin))
  allocate (dostot(0:ne,nspin), pdostot(0:ne,nspin), ldos(0:ne,nspin) )
  pdos(:,:,:) = 0.d0 
  dostot(:,:) = 0.d0 
  pdostot(:,:)= 0.d0 
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1 

  do ik = 1,nkstot 
     if ( nspin == 2 ) current_spin = isk ( ik ) 
     do ibnd = 1, nbnd 
        etev = et(ibnd,ik)
        ie_mid = nint( (etev-Emin)/DeltaE ) 
        do ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne) 
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           ! pdos(:,nwfc,ns) = DOS (states/eV) for spin "ns" 
           !                   projected over atomic wfc "nwfc"
           !
           do nwfc = 1, natomwfc 
              pdos(ie,nwfc,current_spin) = pdos(ie,nwfc,current_spin) + & 
                   wk(ik) * delta * proj (nwfc, ibnd, ik)  
           end do
           !
           ! dostot(:,ns) = total DOS (states/eV) for spin "ns" 
           !
           dostot(ie,current_spin) = dostot(ie,current_spin) + & 
                wk(ik) * delta 
        end do
     end do
  end do
  !
  ! pdostot(:,ns) = sum of all projected DOS
  !
  do is=1,nspin 
     do ie=0,ne 
        pdostot(ie,is) = sum(pdos(ie,:,is)) 
     end do
  end do
  
  do nwfc = 1, natomwfc 
     if (nlmchi(nwfc)%m == 1) then 
        filextension='.pdos_atm#' 
        !             12345678901 
        c_tab = 11 
        if (nlmchi(nwfc)%na < 10) then 
           write (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na 
           c_tab = c_tab + 1 
        else if (nlmchi(nwfc)%na < 100) then 
           write (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na 
           c_tab = c_tab + 2 
        else if (nlmchi(nwfc)%na < 1000) then 
           write (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na 
           c_tab = c_tab + 3 
        else 
           call errore('projwave',& 
                'file extension not supporting so many atoms', nwfc) 
        endif
        write (filextension(c_tab:c_tab+4),'(a1,a)') & 
             '(',trim(atm(ityp(nlmchi(nwfc)%na))) 
        c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1 
        if (nlmchi(nwfc)%n >= 10) & 
             call errore('projwave',& 
             'file extension not supporting so many atomic wfc', nwfc) 
        if (nlmchi(nwfc)%l > 3) & 
             call errore('projwave',& 
             'file extension not supporting so many l', nwfc) 
        write (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  & 
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l) 
        fileout = trim(filpdos)//trim(filextension)
        open (4,file=fileout,form='formatted', & 
             status='unknown') 

        if (nspin == 1) then 
           write (4,'("# E (eV)   ldos(E)  ",$)') 
        else 
           write (4,'("# E (eV)  ldosup(E)  ldosdw(E)",$)') 
        end if
        do m=1,2 * nlmchi(nwfc)%l + 1 
           if (nspin == 1) then 
              write(4,'(" pdos(E)   ",$)')  
           else 
              write(4,'(" pdosup(E) ",$)')  
              write(4,'(" pdosdw(E) ",$)')  
           end if
        end do
        write(4,*) 
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:) = 0.d0
        do ie= 0, ne 
           do is=1, nspin
              do m=1,2 * nlmchi(nwfc)%l + 1 
                 ldos  (ie, is) = ldos  (ie, is) + pdos(ie,nwfc+m-1,is)
              end do
           end do
        end do
        do ie= 0, ne 
           etev = Emin + ie * DeltaE 
           write (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  & 
                (ldos(ie,is), is=1,nspin), & 
                ((pdos(ie,nwfc+m-1,is), is=1,nspin), & 
                m=1,2*nlmchi(nwfc)%l+1)
        end do
        close (4) 
     end if
  end do
  fileout = trim(filpdos)//".pdos_tot"
  open (4,file=fileout,form='formatted', status='unknown') 
  if (nspin == 1) then 
     write (4,'("# E (eV)  dos(E)    pdos(E)")')  
  else 
     write (4,'("# E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)")')  
  end if
  do ie= 0, ne 
     etev = Emin + ie * DeltaE 
     write (4,'(f7.3,4e11.3)') etev*rytoev, (dostot(ie,is), is=1,nspin), & 
          (pdostot(ie,is), is=1,nspin) 
  end do
  close (4) 
  deallocate (ldos, dostot, pdostot)
  deallocate (pdos)
  !
  deallocate (nlmchi) 
  deallocate (proj) 
  !
  return 
end subroutine partialdos
