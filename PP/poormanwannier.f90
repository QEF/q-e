! Copyright (C) 2001-2003 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!----------------------------------------------------------------------- 
program poormanwannier 
  !----------------------------------------------------------------------- 
  ! 
  ! projects wavefunctions onto atomic wavefunctions, 
  ! 
  ! input: namelist "&inputpp", with variables 
  !   prefix      prefix of input files saved by program pwscf 
  !   outdir      temporary directory where files resides 
  ! 
  USE io_global,  ONLY : stdout 
  USE kinds, only : DP 
  use io_files,   only : nd_nmbr, prefix, tmp_dir 
#ifdef __PARA 
  use para,       only : me 
  use mp 
#endif 
  implicit none 
  character(len=256) :: outdir 
  integer :: ios, ionode_id = 0  
  integer :: first_band, last_band
  namelist / inputpp / outdir, prefix, first_band, last_band
  ! 
  call start_postproc (nd_nmbr) 
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf' 
  outdir = './' 
  first_band=-1
  last_band=-1
  ! 
#ifdef __PARA 
  if (me == 1)  then 
#endif 
  read (5, inputpp, err = 200, iostat = ios) 
200 call errore ('pmwannier', 'reading inputpp namelist', abs (ios) ) 
  ! 
  tmp_dir = trim(outdir) 
  ! 
#ifdef __PARA 
  end if 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( first_band, ionode_id ) 
  CALL mp_bcast( last_band, ionode_id ) 
#endif 
  ! 
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  call read_file 
  call openfil 
  !
  call projection( first_band, last_band)
  ! 
  call stop_pp 
  ! 
end program poormanwannier
 
!----------------------------------------------------------------------- 
subroutine projection (first_band, last_band)
  !----------------------------------------------------------------------- 
  ! 
#include "machine.h" 
#define ONE  (1.D0,0.D0)
#define ZERO (0.D0,0.D0)  

  USE io_global,  ONLY : stdout 
  use atom 
  use basis 
  USE cell_base
  use constants, only: rytoev 
  use gvect 
  use klist 
  use ldaU 
  use lsda_mod 
  use symme, only: nsym, irt 
  use wvfct 
  use us 
  use becmod,   only: becp 
  use rbecmod,  only: rbecp => becp 
  use io_files, only: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, &
                      iunat, nwordatwfc, iunigk
  use wavefunctions_module, only: evc 
#ifdef __PARA 
  use para 
#endif 
  implicit none 
  !
  ! I/O variables 
  !
  integer :: first_band, last_band
  !
  ! local variables
  !
  integer :: ik, ia, ib, na, nt, n, m, l, nwfc, lmax_wfc, &
             ldim1, ldim2, lwork, i, j, info, counter, counter_ldau
  logical :: exst 
  complex(kind=DP), allocatable :: proj (:,:,:)
  complex(kind=DP), allocatable :: wfcatom (:,:) 
  integer, allocatable :: index(:) 
  ! 
  complex(kind=DP), allocatable ::  proj0(:,:) 
  ! Some workspace for k-point calculation ... 
  real   (kind=DP), allocatable :: rproj0(:,:) 
  ! ... or for gamma-point. 
  COMPLEX(KIND=DP), ALLOCATABLE :: pp(:,:), u_m(:,:), w_m(:,:), work(:)
  ! the overlap matrix pp
  ! left unitary matrix in the SVD of sp_m
  ! right unitary matrix in the SVD of sp_m
  ! workspace for ZGESVD
  REAL(KIND=DP), ALLOCATABLE :: ew(:), rwork(:)
  ! the eigenvalues of pp
  ! workspace for ZGESVD
  REAL (kind=DP) :: capel
  ! 
  WRITE( stdout, '(/5x,"Calling projection .... ")') 
  if ( gamma_only ) WRITE( stdout, '(5x,"gamma-point specific algorithms are used")') 
  ! 
  allocate(proj (natomwfc, nbnd, nkstot) ) 
  allocate(wfcatom (npwx, natomwfc) ) 
  if ( gamma_only ) then 
     ! Allocate the array rbecp=<beta|wfcatom>, which is not allocated by 
     ! allocate_wfc() in a gamma-point calculation: 
     allocate (rbecp (nkb,natomwfc)) 
  end if 

  if (first_band == -1)  first_band = 1
  if (last_band  == -1)  last_band  = nbnd
  if (first_band > last_band ) call errore ('pmw',' first_band > last_band',1)
  if (first_band < 0         ) call errore ('pmw',' first_band < 0 ',       1)
  if (last_band > nbnd       ) call errore ('pmw',' last_band > nbnd ',     1)
 

  counter = 0
  counter_ldaU = 0
  do na = 1, nat
     nt = ityp (na)
     do n = 1, nchi (nt)
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then
           l = lchi (n, nt)
           if ( (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) .and. &
                                            l.eq.Hubbard_l(nt) )then
               counter_ldaU = counter_ldaU + 2 * l + 1
           end if
           counter = counter + 2 * l + 1
        endif
     enddo
  enddo

  WRITE( stdout, *) "    NBND = ", nbnd
  WRITE( stdout, *) "    NATOMWFC =", natomwfc
  WRITE( stdout, *) "    NKSTOT =", nkstot
  
  ldim1 = counter_ldaU
  ldim2 = last_band + 1 - first_band
  WRITE( stdout, *) ldim1, ldim2

  if (ldim1 > ldim2 ) call errore( 'projection','too few bands',ldim1-ldim2)
  lwork = 5 * max(ldim1,ldim2)
  allocate (pp(ldim1,ldim2), u_m(ldim1,ldim1), w_m(ldim2,ldim2), &
            work(lwork), ew(ldim1), rwork(lwork))
  proj   = 0.d0 
  ! 
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1 
  ! 
  call d_matrix (d1, d2, d3)   
  write (stdout,*) " Hubbard_lmax = ", Hubbard_lmax, lda_plus_u
  nwfc=0 
  lmax_wfc = 0 
  do na = 1, nat 
     nt = ityp (na) 
     do n = 1, nchi (nt) 
        if (oc (n, nt) > 0.d0 .OR..NOT.newpseudo (nt) ) then 
           l = lchi (n, nt) 
           lmax_wfc = max (lmax_wfc, l ) 
           do m = 1, 2 * l + 1 
              nwfc=nwfc+1 
              write(stdout,*) " ATOMIC WFC #", nwfc,":", na,n,l,m
           enddo 
        endif 
     enddo 
  enddo 
 
  if (lmax_wfc > 3) call errore ('projection', 'l > 3 not yet implemented', 1) 
  if (nwfc /= natomwfc) call errore ('projection', 'wrong # of atomic wfcs?', 1)
  ! 
  !    loop on k points 
  ! 
  call init_us_1 
  call init_at_1 
  ! 
  do ik = 1, nks 
     WRITE ( stdout, * ) "KPOINT =", ik
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
 
     ! 
     ! make the projection <psi_i| \hat S | phi_j> 
     ! 
     if ( gamma_only ) then 
        allocate(rproj0(natomwfc,nbnd) ) 
        call pw_gemm ('Y', natomwfc, nbnd, npw, swfcatom, npwx, evc, npwx, &
                       rproj0, natomwfc) 
        proj(:,:,ik) = cmplx(rproj0(:,:),0.d0)
        deallocate (rproj0) 
     else 
        allocate(proj0(natomwfc,nbnd) ) 
        call ccalbec (natomwfc, npwx, npw, nbnd, proj0, swfcatom, evc) 
        proj(:,:,ik) = proj0(:,:)
        deallocate (proj0) 
     end if 

     counter = 0
     counter_ldaU = 0
     do na = 1, nat
        nt = ityp (na)
        do n = 1, nchi (nt)
           if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then
              l = lchi (n, nt)
              if ( (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) .and. &
                                            l.eq.Hubbard_l(nt) )then
                  pp(counter_ldaU+1:counter_ldaU+2*l+1, 1:ldim2) = &
                      proj(counter+1:counter+2*l+1,first_band:last_band,ik)
                  counter_ldaU = counter_ldaU + 2 * l + 1
              end if
              counter = counter + 2 * l + 1
           endif
        enddo
     enddo
     if (counter_ldaU .ne. ldim1) call errore ('projection','wrong counter',1)

     CALL ZGESVD( 'A', 'A', ldim1, ldim2, pp, ldim1, ew, u_m, ldim1, &
                  w_m, ldim2, work, lwork, rwork, info )
     call errore ('projection','Singular Value Deconposition failed', abs(info))
     do i = 1, ldim1
        WRITE ( stdout, * ) ew(i)
        WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        WRITE ( stdout, '(8(2f5.2,2x))') w_m(i,:)
     end do
     !
     ! ... use sp_m to store u_m * w_m
     !
     CALL ZGEMM( 'N', 'N', ldim1, ldim2, ldim1, ONE, u_m, ldim1, w_m, &
                    ldim2, ZERO, pp, ldim1 )
     ! ... check orthogonality
     CALL ZGEMM( 'N', 'C', ldim1, ldim1, ldim2, ONE, pp, ldim1, pp, &
                    ldim1, ZERO, u_m, ldim1 )
     capel = 0.d0
     do i=1,ldim1
        u_m(i,i) = u_m(i,i) -1.d0
        do j=1,ldim1
           capel = capel + abs( u_m(i,j) )
        end do
        u_m(i,i) = u_m(i,i) +1.d0
     end do

     if (capel < 1.d-10) then
        WRITE ( stdout, *) " ORTHOGONALITY CHECK PASSED "
     else
        WRITE ( stdout, *) " ORTHOGONALITY CHECK FAILED"
        WRITE ( stdout, *) " CAPEL = ", capel
        do i=1,ldim1
           WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        end do
     end if
     counter = 0
     counter_ldaU = 0
     do na = 1, nat
        nt = ityp (na)
        do n = 1, nchi (nt)
           if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then
              l = lchi (n, nt)
              if ( (Hubbard_U(nt).ne.0.d0.or.Hubbard_alpha(nt).ne.0.d0) .and. &
                                            l.eq.Hubbard_l(nt) )then
                  CALL ZGEMM( 'N', 'C', npw, 2*l+1, ldim2, ONE, &
                              evc(1,first_band), npwx, &
                              pp(counter_ldaU+1,1), ldim1, ZERO, &
                              wfcatom(1,counter+1), npwx )
                  counter_ldaU = counter_ldaU + 2 * l + 1
              end if
              counter = counter + 2 * l + 1
           endif
        enddo
     enddo
     if ( gamma_only ) then 
        call pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, wfcatom, npwx, rbecp, nkb)   
     else 
        call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom) 
     end if 
     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 

     call davcio (swfcatom, nwordatwfc, iunat, ik, 1)


     ! on k-points 
  enddo 
  ! 
  if ( gamma_only ) then 
     deallocate (rbecp) 
  end if 
  ! 

  deallocate (wfcatom) 
 
  deallocate (proj) 
  return 
end subroutine projection
