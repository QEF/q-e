!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
module restart_module

  USE io_global,  ONLY : stdout

  implicit none
  save

contains

!----------------------------------------------------------------------
subroutine writefile_new( what, ndw, et_g, wg_g, kunit )
  !-----------------------------------------------------------------------
  !
  !     This routine is called at the end of the run to save on a file
  !     the information needed to restart and to other postprocessing
  !     programs.
  !
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : ntyp => nsp
  USE basis,                ONLY : natomwfc
  USE ions_base,            ONLY : nat, ityp, tau, zv, atm
  USE cell_base,            ONLY : at, bg, ibrav, celldm, alat, symm_type
  USE klist,                ONLY : xk, wk, degauss, ngauss, lgauss, nelec, &
                                   ngk, nks, nkstot
  USE ktetra,               ONLY : tetra, ntetra, ltetra, k1, k2, k3, &
                                   nk1, nk2, nk3
  USE lsda_mod,             ONLY : isk, nspin, lsda
  USE force_mod,            ONLY : lforce, lstres, force
  USE gvect,                ONLY : dual, ecutwfc, nrx1, nrx2, nrx3, nr1, nr2, &
                                   nr3, nrxx, ngm, ngm_g, gcutm
  USE gsmooth,              ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
                                   nrxxs, gcutms, doublegrid
  USE wvfct,                ONLY : npwx, nbndx, nbnd, igk, g2kin, &
                                   igk_l2g, gamma_only
  USE char,                 ONLY : title, crystal, sname
  USE dynam,                ONLY : amass
  USE symme,                ONLY : s, irt, ftau, nsym, invsym
  USE ener,                 ONLY : ef
  USE atom,                 ONLY : zmesh, xmin, dx, r, rab, chi, oc, rho_at, &
                                   rho_atc, mesh, msh, nchi, lchi, jchi, &
                                   numeric, nlcc
  USE pseud,                ONLY : cc, alpc, zp, aps, alps, nlc, nnl, lmax, &
                                   lloc, a_nlcc, b_nlcc, alpha_nlcc
  USE us,                   ONLY : okvan
  USE uspp_param,           ONLY : vloc_at, dion, betar, qqq, qfunc, qfcoef, &
                                   rinner, psd, nbeta, kkbeta, nqf, nqlc, &
                                   ifqopt, lll, jjj, iver, nh, tvanp, newpseudo
  USE extfield,             ONLY : tefield, dipfield, edir, emaxpos, eopreg, &
                                   eamp
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE fixed_occ,            ONLY : tfixed_occ 
  USE control_flags,        ONLY : twfcollect, noinv, istep, iswitch, modenum
  USE io_files,             ONLY : prefix, tmp_dir, pseudo_dir, psfile, &
                                   iunwfc, nwordwfc
  USE funct,                ONLY : iexch, icorr, igcx, igcc
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_sum, mp_max, mp_end
  USE mp_global,            ONLY : mpime, nproc, root, me_pool, my_pool_id, &
                                   nproc_pool, intra_pool_comm, root_pool, &
                                   inter_pool_comm, my_image_id
  USE io_base,              ONLY : write_restart_header, write_restart_ions, &
                                   write_restart_cell, write_restart_electrons, &
                                   write_restart_gvec, write_restart_gkvec, &
                                   write_restart_charge, write_restart_wfc, &
                                   write_restart_symmetry, write_restart_xdim, &
                                   write_restart_pseudo, write_restart_ldaU, &
                                   write_restart_tetra
  USE parameters,           ONLY : nacx, nsx, npk
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_alpha
  USE read_pseudo_module
  USE pseudo_types
  USE spin_orb, ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin, npol

  !
  implicit none
  !
  integer, intent(in) :: ndw
  character(len=*), intent(in) :: what
  real(kind=DP), intent(in) :: et_g(:,:), wg_g(:,:)
  integer, intent(in) :: kunit
  !
  integer :: ik, i, ibnd, ia, ispin, ldim, npwt
  logical :: exst
  logical :: twrite
  real(kind=DP) :: trutime
  integer :: nelu
  integer :: neld
  integer :: na(nsx)
  integer :: ngk_g( npk )
  integer :: ngk_l( npk )
  real(kind=DP) :: acc(nacx)
  real(kind=DP) :: ekincm, ecutrho
  real(kind=DP), allocatable :: stau0(:,:), staum(:,:)
  real(kind=DP), allocatable :: svel0(:,:), svelm(:,:), tautmp(:,:)
  real(kind=DP) :: xnosp, xnos0, xnosm, xnosm2, cdmi(3)
  character(len=4) :: atom_label(nsx)
  real(kind=DP), dimension(3,3) :: ht0, htm, htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2
  logical :: tscal, tocc, tlam, teig, tmill
  real(kind=DP) :: xenos0, xenosm, xenosm2, xenosp
  real(kind=DP) :: bg1_(3), bg2_(3), bg3_(3)
  real(kind=DP), allocatable :: occtmp(:), lambda(:,:), g_g(:)
  integer :: strlen
  integer, allocatable :: mill(:,:)
  integer :: ngm_p( nproc_pool )
  character(len=256) :: filename, file_pseudo
  LOGICAL :: twf0, twfm, tupf
  INTEGER :: nkl, nkr, nkbl, iks, ike, npw_g, ipol, j
  INTEGER :: npool, ipmask( nproc ), ipsour
  real(kind=DP) :: wfc_scal

  LOGICAL :: twrhead, twrxdim, twrcell, twrpos, twrpseudo, twrchden
  LOGICAL :: twrocc, twrgvec, twrgkvec, twrwfc, twrsym

  INTEGER :: iunps, ios, flen, ierr
  LOGICAL :: opnd
  LOGICAL :: lgamma
  COMPLEX(KIND=DP), ALLOCATABLE :: wfc_restart(:,:)
  TYPE (pseudo_upf) :: upf

  external DSCAL

! ... end of declarations

  !
  !
  filename = trim(prefix)//'.save'
  ! WRITE( stdout, '(/,5x,"Writing file ",a14)') filename
  !
  if( ionode ) THEN
    call seqopn (ndw, filename, 'unformatted', exst)
    rewind ndw
  end if

  twrhead   = .FALSE.
  twrxdim   = .FALSE.
  twrcell   = .FALSE.
  twrpos    = .FALSE.
  twrsym    = .FALSE.
  twrpseudo = .FALSE.
  twrocc    = .FALSE.
  twrgvec   = .FALSE.
  twrgkvec  = .FALSE.
  twrchden  = .FALSE.
  twrwfc    = .FALSE.

  SELECT CASE ( TRIM( what ) )
     CASE ( 'all' )
      twrhead   = .TRUE.
      twrxdim   = .TRUE.
      twrcell   = .TRUE.
      twrpos    = .TRUE.
      twrsym    = .TRUE.
      twrpseudo = .TRUE.
      twrocc    = .TRUE.
      twrgvec   = .TRUE.
      twrgkvec  = .TRUE.
      twrchden  = .TRUE.
      if( twfcollect ) twrwfc    = .TRUE.
    CASE ( 'config' )
      twrhead   = .TRUE.
      twrxdim   = .TRUE.
      twrcell   = .TRUE.
      twrpos    = .TRUE.
    CASE DEFAULT
      CALL errore( ' writefile_new ', ' unknown value for what ', 1 )
  END SELECT



!  ==--------------------------------------------------------------==
!  ==  WRITE HEADER INFORMATION                                    ==
!  ==--------------------------------------------------------------==
   trutime = 0.0d0
   nelu = 0
   neld = 0
   na = 0.0d0
   acc = 0.0d0
   ekincm = 0.0d0
   ecutrho = dual * ecutwfc

   do i = 1, nat
     na(ityp(i)) = na(ityp(i)) + 1
   end do

   ngk_g = 0
   ngk_l = 0

   IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .OR. ( MOD( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' writefile_new ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .OR. ( MOD( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' writefile_new ',' nproc_pool ', 1 )

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

     ngk_g(iks:ike) = ngk(1:nkl)
     CALL mp_sum( ngk_g, intra_pool_comm )

     !write(400+mpime,*) what, ngk(1:nks), nkstot
     !write(400+mpime,*) what, ngk_g(1:nkstot)

     IF( npool > 1 ) THEN
       CALL mp_sum( ngk_g, inter_pool_comm )
       !write(400+mpime,*) what, ngk_g(1:nkstot)
     END IF

  END IF

   tupf = .FALSE.  ! the pseudopotential are not saved in UPF format
   ! tupf = .TRUE.  ! the pseudopotential are saved in UPF format
   lgamma = gamma_only

   if( twrhead ) then
     CALL write_restart_header(ndw, istep, iswitch, trutime, nr1, nr2, nr3, &
       nr1s, nr2s, nr3s, ngm_g, nkstot, ngk_g, nspin, nbnd, nelec, nelu, neld, &
       nat, ntyp, na, acc, nacx, ecutwfc, ecutrho, alat, ekincm, &
       kunit, k1, k2, k3, nk1, nk2, nk3, degauss, ngauss, lgauss, ntetra, ltetra, &
       natomwfc, gcutm, gcutms, dual, doublegrid, modenum, lforce, lstres, &
       title, crystal, tmp_dir, tupf, lgamma, noncolin, lspinorb, lda_plus_u, &
       tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect)

   else
     CALL write_restart_header(ndw)
   end if

!  ==--------------------------------------------------------------==
!  ==  MAX DIMENSIONS                                              ==
!  ==--------------------------------------------------------------==

   if( twrxdim ) then
     CALL write_restart_xdim( ndw, &
       npwx, nbndx, nrx1, nrx2, nrx3, nrxx, nrx1s, nrx2s, nrx3s, nrxxs )
   else
     CALL write_restart_xdim( ndw )
   end if

!  ==--------------------------------------------------------------==
!  ==  CELL & METRIC                                               ==
!  ==--------------------------------------------------------------==

   xhnosm2 = 0.0d0
   xhnosm  = 0.0d0
   xhnos0  = 0.0d0
   xhnosp  = 0.0d0
   ht0 = TRANSPOSE(at) * alat
   htm = ht0
   htm2 = ht0
   htvel = 0.0d0
   if( twrcell ) then
     CALL write_restart_cell( ndw, ibrav, celldm, ht0, htm, &
          htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2)
   else
     CALL write_restart_cell( ndw )
   end if

!  ==--------------------------------------------------------------==
!  ==  IONS                                                        ==
!  ==--------------------------------------------------------------==

   ALLOCATE( stau0(3, nat) )
   ALLOCATE( staum(3, nat) )
   ALLOCATE( svel0(3, nat) )
   ALLOCATE( svelm(3, nat) )
   ALLOCATE( tautmp(3, nat) )

   DO ia = 1, nat
     stau0(:,ia) = alat * tau(:,ia)
     staum(:,ia) = alat * tau(:,ia)
     svel0(:,ia) = 0.0d0
     svelm(:,ia) = 0.0d0
     tautmp(:,ia) = alat * tau(:,ia)
   END DO
   xnosp = 0.0d0
   xnos0 = 0.0d0
   xnosm = 0.0d0
   xnosm2 = 0.0d0
   atom_label(1:ntyp) = ' '
   cdmi = 0.0d0
   tscal = .FALSE.

   if( twrpos ) then
     CALL write_restart_ions(ndw, atm, tscal, stau0, svel0, &
       staum, svelm, tautmp, force, cdmi, nat, ntyp, ityp, na, amass, xnosp,  &
       xnos0, xnosm, xnosm2)
   else
     CALL write_restart_ions(ndw)
   end if

   DEALLOCATE( stau0, svel0, staum, svelm, tautmp )

!  ==--------------------------------------------------------------==
!  ==  LDA+U                                                       ==
!  ==--------------------------------------------------------------==
     if (lda_plus_u) then
        CALL write_restart_ldaU(ndw, ntyp, Hubbard_lmax, &
             Hubbard_l(1:ntyp), Hubbard_U(1:ntyp), Hubbard_alpha(1:ntyp))
     else
        CALL write_restart_ldaU(ndw)
     endif
  
!  ==--------------------------------------------------------------==
!  ==  SYMMETRIES                                                  ==
!  ==--------------------------------------------------------------==

   IF( twrsym ) THEN
     CALL write_restart_symmetry( ndw, &
        symm_type, sname, s, irt, nat, ftau, nsym, invsym, noinv )
   ELSE
     CALL write_restart_symmetry( ndw )
   END IF

!  ==--------------------------------------------------------------==
!  ==  PSEUDOPOTENTIALS                                            ==
!  ==--------------------------------------------------------------==

   DO i = 1, ntyp
     IF( twrpseudo ) THEN
       IF( tupf ) THEN
         iunps = 10
         flen = len_trim (pseudo_dir)
         if (pseudo_dir (flen:flen) .ne.'/') then
           file_pseudo = pseudo_dir (1:flen) //'/'//psfile (i)
         else
           file_pseudo = pseudo_dir (1:flen) //psfile (i)
         endif
         inquire (unit = iunps, opened = opnd)
         IF( opnd ) &
           CALL errore( ' writefile_new ', ' fortran unit already opened ', 1 )
         open (unit = iunps, file = file_pseudo, status = 'old', form = &
           'formatted', iostat = ios)

         call read_pseudo_upf(iunps, upf, ierr)

         CALL write_restart_pseudo( ndw, &
           upf%generated, upf%date_author, upf%comment, upf%psd, upf%typ, &
           upf%tvanp,  &
           upf%nlcc, upf%dft, upf%zp, upf%etotps, upf%ecutwfc, upf%ecutrho, &
           upf%nv,   &
           upf%lmax, upf%mesh, upf%nwfc, upf%nbeta, upf%els(:), upf%lchi(:), &
           upf%jchi(:), upf%oc(:), &
           upf%r(:), upf%rab(:), upf%rho_atc(:), upf%vloc(:), upf%lll(:), &
           upf%jjj(:), upf%kkbeta(:), &
           upf%beta(:,:), upf%nd, upf%dion(:,:), upf%nqf, upf%nqlc, &
           upf%rinner(:), &
           upf%qqq(:,:), upf%qfunc(:,:,:), upf%qfcoef(:,:,:,:), &
           upf%chi(:,:), upf%rho_at(:) )

         CALL deallocate_pseudo_upf( upf )
         close( iunps )
       ELSE
         CALL write_restart_pseudo( ndw, &
           zmesh(i), xmin(i), dx(i), r(:,i), rab(:,i), vloc_at(:,i), &
           chi(:,:,i), oc(:,i), &
           rho_at(:,i), rho_atc(:,i), mesh(i), msh(i), nchi(i), lchi(:,i), &
           jchi(:,i), &
           numeric(i), cc(:,i), alpc(:,i), zp(i), aps(:,:,i), alps(:,:,i), &
           zv(i), nlc(i), nnl(i), lmax(i), lloc(i), dion(:,:,i), &
           betar(:,:,i), qqq(:,:,i), qfunc(:,:,:,i), qfcoef(:,:,:,:,i), &
           rinner(:,i), nh(i), nbeta(i), kkbeta(i), nqf(i), nqlc(i), &
           ifqopt(i), &
           lll(:,i), jjj(:,i), iver(:,i), tvanp(i), okvan, &
           newpseudo(i), iexch, icorr, &
           igcx, igcc, lsda, a_nlcc(i), b_nlcc(i), alpha_nlcc(i), &
           nlcc(i), psd(i) )
       END IF
     ELSE
       CALL write_restart_pseudo( ndw )
     END IF
   END DO

!  ==--------------------------------------------------------------==
!  ==  OCCUPATION NUMBER                                           ==
!  ==--------------------------------------------------------------==

   xenos0 = 0.0d0
   xenosm = 0.0d0
   xenosm2 = 0.0d0
   xenosp = 0.0d0
   tocc = .FALSE.
   tlam = .FALSE.
   teig = .TRUE.
   ldim = 1

   DO ik = 1, nkstot
     IF( twrocc ) THEN
       ALLOCATE( occtmp(nbnd) )
       ALLOCATE( lambda(1,1) )
       occtmp = 0.0d0
       ispin = isk( ik )
       CALL write_restart_electrons( ndw, occtmp, occtmp, tocc, lambda, lambda, &
         ldim, tlam, nbnd, ispin, nspin, ik, nkstot, nelec, nelu, neld, &
         xenosp, xenos0, xenosm, xenosm2, &
         ef, teig, et_g(:,ik), wg_g(:,ik) )
       DEALLOCATE( occtmp )
       DEALLOCATE( lambda )
     ELSE
       CALL write_restart_electrons( ndw )
     END IF
   END DO

!  ==--------------------------------------------------------------==
!  ==  G-Vectors                                                   ==
!  ==--------------------------------------------------------------==

   IF( twrgvec ) THEN
     ALLOCATE( mill(3,1) )
     mill = 0
     tmill = .FALSE.
     bg1_ = bg(:,1) 
     bg2_ = bg(:,2)
     bg3_ = bg(:,3)
     CALL write_restart_gvec( ndw, ngm_g, bg1_, bg2_, bg3_, &
       bg1_, bg2_, bg3_, tmill, mill )
     DEALLOCATE( mill )
   ELSE
     CALL write_restart_gvec( ndw )
   END IF

!  ==--------------------------------------------------------------==
!  ==  (G+k)-Vectors                                               ==
!  ==--------------------------------------------------------------==

   DO ik = 1, nkstot
     IF( twrgkvec ) THEN
       npwt = npwx
       CALL write_restart_gkvec(ndw, ik, nkstot, ngk_g(ik), &
         xk(:,ik), wk(ik), isk(ik))
     ELSE
       CALL write_restart_gkvec( ndw )
     END IF
   END DO

!  ==--------------------------------------------------------------==
!  ==  Tetrahedra                                                  ==
!  ==--------------------------------------------------------------==

   CALL write_restart_tetra( ndw, ltetra, ntetra, tetra )

!  ==--------------------------------------------------------------==
!  ==  CHARGE DENSITY AND POTENTIALS                               ==
!  ==--------------------------------------------------------------==

   DO ispin = 1, nspin
     CALL write_restart_charge( ndw )
   END DO

!  ==--------------------------------------------------------------==
!  ==  WAVEFUNCTIONS                                               ==
!  ==--------------------------------------------------------------==

   IF( twrwfc ) THEN

     twf0   = .TRUE.
     twfm   = .FALSE.
     wfc_scal = 1.0d0
     npw_g = MAXVAL( igk_l2g(:,:) )
     CALL mp_max( npw_g )

     DO ik = 1, nkstot

       IF( (ik >= iks) .AND. (ik <= ike) ) THEN
         if (noncolin) then
            call davcio (evc_nc, nwordwfc, iunwfc, (ik-iks+1), - 1)
         else
            call davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)
         endif
       END IF
       ispin = isk( ik )
       ! WRITE( stdout,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool ! DEBUG
       if (noncolin) then
         allocate (wfc_restart(npwx, nbnd ))
         wfc_restart=(0.d0,0.d0)
         do ipol = 1, npol
            wfc_restart = evc_nc(:,ipol,:)
            CALL write_restart_wfc(ndw, ik, nkstot, kunit, ispin, nspin, &
            wfc_scal, wfc_restart, twf0, wfc_restart, twfm, npw_g, nbnd, &
            igk_l2g(:,ik-iks+1), ngk(ik-iks+1) )
         enddo
         deallocate(wfc_restart)
       else
         CALL write_restart_wfc(ndw, ik, nkstot, kunit, ispin, nspin, &
            wfc_scal, evc, twf0, evc, twfm, npw_g, nbnd, igk_l2g(:,ik-iks+1), &
            ngk(ik-iks+1) )
       endif

     ENDDO

   ELSE

     DO ik = 1, nkstot
       CALL write_restart_wfc(ndw, nbnd )
     END DO

   END IF

  ! WRITE( stdout, '(5x,"file written")')

  if( ionode ) then
    close (unit = ndw)
  end if
  !
  return
end subroutine writefile_new


!----------------------------------------------------------------------
subroutine readfile_new( what, ndr, et_g, wg_g, kunit, nsizwfc, iunitwfc, ierr )
  !-----------------------------------------------------------------------
  !
  !     This routine is called at the start of the run to read from file
  !     the information needed to restart and to other postprocessing
  !     programs.
  !
  !
  USE parameters,           ONLY : npk, nchix, ndmx, nbrx, lqmax, nqfx
  USE constants,            ONLY : pi
  USE io_files,             ONLY : iunwfc, nwordwfc, prefix, tmp_dir
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE basis,                ONLY : natomwfc
  USE cell_base,            ONLY : at, bg, ibrav, celldm, alat, tpiba, tpiba2, &
                                   omega, symm_type
  USE klist,                ONLY : xk, wk, degauss, ngauss, lgauss, nelec, &
                                   ngk, nks, nkstot
  USE ktetra,               ONLY : tetra, ntetra, ltetra, k1, k2, k3, nk1, nk2, nk3
  USE lsda_mod,             ONLY : isk, nspin, lsda
  USE force_mod,            ONLY : lforce, lstres, force
  USE gvect,                ONLY : dual, ecutwfc, nrx1, nrx2, nrx3, nr1, nr2, &
                                   nr3, nrxx, ngm, ngm_g, gcutm, g, ig_l2g
  USE gsmooth,              ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
                                   nrxxs, gcutms, doublegrid
  USE wvfct,                ONLY :  npwx, nbndx, nbnd, igk, g2kin, igk_l2g, gamma_only
  USE char,                 ONLY : title, crystal, sname
  USE dynam,                ONLY : amass
  USE symme,                ONLY : s, irt, ftau, nsym, invsym
  USE ener,                 ONLY : ef
  USE atom,                 ONLY : zmesh, xmin, dx, r, rab, chi, oc, rho_at, &
                                   rho_atc, mesh, msh, nchi, lchi, jchi, &
                                   numeric, &
                                   nlcc
  USE pseud,                ONLY : cc, alpc, zp, aps, alps, nlc, nnl, lmax, &
                                   lloc, a_nlcc, b_nlcc, alpha_nlcc
  USE us,                   ONLY : okvan
  USE uspp_param,           ONLY : vloc_at, dion, betar, qqq, qfunc, qfcoef, &
                                   rinner, psd, nbeta, kkbeta, nqf, nqlc, &
                                   ifqopt, lll, jjj, iver, nh, tvanp, newpseudo
  USE extfield,             ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE fixed_occ,            ONLY : tfixed_occ 
  USE control_flags,        ONLY : twfcollect, noinv, istep, iswitch, modenum
  USE funct,                ONLY : iexch, icorr, igcx, igcc
  USE pseudo_types,         ONLY : pseudo_upf
  USE mp,                   ONLY : mp_sum, mp_bcast, mp_max, mp_end
  USE mp_global,            ONLY : mpime, nproc, root, me_pool, my_pool_id, &
                                   nproc_pool, intra_pool_comm, root_pool, &
                                   intra_image_comm, my_image_id
  USE io_global,            ONLY : ionode, ionode_id
  USE io_base,              ONLY : read_restart_header, read_restart_ions, &
                                   read_restart_cell, read_restart_electrons, &
                                   read_restart_gvec, read_restart_gkvec, &
                                   read_restart_charge, read_restart_wfc, &
                                   read_restart_symmetry, read_restart_xdim, &
                                   read_restart_pseudo, read_restart_ldaU, &
                                   read_restart_tetra
  USE parameters,           ONLY : nacx, nsx
  USE spin_orb,             ONLY : lspinorb
  USE noncollin_module,     ONLY : noncolin, npol
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_alpha
  USE upf_to_internal

  !
  implicit none
  !
  integer, intent(in) :: ndr, nsizwfc, iunitwfc
  integer, intent(out) :: ierr
  character(len=*), intent(in) :: what
  real(kind=DP), intent(inout) :: et_g(:,:), wg_g(:,:)
  integer, intent(inout) :: kunit
  !
  !
  logical :: exst
  real(kind=DP), allocatable :: stau0(:,:), staum(:,:), amass_(:)
  real(kind=DP), allocatable :: svel0(:,:), svelm(:,:), tautmp(:,:), force_(:,:)
  real(kind=DP) :: xnosp, xnos0, xnosm, xnosm2, cdmi(3)
  character(len=4) :: atom_label(nsx)
  real(kind=DP), dimension(3,3) :: ht0, htm, htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2
  logical :: tscal
  real(kind=DP) :: xenos0, xenosm, xenosm2, xenosp
  integer, allocatable :: ityp_(:)
  character(len=256) :: filename

  integer :: nelu_, neld_, nacx_, kunit_
  integer :: na_(nsx), ngk_l(npk), ngk_g(npk)
  real(kind=DP) :: trutime_, ecutrho_, ekincm_
  real(kind=DP) :: acc_(nacx), bg_(3,3)
  character(len=80) :: title_, crystal_
  character(len=256) :: tmp_dir_
  real(kind=DP), allocatable :: occtmp(:), lambda(:,:)

  integer :: ntau
  integer :: i, ispin, ik, ispin_, nspin_, ik_, nkstot_, nat_, ntyp_
  integer :: nbnd_

  real(kind=DP) :: nelec_
  real(kind=DP) :: xenos0_, xenosm_, xenosm2_, xenosp_
  LOGICAL :: tocc, tlam, teig, tmill, twf0, twfm, tigl
  INTEGER :: ldim, flen
  INTEGER :: nkl, nkr, nkbl, iks, ike
  INTEGER :: npool, ipmask( nproc ), ipsour, ipdest
  LOGICAL :: trdhead, trdxdim, trdcell, trdpos, trdpseudo, trdchden
  LOGICAL :: trdocc, trdgvec, trdgkvec, trdwfc, trdsym, tupf
  INTEGER :: npw, npw_g, ipol, j
  LOGICAL :: lgamma

  real(kind=DP) :: wfc_scal
  COMPLEX(KIND=DP), ALLOCATABLE :: wfc_restart(:,:)

  TYPE( pseudo_upf ) :: upf

  integer, allocatable :: mill(:,:)
  !
  ! ... end of declarations
  !
  !
  ierr = 0
  filename = trim(prefix)//'.save'
  flen = index(filename,' ')-1
  WRITE( stdout, '(/,5x,"Reading file ",a," ... ",$)') filename(1:flen)
  !
  if( ionode ) THEN
    call seqopn (ndr, filename(1:flen), 'unformatted', exst)
    if (.not.exst) then
       close (unit = ndr, status = 'delete')
       ierr = 1
    else
       rewind ndr
    endif   
  end if
  call mp_bcast( ierr, ionode_id, intra_image_comm )
  if( ierr /= 0 ) then
    return
  end if

  trdhead   = .FALSE.
  trdxdim   = .FALSE.
  trdcell   = .FALSE.
  trdpos    = .FALSE.
  trdsym    = .FALSE.
  trdpseudo = .FALSE.
  trdocc    = .FALSE.
  trdgvec   = .FALSE.
  trdgkvec  = .FALSE.
  trdchden  = .FALSE.
  trdwfc    = .FALSE.

  SELECT CASE ( TRIM( what ) )
     CASE ( 'all' )
      WRITE( stdout, '(5x,"all data")') 
      trdhead   = .TRUE.
      trdxdim   = .TRUE.
      trdcell   = .TRUE.
      trdpos    = .TRUE.
      trdsym    = .TRUE.
      trdpseudo = .TRUE.
      trdocc    = .TRUE.
      trdgvec   = .TRUE.
      trdgkvec  = .TRUE.
      trdchden  = .TRUE.
      trdwfc    = .TRUE.
    CASE ( 'nowave' )
      WRITE( stdout, '(5x,"all except wavefuctions")') 
      trdhead   = .TRUE.
      trdxdim   = .TRUE.
      trdcell   = .TRUE.
      trdpos    = .TRUE.
      trdsym    = .TRUE.
      trdpseudo = .TRUE.
      trdocc    = .TRUE.
      trdgvec   = .TRUE.
      trdgkvec  = .TRUE.
      trdchden  = .TRUE.
    CASE ( 'wave' )
      WRITE( stdout, '(5x,"wavefuctions")') 
      trdwfc    = .TRUE.
    CASE ( 'dim' )
      WRITE( stdout, '(5x,"only dimensions")') 
      trdhead   = .TRUE.
      trdxdim   = .TRUE.
    CASE DEFAULT
      WRITE( stdout, '(5x,a)') what
      CALL errore( ' readfile_new ', ' unknown value for what ', 1 )
  END SELECT



!  ==--------------------------------------------------------------==
!  ==  HEADER INFORMATION                                          ==
!  ==--------------------------------------------------------------==

   IF( trdhead ) THEN

     CALL read_restart_header(ndr, istep, iswitch, trutime_, nr1, nr2, nr3, &
       nr1s, nr2s, nr3s, ngm_g, nkstot, ngk_g, nspin, nbnd, nelec, nelu_, &
       neld_, nat, ntyp, na_, acc_, nacx_, ecutwfc, ecutrho_, alat, ekincm_, &
       kunit_, k1, k2, k3, nk1, nk2, nk3, degauss, ngauss, lgauss, ntetra, ltetra, &
       natomwfc, gcutm, gcutms, dual, doublegrid, modenum, lforce, lstres, &
       title_, crystal_, tmp_dir_, tupf, lgamma, noncolin, lspinorb, &
       lda_plus_u,&
       tfixed_occ, tefield, dipfield, edir, emaxpos, eopreg, eamp, twfcollect )

     gamma_only = lgamma
     title = title_
     crystal = crystal_

!
! .. Redistribute k-points according to the new number of processors and pool
! .. this is required in order to calculate the proper value for "nks"
! .. and it is identical to divide_et_impera subroutine, apart that here
! .. we are not dealing with xk and wk ( see G+k vector section )
!

     !  If the blocking factor (kunit) for k points is less than 1 use the
     !  value read from file
     IF( kunit < 1 ) kunit = kunit_

     IF( MOD( nkstot, kunit ) /= 0 ) &
       CALL errore( ' readfile_new ',' inconsistent kunit ',1)

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

     nks = nkl


  ELSE

     CALL read_restart_header(ndr)

   END IF

!  ==--------------------------------------------------------------==
!  ==  MAX DIMENSIONS                                              ==
!  ==--------------------------------------------------------------==

   IF( trdxdim ) THEN

     CALL read_restart_xdim( ndr, npwx, nbndx, nrx1, nrx2, nrx3, &
       nrxx, nrx1s, nrx2s, nrx3s, nrxxs )

   ELSE

     CALL read_restart_xdim( ndr )

   END IF

   IF( what == 'dim' ) GO TO 10

!  ==--------------------------------------------------------------==
!  ==  CELL & METRIC                                               ==
!  ==--------------------------------------------------------------==

   IF( trdcell ) THEN

     CALL read_restart_cell( ndr, ibrav, celldm, ht0, htm, &
          htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2)

     at = TRANSPOSE( ht0 / alat )
     call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2) , bg(1,3) )
     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2
     call volume(alat,at(1,1),at(1,2),at(1,3),omega)

   ELSE

     CALL read_restart_cell( ndr )

   END IF


!  ==--------------------------------------------------------------==
!  ==  IONS                                                        ==
!  ==--------------------------------------------------------------==

   IF( trdpos ) THEN

     ALLOCATE( stau0(3, nat) )
     ALLOCATE( staum(3, nat) )
     ALLOCATE( svel0(3, nat) )
     ALLOCATE( svelm(3, nat) )
     ALLOCATE( tautmp(3, nat) )
     ALLOCATE( force_(3, nat) )
     ALLOCATE( ityp_(nat) )
     ALLOCATE( amass_(ntyp) )

     CALL read_restart_ions(ndr, atom_label(1:ntyp), tscal, stau0, svel0, &
       staum, svelm, tautmp, force_, cdmi, nat_, ntyp_, ityp_, na_, amass_, xnosp,  &
       xnos0, xnosm, xnosm2)
     !
     ! get the right number of "tau" to avoid possible segmentation violations
     !
     ntau  = MIN( nat_, SIZE( tau, 2 ) )
     tau(:,1:ntau)  =  tautmp(:,1:ntau) / alat
     !
     atm(1:ntyp) = atom_label(1:ntyp)
     ityp(1:ntau) = ityp_(1:ntau)
     force(1:3,1:ntau) = force_(1:3,1:ntau)
     DEALLOCATE( stau0, svel0, staum, svelm, tautmp, ityp_, force_, amass_ )

   ELSE

     CALL read_restart_ions(ndr)

   END IF

!  ==--------------------------------------------------------------==
!  ==  LDA+U                                                       ==
!  ==--------------------------------------------------------------==
     IF( lda_plus_u) THEN
        
        CALL read_restart_ldaU(ndr, ntyp, Hubbard_lmax, &
             Hubbard_l(1:ntyp), Hubbard_U(1:ntyp), Hubbard_alpha(1:ntyp))

     ELSE
        
        CALL read_restart_ldaU(ndr)

     ENDIF

!  ==--------------------------------------------------------------==
!  ==  SYMMETRIES                                                  ==
!  ==--------------------------------------------------------------==

   if( trdsym ) then
     CALL read_restart_symmetry( ndr, &
        symm_type, sname, s, irt, nat_, ftau, nsym, invsym, noinv )
   else
     CALL read_restart_symmetry( ndr )
   end if

!  ==--------------------------------------------------------------==
!  ==  PSEUDOPOTENTIALS                                            ==
!  ==--------------------------------------------------------------==

   DO i = 1, ntyp

     if( trdpseudo ) then

       if( tupf ) then

         ALLOCATE( upf%els( nchix ),  upf%lchi( nchix ), upf%jchi(nchix), &
            upf%oc( nchix ), &
            upf%r( ndmx ), upf%rab( ndmx ), upf%rho_atc( ndmx ), &
            upf%vloc( ndmx ), &
            upf%lll( 1:nbrx ), upf%jjj(1:nbrx), upf%kkbeta( 1:nbrx ), &
            upf%beta( ndmx, 1:nbrx ), &
            upf%dion( 1:nbrx, 1:nbrx ), upf%rinner( 1:lqmax ), &
            upf%qqq( 1:nbrx, 1:nbrx ), &
            upf%qfunc( ndmx, 1:nbrx, 1:nbrx ), &
            upf%qfcoef( 1:nqfx, 1:lqmax, 1:nbrx, 1:nbrx ), &
            upf%chi( ndmx, nchix ), upf%rho_at( ndmx ) )

         CALL read_restart_pseudo( ndr, &
           upf%generated, upf%date_author, upf%comment, upf%psd, upf%typ, &
           upf%tvanp,  &
           upf%nlcc, upf%dft, upf%zp, upf%etotps, upf%ecutwfc, &
           upf%ecutrho, upf%nv,   &
           upf%lmax, upf%mesh, upf%nwfc, upf%nbeta, upf%els, upf%lchi, &
           upf%jchi, upf%oc, &
           upf%r, upf%rab, upf%rho_atc, upf%vloc, upf%lll, upf%jjj, &
           upf%kkbeta, &
           upf%beta, upf%nd, upf%dion, upf%nqf, upf%nqlc, upf%rinner, &
           upf%qqq, upf%qfunc, upf%qfcoef, upf%chi, upf%rho_at )

! ...    Here convert to internal scheme
         if (ierr == 0) CALL set_pseudo_upf ( i, upf )

         ! UPF is always numeric
         numeric (i) = .true.
         ! UPF is RRKJ3-like
         newpseudo (i) = .true.

         DEALLOCATE( upf%els,  upf%lchi, upf%jchi, upf%oc, upf%r, upf%rab )
         DEALLOCATE( upf%rho_atc, upf%vloc, upf%lll, upf%jjj, &
             upf%kkbeta, upf%beta, &
             upf%dion, upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef, &
             upf%chi, upf%rho_at )

       else

         CALL read_restart_pseudo( ndr, &
           zmesh(i), xmin(i), dx(i), r(:,i), rab(:,i), vloc_at(:,i), &
           chi(:,:,i), oc(:,i), &
           rho_at(:,i), rho_atc(:,i), mesh(i), msh(i), nchi(i), lchi(:,i), &
           jchi(:,i), &
           numeric(i), cc(:,i), alpc(:,i), zp(i), aps(:,:,i), alps(:,:,i), &
           zv(i), nlc(i), nnl(i), lmax(i), lloc(i), dion(:,:,i), &
           betar(:,:,i), qqq(:,:,i), qfunc(:,:,:,i), qfcoef(:,:,:,:,i), &
           rinner(:,i), nh(i), nbeta(i), kkbeta(i), nqf(i), nqlc(i), &
           ifqopt(i), &
           lll(:,i), jjj(:,i), iver(:,i), tvanp(i), okvan, newpseudo(i), &
           iexch, icorr, &
           igcx, igcc, lsda, a_nlcc(i), b_nlcc(i), alpha_nlcc(i), &
           nlcc(i), psd(i) )
       endif

     else

       CALL read_restart_pseudo( ndr )

     end if

   END DO

!  ==--------------------------------------------------------------==
!  ==  OCCUPATION NUMBER                                           ==
!  ==--------------------------------------------------------------==

   tocc = .FALSE.
   tlam = .FALSE.
   teig = .TRUE.
   ldim = 1

   DO ik = 1, nkstot
     IF ( trdocc ) THEN
       ALLOCATE( occtmp(nbnd) )
       ALLOCATE( lambda(1,1) )
       occtmp = 0.0d0
       CALL read_restart_electrons( ndr, occtmp, occtmp, tocc, lambda, lambda, &
         ldim, tlam, nbnd_, ispin_, nspin_, ik_, nkstot_, nelec_, nelu_, neld_, &
         xenosp_, xenos0_, xenosm_, xenosm2_, ef, teig, et_g(:,ik), wg_g(:,ik) )
       DEALLOCATE( occtmp )
       DEALLOCATE( lambda )
     ELSE
       CALL read_restart_electrons( ndr )
     END IF
   END DO

!  ==--------------------------------------------------------------==
!  ==  G-Vectors                                                   ==
!  ==--------------------------------------------------------------==

   IF ( trdgvec ) THEN
     ALLOCATE( mill(3,1) )
     mill = 0
     tmill = .FALSE.
     CALL read_restart_gvec( ndr, &
       ngm_g, bg_(:,1), bg_(:,2), bg_(:,3), bg_(:,1), bg_(:,2), bg_(:,3), tmill, mill )
     DEALLOCATE( mill )
   ELSE
     CALL read_restart_gvec( ndr )
   END IF

!  ==--------------------------------------------------------------==
!  ==  (G+k)-Vectors                                               ==
!  ==--------------------------------------------------------------==

   DO ik = 1, nkstot
     IF ( trdgkvec ) THEN
       CALL read_restart_gkvec(ndr, &
         ik_, nkstot_, ngk_g(ik), xk(:,ik), wk(ik), isk(ik))
     ELSE
       CALL read_restart_gkvec( ndr )
     END IF
   END DO

   ! .. distribute k points according to the present processor geometry

   nks = nkstot
   IF( nspin == 2 ) lsda = .TRUE.
   CALL divide_et_impera (xk(1,1), wk(1), isk(1), lsda, nkstot, nks)

!  ==--------------------------------------------------------------==
!  ==  Tetrahedra                                                  ==
!  ==--------------------------------------------------------------==

   CALL read_restart_tetra( ndr, ltetra, ntetra, tetra )

!  ==--------------------------------------------------------------==
!  ==  CHARGE DENSITY AND POTENTIALS                               ==
!  ==--------------------------------------------------------------==

   DO ispin = 1, nspin
     CALL read_restart_charge( ndr )
   END DO

!  ==--------------------------------------------------------------==
!  ==  WAVEFUNCTIONS                                               ==
!  ==--------------------------------------------------------------==

   IF( trdwfc ) THEN

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

   END IF

   twf0   = .TRUE.
   twfm   = .FALSE.
   tigl   = .FALSE.

   DO ik = 1, nkstot

     IF ( trdwfc ) THEN

       IF( ( nsizwfc < 1 ) .OR. ( iunitwfc < 1 ) ) &
         CALL errore(' readfile_new ', ' cannot read wave function ', 1 )

       ipmask = 0
       ipdest = root_pool

       !  find out the index of the processor which collect the data in the pool of ik
       IF( npool > 1 ) THEN
         IF( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
           IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
         END IF
         CALL mp_sum( ipmask, intra_image_comm )
         DO i = 1, nproc
           IF( ipmask(i) == 1 ) ipdest = ( i - 1 )
         END DO
       END IF

       IF( (ik >= iks) .AND. (ik <= ike) ) THEN
         CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
         CALL gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,ik-iks+1))
         npw_g = MAXVAL( igk_l2g(:,ik-iks+1) )
         CALL mp_max( npw_g, intra_pool_comm )
       END IF

       CALL mp_bcast( npw_g, ipdest, intra_image_comm )

       if (noncolin) then
         allocate (wfc_restart(npwx, nbnd ))
         wfc_restart = (0.d0,0.d0)
         do ipol = 1, npol
           CALL read_restart_wfc(ndr, ik, nkstot, kunit, ispin_, nspin_, &
                 wfc_scal, wfc_restart, twf0, wfc_restart, twfm, npw_g, &
                 nbnd_, igk_l2g(:,ik-iks+1), npw )
            evc_nc(:,ipol,:)=wfc_restart(:,:)
         enddo
         deallocate (wfc_restart)
       else
         CALL read_restart_wfc(ndr, ik, nkstot, kunit, ispin_, nspin_, &
                        wfc_scal, evc, twf0, evc, twfm, npw_g, nbnd_,  &
                        igk_l2g(:,ik-iks+1), npw )
       endif

       IF( twf0 ) THEN
         IF( (ik >= iks) .AND. (ik <= ike) ) THEN
           IF( wfc_scal /= 1.0d0 ) THEN
             if (noncolin) then
                CALL DSCAL( 2*nsizwfc, wfc_scal, evc_nc, 1)
             else
                CALL DSCAL( 2*nsizwfc, wfc_scal, evc, 1)
             endif
           END IF
             if (noncolin) then
                call davcio (evc_nc(1,1,1), nsizwfc, iunitwfc, (ik-iks+1), 1)
             else
                call davcio (evc(1,1), nsizwfc, iunitwfc, (ik-iks+1), 1)
             endif
         END IF
       ELSE
         ierr = 1
       END IF

     ELSE

       CALL read_restart_wfc(ndr )

     END IF

   END DO


10 continue
   !
   if( ionode ) then
     close (unit = ndr)
   end if
   WRITE( stdout, '(5x,"read complete")')
   !
   return
end subroutine readfile_new


!----------------------------------------------------------------------
subroutine readfile_config( ndr, ibrav, nat, alat, at, tau, ierr )
  !-----------------------------------------------------------------------
  !
  !     This routine is called at the start of the run to read from a file
  !     the information needed to restart and to other postprocessing
  !     programs.
  !
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk, nacx, nsx
  USE io_files,   ONLY : prefix, tmp_dir
  USE io_global,  ONLY : ionode, ionode_id
  USE mp_global,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_bcast
  USE ldaU,       ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, &
                         Hubbard_U, Hubbard_alpha
  USE io_base,    ONLY : read_restart_header, read_restart_ions, &
                         read_restart_cell, read_restart_electrons, &
                         read_restart_gvec, read_restart_gkvec, &
                         read_restart_charge, read_restart_wfc, &
                         read_restart_symmetry, read_restart_xdim, &
                         read_restart_pseudo, read_restart_ldaU
  !
  implicit none
  !
  integer, intent(in) :: ndr
  integer, intent(out)       :: ibrav, nat, ierr
  real(kind=DP), intent(out) :: alat, at(:,:), tau(:,:)
  !
  logical :: tread, exst
  real(kind=DP), allocatable :: stau0(:,:), staum(:,:), amass_(:)
  real(kind=DP), allocatable :: svel0(:,:), svelm(:,:), tautmp(:,:), force_(:,:)
  real(kind=DP) :: xnosp, xnos0, xnosm, xnosm2, cdmi(3)
  character(len=4) :: atom_label( nsx )
  real(kind=DP), dimension(3,3) :: ht0, htm, htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2
  logical :: tscal
  real(kind=DP) :: xenos0, xenosm, xenosm2, xenosp
  integer, allocatable :: ityp_(:)
  character(len=256) :: filename

  integer :: istep_, iswitch_, nr1_, nr2_, nr3_, nr1s_, nr2s_, nr3s_, ngm_, ngmg_, nkstot_
  integer :: nspin_, nbnd_, nelu_, neld_, nat_, ntyp_, nacx_, kunit_, nks_
  integer :: k1_, k2_, k3_, nk1_, nk2_, nk3_, ngauss_
  integer :: na_(nsx), ngk_l(npk), ngk_g(npk)
  integer :: ntetra_, natomwfc_, modenum_
  integer :: Hubbard_lmax_
  integer :: npwx_, nbndx_, nrx1_, nrx2_, nrx3_, nrxx_, nrx1s_, nrx2s_, nrx3s_, nrxxs_
  real(kind=DP) :: trutime_, nelec_, ecutwfc_, ecutrho_, alat_, ekincm_
  real(kind=DP) :: degauss_, gcutm_, gcutms_, dual_
  real(kind=DP) :: acc_(nacx)
  logical :: lgauss_, ltetra_, doublegrid_, lstres_, lforce_, tupf_, &
             lgamma_, noncolin_, lspinorb_, lda_plus_u_
  character(len=80) :: title_, crystal_
  character(len=256) :: tmp_dir_

  real(kind=DP) :: celldm_(6)
  integer :: ibrav_
  integer :: ntau

  logical :: tfixed_occ_, tefield_, dipfield_
  integer :: edir_
  real(kind=DP) :: emaxpos_, eopreg_, eamp_
  logical :: twfcollect_ 
  !
  ! ... end of declarations
  !
  !
  !  read configuration from .save file
  !
  ierr = 0
  filename = trim( prefix )//'.save'
  WRITE( stdout, '(/,5x,"Reading file ",a14)') filename
  !
  if( ionode ) THEN
    call seqopn (ndr, filename, 'unformatted', exst)
    if ( .not. exst ) then
       close (unit = ndr, status = 'delete')
       ierr = 1
    endif
    rewind ndr
  end if
  call mp_bcast( ierr, ionode_id, intra_image_comm )
  !
  !  if the file is not present or unreadable
  !  return immediately
  !
  if( ierr /= 0 ) then
    return
  end if

!  ==--------------------------------------------------------------==
!  ==  HEADER INFORMATION                                          ==
!  ==--------------------------------------------------------------==

   CALL read_restart_header(ndr, istep_, iswitch_, trutime_, nr1_, nr2_, nr3_, &
     nr1s_, nr2s_, nr3s_, ngmg_, nkstot_, ngk_g, nspin_, nbnd_, nelec_, nelu_, &
     neld_, nat_, ntyp_, na_, acc_, nacx_, ecutwfc_, ecutrho_, alat_, ekincm_, &
     kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_, degauss_, ngauss_, lgauss_, &
     ntetra_, ltetra_, &
     natomwfc_, gcutm_, gcutms_, dual_, doublegrid_, modenum_, lstres_, &
     lforce_, title_, crystal_, tmp_dir_, tupf_, lgamma_, noncolin_, lspinorb_, &
     lda_plus_u_, tfixed_occ_, tefield_, dipfield_, edir_, emaxpos_, eopreg_, eamp_, &
     twfcollect_ )

!  ==--------------------------------------------------------------==
!  ==  MAX DIMENSIONS                                              ==
!  ==--------------------------------------------------------------==

   CALL read_restart_xdim( ndr, &
      npwx_, nbndx_, nrx1_, nrx2_, nrx3_, nrxx_, nrx1s_, nrx2s_, nrx3s_, nrxxs_ )


!  ==--------------------------------------------------------------==
!  ==  CELL & METRIC                                               ==
!  ==--------------------------------------------------------------==

   CALL read_restart_cell( ndr, ibrav_, celldm_, ht0, htm, &
          htm2, htvel, xhnosp, xhnos0, xhnosm, xhnosm2)

!  ==--------------------------------------------------------------==
!  ==  IONS                                                        ==
!  ==--------------------------------------------------------------==

   ALLOCATE( stau0(3, nat_) )
   ALLOCATE( staum(3, nat_) )
   ALLOCATE( svel0(3, nat_) )
   ALLOCATE( svelm(3, nat_) )
   ALLOCATE( tautmp(3, nat_) )
   ALLOCATE( force_(3, nat_) )
   ALLOCATE( ityp_(nat_) )
   ALLOCATE( amass_(ntyp_) )

   CALL read_restart_ions(ndr, atom_label(1:ntyp_), tscal, stau0, svel0, &
     staum, svelm, tautmp, force_, cdmi, nat_, ntyp_, ityp_, na_, amass_, xnosp,  &
     xnos0, xnosm, xnosm2)

   ibrav = ibrav_
   nat   = nat_
   alat  = alat_
   at = TRANSPOSE( ht0 / alat )
   !
   ! get the right number of "tau" to avoid possible segmentation violations
   !
   ntau  = MIN( nat_, SIZE( tau, 2 ) )
   tau(:,1:ntau)  = tautmp(:,1:ntau)/alat

   !

   DEALLOCATE( stau0, svel0, staum, svelm, tautmp, ityp_, force_, amass_ )

!  ==--------------------------------------------------------------==
!  ==  LDA+U                                                       ==
!  ==--------------------------------------------------------------==
     IF( lda_plus_u) THEN
        CALL read_restart_ldaU(ndr, ntyp_, Hubbard_lmax_, &
             Hubbard_l(1:ntyp_), Hubbard_U(1:ntyp_), Hubbard_alpha(1:ntyp_))
     ENDIF

   !
   if( ionode ) then
     close (unit = ndr)
   end if
   !
   return
end subroutine readfile_config

end module restart_module
