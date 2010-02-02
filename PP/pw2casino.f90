!
! Copyright (C) 2004-2009 Dario Alfe' and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM pw2casino
  !-----------------------------------------------------------------------

  ! This subroutine writes the file "prefix".pwfn.data containing the
  ! plane wave coefficients and other stuff needed by the QMC code CASINO.
  ! May be useful to anybody desiring to extract data from Quantum ESPRESSO
  ! runs, since the output data is quite easy to understand.
  ! If you want to save the Fermi energy and state occupancies as well,
  ! look at tags KN (courtesy of Karoly Nemeth, Argonne)
  ! Not guaranteed to work in parallel execution ! If you want to read data
  ! written by a parallel run, ensure that the data file was saved in a
  ! portable format (option "wf_collect=.true." for PWscf), run pw2casino
  ! serially. Alternatively: run in the same number of processors and pools
  ! of the previous pw.x calculation.
  ! Usage:
  ! * run first a scf calculation with pw.x
  ! * run pw2casino.x with the following input:
  !      &inputpp prefix='...', outdir ='...' /
  !   where prefix and outdir are the same as those used in the scf calculation
  !   (you may use environment variable ESPRESSO_TMPDIR instead of outdir)
  ! * move all your files named prefix.pwfn.data? to pwfn.data?,
  !   merge the pwfn.data? files using the CASINO utility MERGE_PWFN.
  ! * convert to blips running the BLIP utility.
  ! You do not necessarily have to use casino PP's, but you can if you want;
  ! there is a conversion utility in the upftools directory of the espresso
  ! distribution.

  USE io_files,   ONLY : prefix, outdir, tmp_dir, trimcheck
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, npool, nimage, nogrp, npgrp
  USE environment,ONLY : environment_start
  !
  IMPLICIT NONE
  INTEGER :: ios
  LOGICAL :: casino_gather

  NAMELIST / inputpp / prefix, outdir, casino_gather
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW2CASINO' )

  IF ( npool > 1 .or. nimage > 1) THEN
     CALL errore('pw2casino', 'pool or image parallelization not (yet) implemented')
  ENDIF

  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  casino_gather = .false.
  ios = 0
  IF ( ionode )  THEN
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios/=0 ) CALL errore('pw2casino', 'reading inputpp namelist', ABS(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast(tmp_dir, ionode_id )
  CALL mp_bcast(casino_gather, ionode_id )
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL write_casino_pwfn(casino_gather)
  !
  CALL stop_pp
  STOP

END PROGRAM pw2casino


SUBROUTINE write_casino_pwfn(gather)

  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE printout_base, ONLY: title    ! title of the run
  USE constants, ONLY: tpi, e2
  USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc, demet, ef
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss
  USE lsda_mod, ONLY: lsda, nspin
  USE scf, ONLY: rho, rho_core, rhog_core, vnew
  USE ldaU, ONLY : lda_plus_u, eth, Hubbard_lmax
  USE vlocal, ONLY: vloc, strf
  USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et
  USE control_flags, ONLY : gamma_only
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE becmod, ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_global, ONLY: stdout, ionode, ionode_id
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  USE funct, ONLY : dft_is_meta
  USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm, nproc_pool
  USE mp,                   ONLY: mp_sum, mp_gather, mp_bcast
  USE dfunct,                 only : newd

  IMPLICIT NONE
  LOGICAL, INTENT(in) :: gather
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb, at_num, id, ip
  INTEGER, ALLOCATABLE :: idx(:), igtog(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  REAL(DP), allocatable :: v_h_new(:,:,:,:), kedtaur_new(:,:)
  INTEGER :: ios
  INTEGER, EXTERNAL :: atomic_number
  REAL (DP), EXTERNAL :: ewald, w1gauss

  INTEGER ngtot_g
  INTEGER, ALLOCATABLE :: ngtot_d(:), ngtot_cumsum(:), indx(:)
  REAL(DP), ALLOCATABLE :: g_l(:,:), g_g(:,:), g2(:)
  COMPLEX(DP), ALLOCATABLE :: evc_l(:), evc_g(:)

  CALL init_us_1
  CALL newd

  ALLOCATE (aux(nrxx))
  call allocate_bec_type ( nkb, nbnd, becp )
  ! four times npwx should be enough
  ALLOCATE (idx (4*npwx) )
  ALLOCATE (igtog (4*npwx) )
  if (lda_plus_u) allocate(v_h_new(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
  if (dft_is_meta()) then
     allocate (kedtaur_new(nrxx,nspin))
  else
     allocate (kedtaur_new(1,nspin))
  endif
  idx(:) = 0
  igtog(:) = 0

  IF( lsda ) THEN
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     !     nspin = 2
  ELSE
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     !     nspin = 1
  ENDIF

  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  demet=0.d0
  !
  DO ispin = 1, nspin
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     aux(:) = CMPLX( rho%of_r(:,ispin), 0.d0,kind=DP)
     CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     DO nt=1,ntyp
        DO ig = 1, ngm
           eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                * CONJG(aux(nl(ig)))
        ENDDO
     ENDDO

     DO ik = 1, nk
        ikk = ik + nk*(ispin-1)
        CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, ikk, - 1)
        CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
        CALL calbec ( npw, vkb, evc, becp )
        !
        ! -TS term for metals (if any)
        !
        IF ( degauss > 0.0_dp) THEN
           DO ibnd = 1, nbnd
              demet = demet + wk (ik) * &
                 degauss * w1gauss ( (ef-et(ibnd,ik)) / degauss, ngauss)
           END DO
        END IF
        !
        DO ig =1, npw
           IF( igk(ig) > 4*npwx ) &
                CALL errore ('pw2casino','increase allocation of index', ig)
           idx( igk(ig) ) = 1
        ENDDO
        !
        ! calculate the kinetic energy
        !
        DO ibnd = 1, nbnd
           DO j = 1, npw
              ek = ek +  CONJG(evc(j,ibnd)) * evc(j,ibnd) * &
                              g2kin(j) * wg(ibnd,ikk)
           END DO

           !
           ! Calculate Non-local energy
           !
           ijkb0 = 0
           DO nt = 1, ntyp
              DO na = 1, nat
                 IF (ityp (na) .EQ.nt) THEN
                    DO ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       enl=enl+CONJG(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) &
                            *wg(ibnd,ikk)* dvan(ih,ih,nt)
                       DO jh = ( ih + 1 ), nh(nt)
                          jkb = ijkb0 + jh
                          enl=enl + &
                               (CONJG(becp%k(ikb,ibnd))*becp%k(jkb,ibnd)+&
                                CONJG(becp%k(jkb,ibnd))*becp%k(ikb,ibnd))&
                               * wg(ibnd,ikk) * dvan(ih,jh,nt)

                       END DO

                    ENDDO
                    ijkb0 = ijkb0 + nh (nt)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

#ifdef __PARA
  call mp_sum( eloc,  intra_pool_comm )
  call mp_sum( ek,    intra_pool_comm )
  call mp_sum( ek,    inter_pool_comm )
  call mp_sum( enl,   inter_pool_comm )
  call mp_sum( demet, inter_pool_comm )
#endif
  eloc = eloc * omega
  ek = ek * tpiba2

  ngtot = 0
  DO ig = 1, 4*npwx
     IF( idx(ig) >= 1 ) THEN
        ngtot = ngtot + 1
        igtog(ngtot) = ig
     ENDIF
  ENDDO

  !
  ! compute ewald contribution
  !
  ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
       g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !
  ! compute hartree and xc contribution
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, vnew )
  !
  etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)+demet
  !
  IF (ionode.or..not.gather) THEN

     io = 77

     WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

     CALL seqopn( 77, 'pwfn.data', 'formatted',exst)

     CALL write_header
  ENDIF

  ALLOCATE ( g_l(3,ngtot), evc_l(ngtot) )
  DO ig = 1, ngtot
     g_l(:,ig) = g(:,igtog(ig))
  ENDDO

  IF(gather)THEN
     ALLOCATE ( ngtot_d(nproc_pool), ngtot_cumsum(nproc_pool) )
     CALL mp_gather( ngtot, ngtot_d, ionode_id, intra_pool_comm )
     CALL mp_bcast( ngtot_d, ionode_id, intra_pool_comm )
     id = 0
     DO ip = 1,nproc_pool
        ngtot_cumsum(ip) = id
        id = id + ngtot_d(ip)
     ENDDO
     ngtot_g = id

     ALLOCATE ( g_g(3,ngtot_g), evc_g(ngtot_g) )
     CALL mp_gather( g_l, g_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

     IF (ionode) THEN
        ALLOCATE ( indx(ngtot_g) )
        CALL create_index2(g_g,indx)
        CALL write_gvecs(g_g,indx)
     ENDIF
  ELSE
     ALLOCATE ( indx(ngtot) )
     CALL create_index2(g_l,indx)
     CALL write_gvecs(g_l,indx)
  ENDIF

  IF (ionode.or..not.gather)CALL write_wfn_head

  DO ik = 1, nk
     IF (ionode.or..not.gather)CALL write_kpt_head

     DO ispin = 1, nspin
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, nbnd
           IF (ionode.or..not.gather)CALL write_bnd_head
           evc_l(:) = (0.d0, 0d0)
           DO ig=1, ngtot
              ! now for all G vectors find the PW coefficient for this k-point
              find_ig: DO ig7 = 1, npw
                 IF( igk(ig7) == igtog(ig) )THEN
                    evc_l(ig) = evc(ig7,ibnd)
                    EXIT find_ig
                 ENDIF
              ENDDO find_ig
           ENDDO
           IF(gather)THEN
              CALL mp_gather( evc_l, evc_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

              IF (ionode)CALL write_wfn_data(evc_g,indx)
           ELSE
              CALL write_wfn_data(evc_l,indx)
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  IF (ionode.or..not.gather)CLOSE(io)

  WRITE (stdout,*) 'Kinetic energy   '  , ek/e2
  WRITE (stdout,*) 'Local energy     ', eloc/e2
  WRITE (stdout,*) 'Non-Local energy ', enl/e2
  WRITE (stdout,*) 'Ewald energy     ', ewld/e2
  WRITE (stdout,*) 'xc contribution  ',(etxc-etxcc)/e2
  WRITE (stdout,*) 'hartree energy   ', ehart/e2
  IF ( degauss > 0.0_dp ) &
  WRITE (stdout,*) 'Smearing (-TS)   ', demet/e2
  WRITE (stdout,*) 'Total energy     ', etot/e2

  DEALLOCATE (igtog)
  DEALLOCATE (idx)
  call deallocate_bec_type (becp)
  DEALLOCATE (aux)

  DEALLOCATE ( g_l, evc_l )
  IF(gather) DEALLOCATE ( ngtot_d, ngtot_cumsum, g_g, evc_g )
  IF(ionode.or..not.gather) DEALLOCATE (indx)

CONTAINS

  SUBROUTINE write_header
     WRITE(io,'(a)') title
     WRITE(io,'(a)')
     WRITE(io,'(a)') ' BASIC INFO'
     WRITE(io,'(a)') ' ----------'
     WRITE(io,'(a)') ' Generated by:'
     WRITE(io,'(a)') ' PWSCF'
     WRITE(io,'(a)') ' Method:'
     WRITE(io,'(a)') ' DFT'
     WRITE(io,'(a)') ' DFT Functional:'
     WRITE(io,'(a)') ' unknown'
     WRITE(io,'(a)') ' Pseudopotential'
     WRITE(io,'(a)') ' unknown'
     WRITE(io,'(a)') ' Plane wave cutoff (au)'
     WRITE(io,*) ecutwfc/2
     WRITE(io,'(a)') ' Spin polarized:'
     WRITE(io,*)lsda
     IF ( degauss > 0.0_dp ) THEN
        WRITE(io,'(a)') ' Total energy (au per primitive cell; includes -TS term)'
        WRITE(io,*)etot/e2, demet/e2
     ELSE
        WRITE(io,'(a)') ' Total energy (au per primitive cell)'
        WRITE(io,*)etot/e2
     END IF
     WRITE(io,'(a)') ' Kinetic energy (au per primitive cell)'
     WRITE(io,*)ek/e2
     WRITE(io,'(a)') ' Local potential energy (au per primitive cell)'
     WRITE(io,*)eloc/e2
     WRITE(io,'(a)') ' Non local potential energy(au per primitive cel)'
     WRITE(io,*)enl/e2
     WRITE(io,'(a)') ' Electron electron energy (au per primitive cell)'
     WRITE(io,*)ehart/e2
     WRITE(io,'(a)') ' Ion ion energy (au per primitive cell)'
     WRITE(io,*)ewld/e2
     WRITE(io,'(a)') ' Number of electrons per primitive cell'
     WRITE(io,*)NINT(nelec)
     ! uncomment the following if you want the Fermi energy - KN 2/4/09
     !  WRITE(io,'(a)') ' Fermi energy (au)'
     !  WRITE(io,*) ef/e2
     WRITE(io,'(a)') ' '
     WRITE(io,'(a)') ' GEOMETRY'
     WRITE(io,'(a)') ' -------- '
     WRITE(io,'(a)') ' Number of atoms per primitive cell '
     WRITE(io,*) nat
     WRITE(io,'(a)')' Atomic number and position of the atoms(au) '
     DO na = 1, nat
        nt = ityp(na)
        at_num = atomic_number(TRIM(atm(nt)))
        WRITE(io,'(i6,3f20.14)') at_num, (alat*tau(j,na),j=1,3)
     ENDDO
     WRITE(io,'(a)') ' Primitive lattice vectors (au) '
     WRITE(io,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
     WRITE(io,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
     WRITE(io,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
     WRITE(io,'(a)') ' '

100  FORMAT (3(1x,f20.15))

  END SUBROUTINE


  SUBROUTINE write_gvecs(g,indx)
     REAL(DP),INTENT(in) :: g(:,:)
     INTEGER,INTENT(in) :: indx(:)
     INTEGER ig

     WRITE(io,'(a)') ' G VECTORS'
     WRITE(io,'(a)') ' ---------'
     WRITE(io,'(a)') ' Number of G-vectors'
     WRITE(io,*) SIZE(g,2)
     WRITE(io,'(a)') ' Gx Gy Gz (au)'
     DO ig = 1, SIZE(g,2)
        WRITE(io,100) tpi/alat*g(1,indx(ig)), tpi/alat*g(2,indx(ig)), &
             tpi/alat*g(3,indx(ig))
     ENDDO

100  FORMAT (3(1x,f20.15))

     WRITE(io,'(a)') ' '
  END SUBROUTINE


  SUBROUTINE write_wfn_head
     WRITE(io,'(a)') ' WAVE FUNCTION'
     WRITE(io,'(a)') ' -------------'
     WRITE(io,'(a)') ' Number of k-points'
     WRITE(io,*) nk
  END SUBROUTINE


  SUBROUTINE write_kpt_head
     WRITE(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
          &           k-point coords (au)'
     WRITE(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
          (tpi/alat*xk(j,ik),j=1,3)
  END SUBROUTINE


  SUBROUTINE write_bnd_head
     ! KN: if you want to print occupancies, replace these two lines ...
     WRITE(io,'(a)') ' Band, spin, eigenvalue (au)'
     WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2
     ! ...with the following two - KN 2/4/09
     ! WRITE(io,'(a)') ' Band, spin, eigenvalue (au), occupation number'
     ! WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2, wg(ibnd,ikk)/wk(ikk)
     WRITE(io,'(a)') ' Eigenvectors coefficients'
  END SUBROUTINE


  SUBROUTINE write_wfn_data(evc,indx)
     COMPLEX(DP),INTENT(in) :: evc(:)
     INTEGER,INTENT(in) :: indx(:)
     INTEGER ig

     DO ig=1, SIZE(evc,1)
        WRITE(io,*)evc(indx(ig))
     END DO
  END SUBROUTINE


  SUBROUTINE create_index2(y,x_index)
      DOUBLE PRECISION,INTENT(in) :: y(:,:)
      INTEGER,INTENT(out) :: x_index(SIZE(y,2))
      DOUBLE PRECISION y2(SIZE(y,2))
      INTEGER i
      DO i = 1,SIZE(y,2)
         y2(i) = sum(y(:,i)**2)
      ENDDO
      CALL create_index(y2,x_index)
  END SUBROUTINE


  SUBROUTINE create_index(y,x_index)
 !-----------------------------------------------------------------------------!
 ! This subroutine creates an index array x_index for the n items of data in   !
 ! the array y.  Adapted from Numerical Recipes.                               !
 ! Copied from merge_pwfn.f90, included with CASINO distribution               !
 !-----------------------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: y(:)
  INTEGER,INTENT(out) :: x_index(:)
  INTEGER,PARAMETER :: ins_sort_thresh=7,stacksize=80
  INTEGER n,i,x_indexj,ir,itemp,j,jstack,k,l,lp1,istack(stacksize)
  DOUBLE PRECISION yj
  n=size(x_index)
  do j=1,n
   x_index(j)=j
  enddo ! j
  if(n<=1)return
  jstack=0
  l=1
  ir=n
  do
   if(ir-l<ins_sort_thresh)then
jloop : do j=l+1,ir
     x_indexj=x_index(j) ; yj=y(x_indexj)
     do i=j-1,l,-1
      if(y(x_index(i))<=yj)then
       x_index(i+1)=x_indexj
       cycle jloop
      endif ! y(x_index(i))<=yj
      x_index(i+1)=x_index(i)
     enddo ! i
     x_index(l)=x_indexj
    enddo jloop ! j
    if(jstack==0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
   else
    k=(l+ir)/2
    lp1=l+1
    itemp=x_index(k)    ; x_index(k)=x_index(lp1)  ; x_index(lp1)=itemp
    if(y(x_index(l))>y(x_index(ir)))then
     itemp=x_index(l)   ; x_index(l)=x_index(ir)   ; x_index(ir)=itemp
    endif
    if(y(x_index(lp1))>y(x_index(ir)))then
     itemp=x_index(lp1) ; x_index(lp1)=x_index(ir) ; x_index(ir)=itemp
    endif
    if(y(x_index(l))>y(x_index(lp1)))then
     itemp=x_index(l)   ; x_index(l)=x_index(lp1)  ; x_index(lp1)=itemp
    endif
    i=lp1
    j=ir
    x_indexj=x_index(lp1)
    yj=y(x_indexj)
    do
     do
      i=i+1
      if(y(x_index(i))>=yj)exit
     enddo ! i
     do
      j=j-1
      if(y(x_index(j))<=yj)exit
     enddo ! j
     if(j<i)exit
     itemp=x_index(i) ; x_index(i)=x_index(j) ; x_index(j)=itemp
    enddo
    x_index(lp1)=x_index(j)
    x_index(j)=x_indexj
    jstack=jstack+2
    if(jstack>stacksize)then
     write(6,*)'stacksize is too small.'
     stop
    endif ! jstack>stacksize
    if(ir-i+1>=j-l)then
     istack(jstack)=ir
     istack(jstack-1)=i
     ir=j-1
    else
     istack(jstack)=j-1
     istack(jstack-1)=l
     l=i
    endif ! ir-i+1>=j-l
   endif ! ir-l<ins_sort_thresh
  enddo
  END SUBROUTINE create_index


END SUBROUTINE write_casino_pwfn
