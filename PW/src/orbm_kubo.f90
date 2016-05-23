!==============================================================================!
!
! Copyright (C) 2001-2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! this routine is used to calculate the Kubo terms
! of orbital magnetization (in SI, i.e., A/m)
! written by Andrei Malashevich at UC Berkeley
! For details see
! New. J. Phys. 12, 053032 (2010)
! Many parts from bp_c_phase.f90 and 
! h_epsi_her_set.f90 are reused

! NOTES:
! 
! In order to compute Kubo terms one must first perform a usual SCF 
! calculation, then NSCF calculation with flag 
! lorbm=.true. (in the "control" section) using a UNIFORM grid 
! of k-points.
!
!==============================================================================!

SUBROUTINE orbm_kubo()

!------------------------------------------------------------------------------!

  !  --- Make use of the module with common information ---
  USE ener, ONLY : ef
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : noncolin, npol
  USE wvfct,                ONLY : npwx, nbnd, et, current_k
  USE lsda_mod,             ONLY : nspin
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm,ngm_g,g,gcutm,ig_l2g
  USE start_k,              ONLY : nk1, nk2, nk3
  USE klist,                ONLY : nks,xk,ngk, igk_k
  USE cell_base,            ONLY : tpiba,gpar=>bg,at,alat,omega
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE constants,            ONLY : pi, tpi,rytoev
  USE bp,                   ONLY : lelfield,mapgp_global,mapgm_global,nx_el
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE uspp,                 ONLY : nkb,vkb
  USE scf,                  ONLY : vrs, vltot, v, kedtau
  USE gvecs,                ONLY : doublegrid
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp_world,             ONLY : world_comm

  !  --- Avoid implicit definitions ---
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: evc_k(:,:)!for wavefunctios at k
  COMPLEX(DP), ALLOCATABLE :: evc_kp(:,:)!for wavefunctios at k'
  COMPLEX(DP), ALLOCATABLE :: aux_k(:)
  COMPLEX(DP), ALLOCATABLE :: aux_kp(:)
  COMPLEX(DP), ALLOCATABLE :: aux_kp_g(:)
  COMPLEX(DP), ALLOCATABLE :: evcpm(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: H_evc(:,:)
  COMPLEX(DP), ALLOCATABLE :: temp(:),temp2(:)
  COMPLEX(DP) :: store1, store2
  COMPLEX(DP) :: sca
  ! map g-space global to g-space k-point dependent
  INTEGER, ALLOCATABLE :: ln(:,:,:)
  INTEGER, ALLOCATABLE  :: map_g(:)
  INTEGER :: i,j,k,n,np ! Numbering of k-points
  ! np (n') is used in the loop over neigboring k-points
  INTEGER :: signum
  INTEGER :: tmp
  INTEGER :: ipol ! from 1 to npol
  INTEGER :: istart, iend ! ranges of some arrays
  LOGICAL :: inbz(6) ! if true k' is in BZ
  REAL(DP) :: gtr(3) ! G+G_0
  INTEGER :: npw_k, npw_kp
  INTEGER :: igk_kp(npwx)
  INTEGER :: nb, mb
  INTEGER :: ig
  INTEGER :: n1,n2,n3,ng
  COMPLEX(DP) :: mat(nbnd,nbnd)
  REAL(DP) :: eps ! small number
  ! For matrix inversion (zgedi)
  INTEGER :: ivpt(nbnd)
  INTEGER :: info
  COMPLEX(DP) :: cdet(2)
  COMPLEX(DP) :: cdwork(nbnd)
  REAL(DP) :: mlc(3),mic(3) ! orbital magnetization (LC and IC terms)
  COMPLEX(DP) :: zdotc
  INTEGER :: kpt_arr(3) ! k-point mesh
  INTEGER :: eps_i(3) ! these play role of the antisymmetric tensor e_ijk
  INTEGER :: eps_j(3)
  INTEGER :: sig, sigp
  INTEGER :: l
  REAL(DP) :: pref ! prefactor for MAGNETIZATION in SI
  REAL(DP) :: pref_bm ! prefactor for MAGNETIC MOMENT per cell in Bohr magnetons
  REAL(DP) :: pbm ! prefactor for MAGNETIC MOMENT per cell in Bohr magnetons
  REAL(DP), PARAMETER :: el_si=1.60217646E-19 ! electron charge (SI)
  REAL(DP), PARAMETER :: hbar_si=1.054571628E-34 ! hbar (SI)
  REAL(DP), PARAMETER :: bohr_si=5.2917720859E-11 ! Bohr radius in m
  REAL(DP), PARAMETER :: ry_si=2.179871993E-18 ! Rydberg in J (energy)
  REAL(DP), PARAMETER :: ry_ev=13.6056923 ! Rydberg in eV (energy)
  LOGICAL :: store_flag
  INTEGER :: nbr(6) ! map for 6 neighboring k-points

  ! Allocate necessary arrays
  ALLOCATE(evc_k(npwx*npol,nbnd))
  ALLOCATE(evc_kp(npwx*npol,nbnd))
  ALLOCATE(map_g(npwx))
  ALLOCATE(ln(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3) )
  ALLOCATE(aux_k(ngm*npol))
  ALLOCATE(aux_kp(ngm*npol))
  ALLOCATE(evcpm(npol*npwx,nbnd,6))
  ALLOCATE(H_evc(npol*npwx,nbnd))
  ALLOCATE(temp(ngm))

  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  CALL allocate_bec_type ( nkb, nbnd, becp )
  ! Initializations

  ! Define small number
  eps=1.0d-6

  mlc=0.0d0
  mic=0.0d0

  kpt_arr(1)=nk1
  kpt_arr(2)=nk2
  kpt_arr(3)=nk3

  eps_i(1)=2
  eps_i(2)=3
  eps_i(3)=1
  eps_j(1)=3
  eps_j(2)=1
  eps_j(3)=2

  ! convert energy from Ry to J
  ! alat is in a.u. (Bohr) need to convert to SI
  pref=ry_si*el_si/hbar_si/4.0_dp/(tpi**3)*tpiba/bohr_si
  ! magnetic moment in Bohr magnetons
  ! need to multiply by unit-cell volume omega
  ! convert Ry to Ha by dividing by 2.0
  ! Bohr magneton in a.u. is 1/2
  ! so these two factors cancel out
  ! e=hbar=1 so forget about it
  ! the rest is in atomic units already
  pref_bm=omega/4.0_dp/(tpi**3)*tpiba
  pbm=pref_bm/pref

  !--- Recalculate FFT correspondence (see ggen.f90) ---
  ln=0
  DO ng=1,ngm
    n1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
    n2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
    n3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
    ln(n1,n2,n3) = ng
  END DO

  DO i=1,nk1  ! x
    DO j=1,nk2 ! y
      DO k=1,nk3 ! z
        ! Consecutive ordering of k-points
        n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1

        ! Read wavefunction at k
        CALL get_buffer ( evc_k, nwordwfc, iunwfc, n )
        ! set current k-point (n), recompute projectors, kinetic energy
        ! needed by h_psi
        npw_k = ngk(n)
        current_k=n
        CALL init_us_2(npw_k,igk_k(1,n),xk(1,n),vkb)
        CALL g2_kin( n )

        evcpm=(0.0d0,0.0d0)

        !====================================================!
        !=== Compute dual vectors ===========================!
        !====================================================!

        inbz=.false.
        ! Find indices of neighboring k-points
        ! Current point
        ! n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
        !
        ! k'=k-dx
        IF(i>1) THEN
          !np = (k-1) + (j-1)*nk3 + (i-2)*nk2*nk3 + 1
          nbr(1)=n-nk2*nk3
          inbz(1)=.true.
        ELSE
          !np = (k-1) + (j-1)*nk3 + (nk1-1)*nk2*nk3 + 1
          nbr(1)=n+(nk1-1)*nk2*nk3
        END IF
        ! k'=k+dx
        IF(i<nk1) THEN
          !np = (k-1) + (j-1)*nk3 + i*nk2*nk3 + 1
          nbr(2)=n+nk2*nk3
          inbz(2)=.true.
        ELSE
          !np = (k-1) + (j-1)*nk3 + 1
          nbr(2)=n-(nk1-1)*nk2*nk3
        END IF
        ! k'=k-dy
        IF(j>1) THEN
          !np = (k-1) + (j-2)*nk3 + (i-1)*nk2*nk3 + 1
          nbr(3)=n-nk3
          inbz(3)=.true.
        ELSE
          !np = (k-1) + (nk2-1)*nk3 + (nk1-1)*nk2*nk3 + 1
          nbr(3)=n+(nk2-1)*nk3
        END IF
        ! k'=k+dy
        IF(j<nk2) THEN
          !np = (k-1) + j*nk3 + (i-1)*nk2*nk3 + 1
          nbr(4)=n+nk3
          inbz(4)=.true.
        ELSE
          !np = (k-1) + (i-1)*nk2*nk3 + 1
          nbr(4)=n-(nk2-1)*nk3
        END IF
        ! k'=k-dz
        IF(k>1) THEN
          !np = (k-2) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
          nbr(5)=n-1
          inbz(5)=.true.
        ELSE
          !np = (nk3-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
          nbr(5)=n+(nk3-1)
        END IF
        ! k'=k+dz
        IF(k<nk3) THEN
          !np = k + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
          nbr(6)=n+1
          inbz(6)=.true.
        ELSE
          !np = (j-1)*nk3 + (i-1)*nk2*nk3 + 1
          nbr(6)=n-(nk3-1)
        END IF
        ! NOTE: nbr(5) and nbr(6) are neighbors along z 
        !       and require special treatment in parallel case

        ! loop over neighbors
        signum=1
        DO np=1,6
          signum=-signum

          npw_kp = ngk(nbr(np))
          igk_kp(:)= igk_k(:,nbr(np))
          CALL get_buffer ( evc_kp, nwordwfc, iunwfc, nbr(np) )

          ! Calculate S-1(k,k-dx)

          ! --- Matrix elements calculation ---
          mat=(0.d0,0.d0)

          IF (inbz(np)) THEN

            DO nb=1,nbnd
              DO mb=1,nbnd
                aux_k=(0.d0,0.d0)
                aux_kp=(0.d0,0.d0)

                DO ipol=1,npol

                  istart = (ipol-1)*npwx+1
                  iend = istart+npw_k-1
                  aux_k(igk_k(1:npw_k,n)+ngm*(ipol-1))=evc_k(istart:iend,nb)

                  iend = istart+npw_kp-1
                  aux_kp(igk_kp(1:npw_kp)+ngm*(ipol-1))=evc_kp(istart:iend,mb)

                END DO

                mat(nb,mb) = zdotc(ngm*npol,aux_k,1,aux_kp,1)

              END DO
            END DO
            CALL mp_sum( mat, intra_pool_comm )

          ELSE ! (.not. inbz)

            IF((ngm==ngm_g).OR.(np>4)) THEN ! regular treatment, same for serial and parallel
              map_g=0
              DO ig=1,npw_kp
                !--- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
                !--- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---
                !--- or psi_k'(G)=psi_k(G+G_0)
                ! np=1,3,5 sign "+"
                ! np=2,4,6 sign "-" (use signum for this purpose)
                ! np=1,2 gpar(:,1)
                ! np=3,4 gpar(:,2)
                ! np=5,6 gpar(:,3) 
                ! (use (np+1)/2 for this purpose, note integer arithmetic)
                gtr(1)=g(1,igk_kp(ig)) - DBLE(signum) * gpar(1,(np+1)/2)
                gtr(2)=g(2,igk_kp(ig)) - DBLE(signum) * gpar(2,(np+1)/2)
                gtr(3)=g(3,igk_kp(ig)) - DBLE(signum) * gpar(3,(np+1)/2)
                !--- Find crystal coordinates of gtr, n1,n2,n3 ---
                !--- and the position ng in the ngm array ---
                IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                  n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1)+gtr(3)*at(3,1))
                  n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2)+gtr(3)*at(3,2))
                  n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3)+gtr(3)*at(3,3))
                  ng=ln(n1,n2,n3)
                  IF ((ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                      (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                      (ABS(g(3,ng)-gtr(3)) > eps)) THEN
                    WRITE(6,*) ' error hepsiher: translated G=', &
                      gtr(1),gtr(2),gtr(3), &
                      ' with crystal coordinates',n1,n2,n3, &
                      ' corresponds to ng=',ng,' but G(ng)=', &
                      g(1,ng),g(2,ng),g(3,ng)
                    WRITE(6,*) ' probably because G_par is NOT', &
                      ' a reciprocal lattice vector '
                    !WRITE(6,*) 'DBGG: n,np=',n,np
                    STOP
                  ENDIF
                ELSE
                  WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                    gtr(1),gtr(2),gtr(3)
                  STOP
                END IF
                map_g(ig)=ng
              END DO
  
            END IF ! regular treatment

              DO mb=1,nbnd
                IF((ngm==ngm_g).OR.(np>4)) THEN ! regular treatment, same for serial and parallel
            DO nb=1,nbnd
                  aux_k=(0.d0,0.d0)
                  aux_kp=(0.d0,0.d0)

                  DO ipol=1,npol

                    istart = (ipol-1)*npwx+1
                    iend = istart+npw_k-1
                    aux_k(igk_k(1:npw_k,n)+ngm*(ipol-1))=evc_k(istart:iend,nb)

                    iend = istart+npw_kp-1
                    aux_kp(map_g(1:npw_kp)+ngm*(ipol-1))=evc_kp(istart:iend,mb)

                  END DO

                  mat(nb,mb) = zdotc(ngm*npol,aux_k,1,aux_kp,1)
            END DO
  
                ELSE ! Special parallel treatment
                  ! allocate global array
  
                  ALLOCATE(aux_kp_g(ngm_g*npol))
                  aux_kp_g=(0.0d0,0.0d0)
                  DO ipol=1,npol
                    istart = (ipol-1)*npwx+1
                    iend = istart+npw_kp-1
                    IF(np==1) THEN
                      aux_kp_g(mapgp_global(ig_l2g(igk_kp(1:npw_kp)),1)+ngm_g*(ipol-1))= &
                        evc_kp(istart:iend,mb)
                    END IF
                    IF(np==2) THEN
                      aux_kp_g(mapgm_global(ig_l2g(igk_kp(1:npw_kp)),1)+ngm_g*(ipol-1))= &
                        evc_kp(istart:iend,mb)
                    END IF
                    IF(np==3) THEN
                      aux_kp_g(mapgp_global(ig_l2g(igk_kp(1:npw_kp)),2)+ngm_g*(ipol-1))= &
                        evc_kp(istart:iend,mb)
                    END IF
                    IF(np==4) THEN
                      aux_kp_g(mapgm_global(ig_l2g(igk_kp(1:npw_kp)),2)+ngm_g*(ipol-1))= &
                        evc_kp(istart:iend,mb)
                    END IF
                  END DO
                  CALL mp_sum(aux_kp_g(:),world_comm)
            DO nb=1,nbnd
                  sca=(0.0d0,0.0d0)
                  DO ipol=1,npol
                    DO ig=1,npw_k
                      sca=sca+CONJG(evc_k(ig+npwx*(ipol-1),nb))*&
                              aux_kp_g(ig_l2g(igk_k(ig,n))+ngm_g*(ipol-1))
                    END DO
                  END DO
                  mat(nb,mb)=sca
            END DO
                  DEALLOCATE(aux_kp_g)
                END IF ! parallel treatment
              END DO
            CALL  mp_sum( mat, intra_pool_comm )
          END IF

          !--- Calculate matrix inverse ---
          CALL zgefa(mat,nbnd,nbnd,ivpt,info)
          CALL errore('orbm_kubo','error in zgefa',abs(info))
          CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)

          DO nb=1,nbnd
            DO ipol=1,npol
              temp=(0.0d0,0.0d0)
              istart = (ipol-1)*npwx+1
              iend = istart+npw_kp-1
              ! map_g is needed only if kp is outside of BZ
              ! otherwise use igk_kp
              IF (inbz(np)) THEN
                temp(igk_kp(1:npw_kp))=evc_kp(istart:iend,nb)
              ELSE
                IF((ngm==ngm_g).OR.(np>4)) THEN ! regular treatment
                  temp(map_g(1:npw_kp))=evc_kp(istart:iend,nb)
                ELSE ! map_g is not defined
                  ALLOCATE(temp2(ngm_g))
                  temp2=(0.0d0,0.0d0) !playing role of temp above
                  IF(np==1) THEN
                    temp2(mapgp_global(ig_l2g(igk_kp(1:npw_kp)),1))=&
                      evc_kp(istart:iend,nb)
                  END IF
                  IF(np==2) THEN
                    temp2(mapgm_global(ig_l2g(igk_kp(1:npw_kp)),1))=&
                      evc_kp(istart:iend,nb)
                  END IF
                  IF(np==3) THEN
                    temp2(mapgp_global(ig_l2g(igk_kp(1:npw_kp)),2))=&
                      evc_kp(istart:iend,nb)
                  END IF
                  IF(np==4) THEN
                    temp2(mapgm_global(ig_l2g(igk_kp(1:npw_kp)),2))=&
                      evc_kp(istart:iend,nb)
                  END IF
                  CALL mp_sum(temp2(:),world_comm)
                END IF
              END IF
              iend = istart+npw_k-1
              DO mb=1,nbnd
                IF(inbz(np).OR.(np>4).OR.(ngm==ngm_g)) THEN
                  evcpm(istart:iend,mb,np)=evcpm(istart:iend,mb,np)+&
                      mat(nb,mb)*temp(igk_k(1:npw_k,n))
                ELSE ! special parallel case
                  evcpm(istart:iend,mb,np)=evcpm(istart:iend,mb,np)+&
                      mat(nb,mb)*temp2(ig_l2g(igk_k(1:npw_k,n)))
                END IF
              END DO
              IF(ALLOCATED(temp2)) DEALLOCATE(temp2)
            END DO
          END DO

        END DO ! loop over neighbors

        !====================================================!
        !=== Compute orbital magnetization ==================!
        !====================================================!
        npw_kp = ngk(nbr(1))
        igk_kp(:) = igk_k(:,nbr(1))

        ! LC TERM

        DO l=1,3 ! loop over gpar's
          DO sig=0,1 ! i -/+
            DO sigp=0,1 ! j -/+
              IF(sig==sigp) THEN
                signum=1
              ELSE
                signum=-1
              END IF

              ! H | u_{nk j sigp} >
              H_evc=(0.0d0,0.0d0)
              store_flag=lelfield
              lelfield=.false.
              CALL h_psi(npwx, npw_k, nbnd, evcpm(:,:,2*eps_j(l)+sigp-1), H_evc)
              lelfield=store_flag

              DO nb=1,nbnd ! loop over bands
                  mlc=mlc+DBLE(signum)*pref*gpar(:,l)/kpt_arr(l)* &
                    AIMAG( zdotc(npwx*npol,evcpm(:,nb,2*eps_i(l)+sig-1),1,H_evc(:,nb),1) )
              END DO

            END DO
          END DO
        END DO

        ! IC TERM

        DO l=1,3 ! loop over gpar's
          DO sig=0,1 ! i +/-
            DO sigp=0,1 ! j +/-
              IF(sig==sigp) THEN
                signum=1
              ELSE
                signum=-1
              END IF

              H_evc=(0.0d0,0.0d0)
              store_flag=lelfield
              lelfield=.false.
              CALL h_psi(npwx, npw_k, nbnd, evc_k, H_evc)
              lelfield=store_flag

              DO nb=1,nbnd ! loop over bands
                DO mb=1,nbnd ! loop over bands
                  store1=zdotc(npw_k,evc_k(1:npw_k,nb),1,H_evc(1:npw_k,mb),1)
                  store2=zdotc(npw_k,evcpm(1:npw_k,mb,2*eps_i(l)+sig-1),1, &
                                     evcpm(1:npw_k,nb,2*eps_j(l)+sigp-1),1)
                  IF(noncolin) THEN
                    store1=store1+zdotc(npw_k,evc_k(npwx+1:npwx+npw_k,nb),1,H_evc(npwx+1:npwx+npw_k,mb),1)
                    store2=store2+zdotc(npw_k,evcpm(npwx+1:npwx+npw_k,mb,2*eps_i(l)+sig-1),1, &
                                              evcpm(npwx+1:npwx+npw_k,nb,2*eps_j(l)+sigp-1),1)
                  END IF
                  CALL mp_sum(store1,world_comm)
                  CALL mp_sum(store2,world_comm)
                    mic=mic+DBLE(signum)*pref*gpar(:,l)/kpt_arr(l)*AIMAG(store1*store2)
                END DO
              END DO
            END DO
          END DO
        END DO

      END DO ! loop over k-points
    END DO
  END DO
  CALL mp_sum(mlc,world_comm)

  WRITE (stdout,*) ' '
  WRITE (stdout,*) '=============================================='
  WRITE (stdout,*) '=     ORBITAL MAGNETIZATION (KUBO TERMS)     ='
  WRITE (stdout,*) '=============================================='
  WRITE (stdout,*) ' '
  WRITE (stdout,*) '=     Local circulation term                 ='
  WRITE (stdout,*) 'M_LC = ', mlc(1), mlc(2), mlc(3),' (A/m)'
  WRITE (stdout,*) 'M_LC = ', mlc(1)*pbm, mlc(2)*pbm, mlc(3)*pbm,' (Bohr mag/cell)'
  WRITE (stdout,*) ' '
  WRITE (stdout,*) '=     Itinerant circulation term             ='
  WRITE (stdout,*) 'M_IC = ', mic(1), mic(2), mic(3),' (A/m)'
  WRITE (stdout,*) 'M_IC = ', mic(1)*pbm, mic(2)*pbm, mic(3)*pbm,' (Bohr mag/cell)'
  WRITE (stdout,*) ' '
  WRITE (stdout,*) '=============================================='
  WRITE (stdout,*) ' '

  ! Deallocate arrays
  CALL deallocate_bec_type ( becp )
  DEALLOCATE(temp)
  DEALLOCATE(evc_k)
  DEALLOCATE(evc_kp)
  DEALLOCATE(aux_k)
  DEALLOCATE(aux_kp)
  DEALLOCATE(ln)
  DEALLOCATE(map_g)
  DEALLOCATE(evcpm)
  DEALLOCATE(H_evc)

END SUBROUTINE orbm_kubo
!==============================================================================!
