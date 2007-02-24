!
! Copyright (C) 2006 quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE sym_band(filband, spin_component, firstk, lastk)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba2, at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : nrx1, nrx2, nrx3, nrxx, nr1, nr2, &
                                   nr3, ngm, nl, g, ecutwfc
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE symme,                ONLY : s, ftau, nsym, t_rev
  USE char,                 ONLY : sname
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
                                   char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
                                   which_irr_so, char_mat_so, name_rap_so, &
                                   name_class_so, d_spin, name_class_so1
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ftau_is, gname_is, &
                                   sname_is, code_group_is
  USE uspp,                 ONLY : nkb, vkb
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j, irot, iclass, code_group_old, ig, ibnd
  INTEGER :: spin_component, nks1, nks2, firstk, lastk
  INTEGER :: iunout, ios
  INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), sk_is(3,3,48), &
             gk_is(3,48), t_revk(48), nsymk, isym, ipol, jpol
  LOGICAL :: is_complex, is_complex_so, search_sym, is_symmorphic, is_high_sym
  REAL(DP) :: sr(3,3,48), ft1, ft2, ft3
  COMPLEX(DP) :: d_spink(2,2,48), d_spin_is(2,2,48), ZDOTC
  INTEGER, ALLOCATABLE :: rap_et(:)
  CHARACTER(LEN=45) :: snamek(48)
  CHARACTER (LEN=256) :: filband, namefile
  !
  ALLOCATE(rap_et(nbnd))
  IF (nspin==1.OR.nspin==4) THEN
     nks1=MAX(1,firstk)
     nks2=MIN(nkstot, lastk)
     IF (spin_component .ne. 1)  &
        CALL errore('punch_bands','uncorrect spin_component',1)
  ELSE IF (nspin.eq.2) THEN
     IF (spin_component == 1) THEN
        nks1=MAX(1,firstk)
        nks2=MIN(nks/2,lastk)
     ELSE IF (spin_component==2) THEN
        nks1=nks/2 + MAX(1,firstk)
        nks2=nks/2 + MIN(nks/2,lastk)
     ELSE
        CALL errore('punch_bands','uncorrect spin_component',1)
     END IF
  END IF

  IF ( ionode ) THEN
     iunout=58
     namefile=TRIM(filband)//".rap"
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', err = 200, iostat = ios)
200  CALL errore ('sym_band', 'Opening filband file', ABS (ios) )
     REWIND (iunout)
  ENDIF

  code_group_old=0
  is_high_sym=.FALSE.
  DO ik = nks1, nks2
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     ! 
     ! Find the small group of k
     !
     CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, nsym, sk, ftauk, &
                   gk, t_revk, snamek, nsymk)

     is_symmorphic=.TRUE.
     DO isym=1,nsymk
        is_symmorphic=( is_symmorphic.AND.(ftauk(1,isym)==0).AND.  &
                                          (ftauk(2,isym)==0).AND.  &
                                          (ftauk(3,isym)==0) )
     END DO
     search_sym=.TRUE.
     IF (.NOT.is_symmorphic) THEN
        DO isym=1,nsymk
           search_sym=( search_sym.AND.(gk(1,isym)==0).AND.  &
                                       (gk(2,isym)==0).AND.  &
                                       (gk(3,isym)==0) )
        END DO
     END IF
!
!  Set the group name, divide it in classes and set the
!  character of the irreducible representations
!
     nsym_is=0
     DO isym=1,nsymk
        CALL s_axis_to_cart (sk(1,1,isym), sr(1,1,isym), at, bg)
        IF (noncolin) THEN
!
!  In the noncollinear magnetic case finds the invariant subgroup of the point
!  group of k. Presently we use only this subgroup to classify the levels.
!
           IF (domag) THEN
              IF (t_revk(isym)==0) THEN
                 nsym_is=nsym_is+1
                 CALL s_axis_to_cart (sk(1,1,isym), sr_is(1,1,nsym_is), at, bg)
                 CALL find_u(sr_is(1,1,nsym_is),d_spin_is(1,1,nsym_is))
                 sk_is(:,:,nsym_is)=sk(:,:,isym)
                 gk_is(:,nsym_is)=gk(:,isym)
                 ftau_is(:,nsym_is)=ftauk(:,isym)
                 sname_is(nsym_is)=snamek(isym)
              END IF
           ELSE
              CALL find_u(sr(1,1,isym),d_spink(1,1,isym))
           END IF
        END IF
     END DO
     CALL find_group(nsymk,sr,gname,code_group)
     IF (noncolin) THEN
        IF (domag) THEN
           CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
           CALL set_irr_rap_so(code_group_is,nclass,nrap,char_mat_so, &
                              name_rap_so,name_class_so,name_class_so1)
           CALL divide_class_so(code_group_is,nsym_is,sr_is,d_spin_is,&
                                has_e,nclass,nelem_so,elem_so,which_irr_so)
        ELSE
           CALL set_irr_rap_so(code_group,nclass,nrap,char_mat_so, &
                           name_rap_so,name_class_so,name_class_so1)
           CALL divide_class_so(code_group,nsymk,sr,d_spink,has_e,nclass,  &
                            nelem_so,elem_so,which_irr_so)
        END IF
     ELSE
        CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
        CALL divide_class(code_group,nsymk,sr,nclass,nelem,elem,which_irr)
     END IF

     WRITE(stdout, '(/,1x,74("*"))')
     WRITE(stdout, '(/,20x,"xk=(",2(f10.5,","),f10.5,"  )")') &
                                       xk(1,ik), xk(2,ik), xk(3,ik)
     IF (.not.search_sym) THEN
        WRITE(stdout,'(/,5x,"zone border point and non-symmorphic group ")') 
        WRITE(stdout,'(5x,"symmetry decomposition not available")') 
        WRITE( stdout, '(/,1x,74("*"))')
        rap_et=-1
        GOTO 100
     ENDIF

     IF (code_group.ne.code_group_old) THEN
        CALL write_group_info(.true.)
        is_high_sym=.NOT.is_high_sym
     ELSE
        is_high_sym=.FALSE.
     ENDIF
     IF (noncolin) THEN
        IF (domag) THEN
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsym_is, &
              ngm,sk_is,ftau_is,d_spin_is,gk_is,xk(1,ik),igk,nl,nr1,nr2,& 
              nr3,nrx1,nrx2,nrx3,nrxx,npwx,rap_et )
        ELSE
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsymk,ngm, &
              sk,ftauk,d_spink,gk,xk(1,ik),igk,nl,nr1,nr2,nr3,nrx1, &
              nrx2,nrx3,nrxx,npwx,rap_et)
        ENDIF
     ELSE
        CALL find_band_sym (evc, et(1,ik), at, nbnd, npw, nsymk, ngm, &
           sk, ftauk, gk, xk(1,ik), igk, nl, nr1, nr2, nr3, nrx1, &
           nrx2, nrx3, nrxx, npwx, rap_et )
     END IF

     code_group_old=code_group
100  CONTINUE
     IF (ionode) THEN
        IF (ik == nks1) &
           WRITE (iunout, '(" &plot_rap nbnd_rap=",i4,", nks_rap=",i4," /")') &
                 nbnd, nks2-nks1+1
        WRITE (iunout, '(10x,3f10.6,l5)') xk(1,ik),xk(2,ik),xk(3,ik),is_high_sym
        WRITE (iunout, '(10i8)') (rap_et(ibnd), ibnd=1,nbnd)
     ENDIF
  END DO

  IF (ionode) THEN
     CLOSE(iunout)
  END IF
  !
  DEALLOCATE(rap_et)
  RETURN
END SUBROUTINE sym_band
!
SUBROUTINE find_band_sym (evc,et,at,nbnd,npw,nsym,ngm,s,ftau,gk, &
                     xk,igk,nl,nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,npwx,rap_et)
!
!   This subroutine finds the irreducible representations which give
!   the transformation properties of the wavefunctions. 
!   Presently it does NOT work at zone border if the space group of
!   the crystal has fractionary translations (non-symmorphic space groups).
!  
!
USE io_global,       ONLY : stdout
USE kinds,           ONLY : DP
USE constants,       ONLY : rytoev
USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
                            char_mat, name_rap, name_class, gname
USE uspp,            ONLY : vkb, nkb, okvan
USE becmod,          ONLY : becp
IMPLICIT NONE

INTEGER :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, npw, npwx

INTEGER ::                  &
          nsym,             & 
          nbnd,             &
          rap_et(nbnd),     &
          igk(npwx),        &
          nl(ngm),          &
          ftau(3,48),       &
          gk(3,48),         &
          s(3,3,48)  

REAL(DP) ::                 &
            at(3,3),        &
            xk(3),          &
            et(nbnd)       

COMPLEX(DP) ::  &
            evc(npwx, nbnd)       

REAL(DP), PARAMETER :: eps=1.d-5

INTEGER ::         &
        ngroup,    &  ! number of different frequencies groups
        ibnd,      &
        igroup,    &
        dim_rap,   &
        irot,      &
        irap,      &
        iclass,    &
        shift,     &
        na, i, j, ig, dimen

INTEGER, ALLOCATABLE :: istart(:)

COMPLEX(DP) :: ZDOTC, times              ! safe dimension 
                                         ! in case of accidental degeneracy 
REAL(DP), ALLOCATABLE ::  w1(:)
COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), trace(:,:)
!
!    Divide the bands on the basis of the band degeneracy.
!
ALLOCATE(istart(nbnd+1))
ALLOCATE(w1(nbnd))
ALLOCATE(evcr(npwx,nbnd))
ALLOCATE(trace(48,nbnd))
IF (okvan) ALLOCATE(becp(nkb,nbnd))

rap_et=-1
w1=et*rytoev

ngroup=1
istart(ngroup)=1
DO ibnd=2,nbnd
   IF (ABS(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
      ngroup=ngroup+1
      istart(ngroup)=ibnd
   END IF
END DO
istart(ngroup+1)=nbnd+1
!
!  Find the character of one symmetry operation per class
!
DO iclass=1,nclass
   irot=elem(1,iclass)
!
!   Rotate all the bands together.
!   NB: rotate_psi assume that s is in the small group of k. It does not
!       rotate the k point.
!
!
   DO ibnd=1,nbnd
      CALL rotate_psi(evc(1,ibnd),evcr(1,ibnd),s(1,1,irot), &
               ftau(1,irot),gk(1,irot),nl,igk,nr1,nr2,nr3,nrx1, &
               nrx2,nrx3,nrxx,ngm,npw)

   END DO
!
!   and apply S if necessary
!
   IF ( okvan ) THEN
      CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evcr ) 
      CALL s_psi( npwx, npw, nbnd, evcr, evcr )
   ENDIF
!
!  Compute the trace of the representation for each group of bands
!
   DO igroup=1,ngroup
      dim_rap=istart(igroup+1)-istart(igroup)
      trace(iclass,igroup)=(0.d0,0.d0)
      DO i=1,dim_rap
         ibnd=istart(igroup)+i-1
         trace(iclass,igroup)=trace(iclass,igroup) + &
                            ZDOTC(npw,evc(1,ibnd),1,evcr(1,ibnd),1)
      END DO
!      write(6,*) igroup, iclass, trace(iclass,igroup)
   END DO
END DO
!
CALL reduce(2*48*nbnd,trace)
!
!  And now use the character table to identify the symmetry representation
!  of each group of bands
!
WRITE(stdout,'(/,5x,"Band symmetry, ",a11," point group:",/)') gname

DO igroup=1,ngroup
   dim_rap=istart(igroup+1)-istart(igroup)
   shift=0
   DO irap=1,nclass
      times=(0.d0,0.d0)
      DO iclass=1,nclass
         times=times+trace(iclass,igroup)*char_mat(irap, &
                     which_irr(iclass))*nelem(iclass)
      ENDDO
      times=times/nsym
      IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
          (ABS(AIMAG(times)) > eps) ) THEN
            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
                      & "-->   ?")') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), dim_rap
         ibnd=istart(igroup)
         IF (rap_et(ibnd)==-1) THEN
            DO i=1,dim_rap
               ibnd=istart(igroup)+i-1
               rap_et(ibnd)=0
            END DO
         END IF
         GOTO 300
      ELSE IF (ABS(times) > eps) THEN
         ibnd=istart(igroup)+shift
         dimen=NINT(DBLE(char_mat(irap,1)))
         IF (rap_et(ibnd)==-1) THEN
            DO i=1,dimen*NINT(DBLE(times))
               ibnd=istart(igroup)+shift+i-1
               rap_et(ibnd)=irap
            ENDDO
            shift=shift+dimen*NINT(DBLE(times))
         ENDIF
         IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
            WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,&
                      & 3x,"--> ",a15)') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                              dim_rap, name_rap(irap)
         ELSE
            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
                      & "--> ",i3," ",a15)') &
              istart(igroup), istart(igroup+1)-1, &
              w1(istart(igroup)), dim_rap, NINT(DBLE(times)), name_rap(irap)
         END IF
      END IF
   END DO
300 CONTINUE
END DO
WRITE( stdout, '(/,1x,74("*"))')

DEALLOCATE(trace)
DEALLOCATE(w1)
DEALLOCATE(evcr)
DEALLOCATE(istart)
IF (okvan) DEALLOCATE(becp)

RETURN
END SUBROUTINE find_band_sym

SUBROUTINE rotate_psi(evc,evcr,s,ftau,gk,nl,igk,nr1,nr2,nr3, &
                      nrx1,nrx2,nrx3,nrxx,ngm,npw)

USE kinds,     ONLY : DP
USE constants, ONLY : tpi
IMPLICIT NONE

INTEGER :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, npw, nbnd
INTEGER :: s(3,3), ftau(3), gk(3), nl(ngm), igk(npw)

COMPLEX(DP), ALLOCATABLE :: psic(:), psir(:)
COMPLEX(DP) :: evc(npw), evcr(npw)
COMPLEX(DP) :: phase
REAL(DP) :: arg
INTEGER :: i, j, k, ri, rj, rk, ir, rir, ipol
LOGICAL :: zone_border
#if defined  (__PARA)
!
COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
COMPLEX (DP), ALLOCATABLE :: psic_collect(:)
!
#endif
!
ALLOCATE(psic(nrxx))
ALLOCATE(psir(nrxx))
!
zone_border=(gk(1).ne.0.or.gk(2).ne.0.or.gk(3).ne.0)
!
psic = ( 0.D0, 0.D0 )
!
psic(nl(igk(1:npw))) = evc(1:npw)
!
CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
!
#if defined  (__PARA)
  !
  ALLOCATE (psic_collect(nrx1*nrx2*nrx3))
  ALLOCATE (psir_collect(nrx1*nrx2*nrx3))
  !
  CALL cgather_sym( psic, psic_collect )
  !
  psir_collect=(0.d0,0.d0)
  !
  IF (zone_border) THEN
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
              rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
              arg=tpi*( (gk(1)*(i-1))/DBLE(nr1)+(gk(2)*(j-1))/DBLE(nr2)+ &
                        (gk(3)*(k-1))/DBLE(nr3) )
              phase=CMPLX(cos(arg),sin(arg))
              psir_collect(ir)=psic_collect(rir)*phase
           END DO
        END DO
     END DO
  ELSE
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
              rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
              psir_collect(ir)=psic_collect(rir)
           END DO
        END DO
     END DO
  END IF
  !
  CALL cscatter_sym( psir_collect, psir )
  !
  DEALLOCATE (psic_collect)
  DEALLOCATE (psir_collect)
  !
#else
  psir=(0.d0,0.d0)
  IF (zone_border) THEN
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
              rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
              arg=tpi*( (gk(1)*(i-1))/DBLE(nr1)+(gk(2)*(j-1))/DBLE(nr2)+ &
                        (gk(3)*(k-1))/DBLE(nr3) )
              phase=CMPLX(cos(arg),sin(arg))
              psir(ir)=psic(rir)*phase
           END DO
        END DO
     END DO
  ELSE
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
              rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
              psir(ir)=psic(rir)
           END DO
        END DO
     END DO
  END IF
  !
#endif
!
CALL cft3( psir, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
!
evcr(1:npw) = psir(nl(igk(1:npw))) 
!
DEALLOCATE(psic)
DEALLOCATE(psir)
!
RETURN
END SUBROUTINE rotate_psi

SUBROUTINE find_band_sym_so (evc,et,at,nbnd,npw,nsym,ngm,s,ftau,d_spin,gk, &
                     xk,igk,nl,nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,npwx,rap_et)
!
!   This subroutine finds the irreducible representations of the 
!   double group which give the transformation properties of the 
!   spinor wavefunctions evc. 
!   Presently it does NOT work at zone border if the space group of
!   the crystal has fractionary translations (non-symmorphic space groups).
!
!
USE io_global,          ONLY : stdout
USE kinds,              ONLY : DP
USE constants,          ONLY : rytoev
USE rap_point_group,    ONLY : code_group, nclass, gname
USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, which_irr_so, &
                               char_mat_so, name_rap_so, name_class_so,      &  
                               name_class_so1
USE rap_point_group_is, ONLY : gname_is
USE spin_orb,           ONLY : domag
USE uspp,               ONLY : vkb, nkb, okvan
USE noncollin_module,   ONLY : npol
USE becmod,             ONLY : becp_nc
IMPLICIT NONE

INTEGER :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, npw, npwx

INTEGER ::                  &
          nsym,             & 
          nbnd,             &
          rap_et(nbnd),     &
          igk(npwx),        &
          nl(ngm),          &
          ftau(3,48),       &
          gk(3,48),         &
          s(3,3,48)  

REAL(DP) ::                 &
            at(3,3),        &
            xk(3),          &
            et(nbnd)       

COMPLEX(DP) ::  &
            d_spin(2,2,48), & 
            evc(npwx*npol, nbnd)       

REAL(DP), PARAMETER :: eps=1.d-5

INTEGER ::         &
        ngroup,    &  ! number of different energy groups
        ibnd,      &
        igroup,    &
        dim_rap,   &  ! counters
        irot,      &
        irap,      &
        shift,     &
        iclass,    &
        na, i, j, ig, ipol, jpol, jrap, dimen

INTEGER, ALLOCATABLE :: istart(:)    ! point in list of energies where the
                                     ! new group starts
COMPLEX(DP) :: ZDOTC, times          ! moltiplication factors     
                                        
REAL(DP), ALLOCATABLE ::  w1(:)      ! list of energy eigenvalues in eV
COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), & ! the rotated of each wave function
                             trace(:,:)   ! the trace of the symmetry matrix
                                          ! within a given group
!
!    Divide the bands on the basis of the band degeneracy.
!
ALLOCATE(istart(nbnd+1))
ALLOCATE(w1(nbnd))
ALLOCATE(evcr(npwx*npol,nbnd))
ALLOCATE(trace(48,nbnd))
IF (okvan) ALLOCATE(becp_nc(nkb,npol,nbnd))

rap_et=-1
w1=et*rytoev
!
!  divide the energies in groups of degenerate eigenvalues. Two eigenvalues
!  are assumed to be degenerate if their difference is less than 0.0001 eV.
!
ngroup=1
istart(ngroup)=1
DO ibnd=2,nbnd
   IF (ABS(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
      ngroup=ngroup+1
      istart(ngroup)=ibnd
   END IF
END DO
istart(ngroup+1)=nbnd+1
!
!  Find the character of one symmetry operation per class
!
trace=(0.d0,0.d0)
DO iclass=1,nclass
   irot=elem_so(1,iclass)
!
!   Rotate all the bands together.
!   NB: rotate_psi assumes that s is in the small group of k. It does not
!       rotate the k point.
!
   DO ibnd=1,nbnd
      CALL rotate_psi_so(evc(1,ibnd),evcr(1,ibnd),s(1,1,irot),        &
               ftau(1,irot),d_spin(1,1,irot),has_e(1,iclass),gk(1,irot), &
               nl,igk,npol,nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,ngm,npw,npwx)
   END DO
!
!   and apply S in the US case.
!
   IF ( okvan ) THEN
      CALL ccalbec_nc( nkb, npwx, npw, npol, nbnd, becp_nc, vkb, evcr ) 
      CALL s_psi_nc( npwx, npw, nbnd, evcr, evcr )
   ENDIF
!
!  Compute the trace of the representation for each group of bands
!
   DO igroup=1,ngroup
      dim_rap=istart(igroup+1)-istart(igroup)
      DO i=1,dim_rap
         ibnd=istart(igroup)+i-1
         trace(iclass,igroup)=trace(iclass,igroup) +            &
                ZDOTC(2*npwx,evc(1,ibnd),1,evcr(1,ibnd),1)
      END DO
!      write(6,*) igroup, iclass, dim_rap, trace(iclass,igroup)
   END DO
END DO
!
CALL reduce(2*48*nbnd,trace)
!
!DO iclass=1,nclass
!   write(6,'(i5,3(2f10.5,3x))') iclass,trace(iclass,1),trace(iclass,2), &
!                                       trace(iclass,3)
!ENDDO
!
!  And now use the character table to identify the symmetry representation
!  of each group of bands
!
IF (domag) THEN
   WRITE(stdout,'(/,5x,"Band symmetry, ",a11," [",a11, &
          & "] magnetic double point group,")') gname, gname_is
   WRITE(stdout,'(5x,"using ",a11,/)') gname_is
ELSE
   WRITE(stdout,'(/,5x,"Band symmetry, ",a11," double point group:",/)') gname
ENDIF

DO igroup=1,ngroup
   dim_rap=istart(igroup+1)-istart(igroup)
   shift=0
   DO irap=1,nrap
      times=(0.d0,0.d0)
      DO iclass=1,nclass
         times=times+CONJG(trace(iclass,igroup))*char_mat_so(irap, &
                     which_irr_so(iclass))*DBLE(nelem_so(iclass))
      ENDDO
      times=times/2/nsym
      IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
          (ABS(AIMAG(times)) > eps) ) THEN
            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
                      & "-->   ?")') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), dim_rap
         ibnd=istart(igroup)
         IF (rap_et(ibnd)==-1) THEN
            DO i=1,dim_rap
               ibnd=istart(igroup)+i-1
               rap_et(ibnd)=0
            END DO
         END IF
         GOTO 300
      END IF
      IF (ABS(times) > eps) THEN
         dimen=NINT(DBLE(char_mat_so(irap,1)))
         ibnd=istart(igroup) + shift
         IF (rap_et(ibnd)==-1) THEN
            DO i=1,dimen*NINT(DBLE(times))
               ibnd=istart(igroup)+shift+i-1
               rap_et(ibnd)=irap
            END DO
            shift=shift+dimen*NINT(DBLE(times))
         END IF
         IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
            WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
                      & "--> ",a15)') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                             dim_rap, name_rap_so(irap)
         ELSE
            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,&
                      & 3x,"--> ",i3," ",a15)') &
              istart(igroup), istart(igroup+1)-1, &
              w1(istart(igroup)), dim_rap, NINT(DBLE(times)), name_rap_so(irap)
         END IF
      END IF
   END DO
300 CONTINUE
END DO
WRITE( stdout, '(/,1x,74("*"))')

DEALLOCATE(trace)
DEALLOCATE(w1)
DEALLOCATE(evcr)
DEALLOCATE(istart)
IF (okvan) DEALLOCATE(becp_nc)
RETURN
END SUBROUTINE find_band_sym_so

SUBROUTINE rotate_psi_so(evc_nc,evcr,s,ftau,d_spin,has_e,gk,nl,igk,npol, &
                         nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,ngm,npw,npwx)
!
!  This subroutine rotates a spinor wavefunction according to the symmetry
!  s. d_spin contains the 2x2 rotation matrix in the spin space.
!  has_e=-1 means that also a 360 degrees rotation is applied in spin space.
!
USE kinds,     ONLY : DP
USE constants, ONLY : tpi
IMPLICIT NONE

INTEGER :: npol, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, npw, nbnd, npwx
INTEGER :: s(3,3), ftau(3), gk(3), nl(ngm), igk(npw), has_e

COMPLEX(DP), ALLOCATABLE :: psic(:,:), psir(:,:), evcr_save(:,:)
COMPLEX(DP) :: evc_nc(npwx,2), evcr(npwx,2), d_spin(2,2)
COMPLEX(DP) :: phase
REAL(DP) :: arg, sum
INTEGER :: i, j, k, ri, rj, rk, ir, rir, ipol, jpol
LOGICAL :: zone_border
!
#if defined  (__PARA)
!
COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
COMPLEX (DP), ALLOCATABLE :: psic_collect(:)
!
ALLOCATE (psic_collect(nrx1*nrx2*nrx3))
ALLOCATE (psir_collect(nrx1*nrx2*nrx3))
#endif
!
ALLOCATE(psic(nrxx,npol))
ALLOCATE(psir(nrxx,npol))
ALLOCATE(evcr_save(npwx,npol))
!
zone_border=(gk(1).ne.0.or.gk(2).ne.0.or.gk(3).ne.0)
!
psic = ( 0.D0, 0.D0 )
psir = ( 0.D0, 0.D0 )
!
DO ipol=1,npol
   !
   psic(nl(igk(1:npw)),ipol) = evc_nc(1:npw,ipol)
   CALL cft3( psic(1,ipol), nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
   !
#if defined  (__PARA)
   !
   !
   CALL cgather_sym( psic(1,ipol), psic_collect )
   !
   psir_collect=(0.d0,0.d0)
   !
   IF (zone_border) THEN
      DO k = 1, nr3
         DO j = 1, nr2
            DO i = 1, nr1
               CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
               ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
               rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
               arg=tpi*( (gk(1)*(i-1))/DBLE(nr1)+(gk(2)*(j-1))/DBLE(nr2)+ &
                         (gk(3)*(k-1))/DBLE(nr3) )
               phase=CMPLX(cos(arg),sin(arg))
               psir_collect(ir)=psic_collect(rir)*phase
            END DO
         END DO
      END DO
   ELSE
      DO k = 1, nr3
         DO j = 1, nr2
            DO i = 1, nr1
               CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
               ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
               rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
               psir_collect(ir)=psic_collect(rir)
            END DO
         END DO
      END DO
   END IF
   !
   CALL cscatter_sym( psir_collect, psir(1,ipol) )
   !
#else
   IF (zone_border) THEN
      DO k = 1, nr3
         DO j = 1, nr2
            DO i = 1, nr1
               CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
               ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
               rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
               arg=tpi*( (gk(1)*(i-1))/DBLE(nr1)+(gk(2)*(j-1))/DBLE(nr2)+ &
                         (gk(3)*(k-1))/DBLE(nr3) )
               phase=CMPLX(COS(arg),SIN(arg))
               psir(ir,ipol)=psic(rir,ipol)*phase
            END DO
         END DO
      END DO
   ELSE
      DO k = 1, nr3
         DO j = 1, nr2
            DO i = 1, nr1
               CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
               ir=i+(j-1)*nrx1+(k-1)*nrx1*nrx2 
               rir=ri+(rj-1)*nrx1+(rk-1)*nrx1*nrx2 
               psir(ir,ipol)=psic(rir,ipol)
            END DO
         END DO
      END DO
   END IF
   !
#endif
   !
   CALL cft3( psir(1,ipol), nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
   !
   evcr_save(1:npw,ipol) = psir(nl(igk(1:npw)),ipol) 
   !
ENDDO
evcr=(0.d0,0.d0)
DO ipol=1,npol
   DO jpol=1,npol
      evcr(:,ipol)=evcr(:,ipol)+CONJG(d_spin(jpol,ipol))*evcr_save(:,jpol)
   END DO
END DO
IF (has_e==-1) evcr=-evcr
!
DEALLOCATE(evcr_save)
DEALLOCATE(psic)
DEALLOCATE(psir)
#if defined (__PARA)
   DEALLOCATE (psic_collect)
   DEALLOCATE (psir_collect)
#endif
RETURN
END SUBROUTINE rotate_psi_so
