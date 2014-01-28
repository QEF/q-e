!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE sym_band(filband, spin_component, firstk, lastk)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba2, at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, nl, g
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin, ecutwfc
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE symm_base,            ONLY : s, ftau, nsym, t_rev, sname
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
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j, irot, iclass, ig, ibnd
  INTEGER :: spin_component, nks1, nks2, firstk, lastk
  INTEGER :: nks1tot, nks2tot
  INTEGER :: iunout, igroup, irap, dim_rap, ios
  INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), sk_is(3,3,48), &
       gk_is(3,48), t_revk(48), nsymk, isym, ipol, jpol
  LOGICAL :: is_complex, is_complex_so, is_symmorphic, is_high_sym, search_sym
  REAL(DP), PARAMETER :: accuracy=1.d-4
  COMPLEX(DP) :: d_spink(2,2,48), d_spin_is(2,2,48), zdotc
  COMPLEX(DP),ALLOCATABLE :: times(:,:,:)
  REAL(DP) :: dxk(3), dkmod, dkmod_save
  INTEGER, ALLOCATABLE :: rap_et(:,:), code_group_k(:)
  INTEGER, ALLOCATABLE :: ngroup(:), istart(:,:)
  CHARACTER(len=11) :: group_name
  CHARACTER(len=45) :: snamek(48)
  CHARACTER (len=256) :: filband, namefile
  !
  IF (spin_component/=1.and.nspin/=2) &
       CALL errore('sym_band','incorrect spin_component',1)
  IF (spin_component<1.or.spin_component>2) &
       CALL errore('sym_band','incorrect lsda spin_component',1)

  ALLOCATE(rap_et(nbnd,nkstot))
  ALLOCATE(code_group_k(nkstot))
  ALLOCATE(times(nbnd,24,nkstot))
  ALLOCATE(ngroup(nkstot))
  ALLOCATE(istart(nbnd+1,nkstot))
  code_group_k=0
  rap_et=-1
  times=(0.0_DP,0.0_DP)
  CALL find_nks1nks2(firstk,lastk,nks1tot,nks1,nks2tot,nks2,spin_component)

  ios=0
  IF ( ionode ) THEN
     iunout=58
     namefile=trim(filband)//".rap"
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0) CALL errore ('sym_band', 'Opening filband file', abs (ios) )

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
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     ! Find the small group of k
     !
     CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, nsym, sk, ftauk, &
          gk, t_revk, snamek, nsymk)
     !
     !  character of the irreducible representations
     !
     CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
          sk_is,d_spin_is,gk_is, &
          is_symmorphic,search_sym)
     code_group_k(ik)=code_group
     !
     IF (.not.search_sym) THEN
        rap_et(:,ik)=-1
        GOTO 100
     ENDIF
     !
     !  Find the symmetry of each state
     !
     IF (noncolin) THEN
        IF (domag) THEN
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsym_is, &
                ngm,sk_is,ftau_is,d_spin_is,gk_is,xk(1,ik),igk,nl,dfftp%nr1,dfftp%nr2,&
                dfftp%nr3,dfftp%nr1x,dfftp%nr2x,dfftp%nr3x,dfftp%nnr,npwx,rap_et(1,ik),times(1,1,ik), &
                ngroup(ik),istart(1,ik),accuracy)
        ELSE
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsymk,ngm, &
                sk,ftauk,d_spink,gk,xk(1,ik),igk,nl,dfftp%nr1,dfftp%nr2,dfftp%nr3,dfftp%nr1x, &
                dfftp%nr2x,dfftp%nr3x,dfftp%nnr,npwx,rap_et(1,ik),times(1,1,ik),ngroup(ik),&
                istart(1,ik),accuracy)
        ENDIF
     ELSE
        CALL find_band_sym (evc, et(1,ik), at, nbnd, npw, nsymk, ngm, &
             sk, ftauk, gk, xk(1,ik), igk, nl, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, &
             dfftp%nr2x, dfftp%nr3x, dfftp%nnr, npwx, rap_et(1,ik), times(1,1,ik), ngroup(ik),&
             istart(1,ik),accuracy)
     ENDIF

100  CONTINUE
  ENDDO

#ifdef __MPI
  !
  !  Only the symmetry of a set of k points is calculated by this
  !  processor with pool. Here we collect the results into ionode
  !
  CALL ipoolrecover(code_group_k,1,nkstot,nks)
  CALL ipoolrecover(rap_et,nbnd,nkstot,nks)
  CALL poolrecover(times,2*24*nbnd,nkstot,nks)
  CALL ipoolrecover(ngroup,1,nkstot,nks)
  CALL ipoolrecover(istart,nbnd+1,nkstot,nks)
#endif
  IF (ionode) THEN
     is_high_sym=.false.
     DO ik=nks1tot, nks2tot
        CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, &
             nsym, sk, ftauk, gk, t_revk, snamek, nsymk)
        CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
             sk_is,d_spin_is,gk_is, &
             is_symmorphic,search_sym)
        IF (code_group_k(ik) /= code_group) &
             CALL errore('sym_band','problem with code_group',1)
        WRITE(stdout, '(/,1x,74("*"))')
        WRITE(stdout, '(/,20x,"xk=(",2(f10.5,","),f10.5,"  )")') &
             xk(1,ik), xk(2,ik), xk(3,ik)
        IF (.not.search_sym) THEN
           WRITE(stdout,'(/,5x,"zone border point and non-symmorphic group ")')
           WRITE(stdout,'(5x,"symmetry decomposition not available")')
           WRITE(stdout, '(/,1x,74("*"))')
        ENDIF
        IF (ik == nks1tot) THEN
           WRITE (iunout, '(" &plot_rap nbnd_rap=",i4,", nks_rap=",i4," /")') &
                nbnd, nks2tot-nks1tot+1
           IF (search_sym) CALL write_group_info(.true.)
           is_high_sym=.true.
           dxk(:) = xk(:,nks1tot+1) - xk(:,nks1tot)
           dkmod_save = sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
        ELSE
           IF (code_group_k(ik)/=code_group_k(ik-1).and.search_sym) &
                CALL write_group_info(.true.)
!
!    When the symmetry changes the point must be considered a high
!    symmetry point. If the previous point was also high_symmetry, there
!    are two possibilities. The two points are distant and in this case
!    both of them must be considered high symmetry. If they are close only
!    the first point is a high symmetry point. First compute the distance
!
           dxk(:) = xk(:,ik) - xk(:,ik-1)
           dkmod= sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
           IF (dkmod < 5.0_DP * dkmod_save) THEN
!
!    In this case the two points are considered close
!
              is_high_sym= ((code_group_k(ik)/=code_group_k(ik-1)) &
                 .and..not.is_high_sym) 
              IF (dkmod > 1.d-3) dkmod_save=dkmod
           ELSE
!
!    Points are distant. They are all high symmetry
!
              is_high_sym= .TRUE.
           ENDIF
        ENDIF
        WRITE (iunout, '(10x,3f10.6,l5)') xk(1,ik),xk(2,ik),xk(3,ik), &
             is_high_sym
        WRITE (iunout, '(10i8)') (rap_et(ibnd,ik), ibnd=1,nbnd)
        IF (.not.search_sym) CYCLE
        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11," [",a11, &
                   & "] magnetic double point group,")') gname, gname_is
              WRITE(stdout,'(5x,"using ",a11,/)') gname_is
           ELSE
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11,&
                   & " double point group:",/)') gname
           ENDIF
        ELSE
           WRITE(stdout,'(/,5x,"Band symmetry, ",a11," point group:",/)') &
                group_name(code_group_k(ik))
        ENDIF

        DO igroup=1,ngroup(ik)
           dim_rap=istart(igroup+1,ik)-istart(igroup,ik)
           DO irap=1,nclass
              IF (noncolin) THEN
                 IF ((abs(nint(dble(times(igroup,irap,ik)))-  &
                      dble(times(igroup,irap,ik))) > accuracy).or. &
                      (abs(aimag(times(igroup,irap,ik))) > accuracy) ) THEN
                    WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,&
                         &"eV",3x,i3,3x, "-->   ?")') &
                         istart(igroup,ik), istart(igroup+1,ik)-1, &
                         et(istart(igroup,ik),ik)*rytoev, dim_rap
                    exit
                 ELSEIF (abs(times(igroup,irap,ik)) > accuracy) THEN
                    IF (abs(nint(dble(times(igroup,irap,ik))-1.d0)) < &
                         accuracy) THEN
                       WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, &
                            dim_rap, name_rap_so(irap)
                    ELSE
                       WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",i3," ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, dim_rap, &
                            nint(dble(times(igroup,irap,ik))), name_rap_so(irap)
                    ENDIF
                 ENDIF
              ELSE
                 IF ((abs(nint(dble(times(igroup,irap,ik)))-  &
                      dble(times(igroup,irap,ik))) > accuracy).or. &
                      (abs(aimag(times(igroup,irap,ik))) > accuracy) ) THEN
                    WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,&
                         &"eV",3x,i3,3x, "-->   ?")') &
                         istart(igroup,ik), istart(igroup+1,ik)-1, &
                         et(istart(igroup,ik),ik)*rytoev, dim_rap
                    exit
                 ELSEIF (abs(times(igroup,irap,ik)) > accuracy) THEN
                    IF (abs(nint(dble(times(igroup,irap,ik))-1.d0)) < &
                         accuracy) THEN
                       WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, &
                            dim_rap, name_rap(irap)
                    ELSE
                       WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",i3," ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, dim_rap, &
                            nint(dble(times(igroup,irap,ik))), name_rap(irap)
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
        WRITE( stdout, '(/,1x,74("*"))')
     ENDDO
     CLOSE(iunout)
  ENDIF

  !
  DEALLOCATE(times)
  DEALLOCATE(code_group_k)
  DEALLOCATE(rap_et)
  DEALLOCATE(ngroup)
  DEALLOCATE(istart)
  !
  RETURN
END SUBROUTINE sym_band
!
SUBROUTINE find_band_sym (evc,et,at,nbnd,npw,nsym,ngm,s,ftau,gk, &
     xk,igk,nl,nr1,nr2,nr3,nr1x,nr2x,nr3x,nrxx,npwx, &
     rap_et,times,ngroup,istart,accuracy)
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
  USE becmod,          ONLY : bec_type, becp, calbec, &
       allocate_bec_type, deallocate_bec_type
  USE mp_global,       ONLY : intra_pool_comm
  USE mp,              ONLY : mp_sum

  IMPLICIT NONE

  REAL(DP), INTENT(in) :: accuracy

  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, ngm, npw, npwx

  INTEGER ::                  &
       nsym,             &
       nbnd,             &
       rap_et(nbnd),     &
       igk(npwx),        &
       nl(ngm),          &
       ftau(3,48),       &
       gk(3,48),         &
       s(3,3,48),        &
       ngroup,           &  ! number of different frequencies groups
       istart(nbnd+1)

  REAL(DP) ::                 &
       at(3,3),        &
       xk(3),          &
       et(nbnd)

  COMPLEX(DP) ::  &
       times(nbnd,24), &
       evc(npwx, nbnd)

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::         &
       ibnd,      &
       igroup,    &
       dim_rap,   &
       irot,      &
       irap,      &
       iclass,    &
       shift,     &
       na, i, j, ig, dimen

  COMPLEX(DP) :: zdotc

  REAL(DP), ALLOCATABLE ::  w1(:)
  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), trace(:,:)
  !
  !    Divide the bands on the basis of the band degeneracy.
  !
  ALLOCATE(w1(nbnd))
  ALLOCATE(evcr(npwx,nbnd))
  ALLOCATE(trace(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )

  rap_et=-1
  w1=et*rytoev

  ngroup=1
  istart(ngroup)=1
  DO ibnd=2,nbnd
     IF (abs(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd
     ENDIF
  ENDDO
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
             ftau(1,irot),gk(1,irot),nl,igk,nr1,nr2,nr3,nr1x, &
             nr2x,nr3x,nrxx,ngm,npw)

     ENDDO
     !
     !   and apply S if necessary
     !
     IF ( okvan ) THEN
        CALL calbec( npw, vkb, evcr, becp )
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
                zdotc(npw,evc(1,ibnd),1,evcr(1,ibnd),1)
        ENDDO
        !      write(6,*) igroup, iclass, trace(iclass,igroup)
     ENDDO
  ENDDO
  !
  CALL mp_sum( trace, intra_pool_comm )

  !DO iclass=1,nclass
  !   write(6,'(i5,3(2f11.8,1x))') iclass,trace(iclass,4),trace(iclass,5), &
  !                                       trace(iclass,6)
  !ENDDO

  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of bands
  !
  !WRITE(stdout,'(/,5x,"Band symmetry, ",a11," point group:",/)') gname

  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     shift=0
     DO irap=1,nclass
        times(igroup,irap)=(0.d0,0.d0)
        DO iclass=1,nclass
           times(igroup,irap)=times(igroup,irap) &
                +trace(iclass,igroup)*char_mat(irap,which_irr(iclass))&
                *nelem(iclass)
        ENDDO
        times(igroup,irap)=times(igroup,irap)/nsym
        IF ((abs(nint(dble(times(igroup,irap)))-dble(times(igroup,irap))) &
             > accuracy).or. (abs(aimag(times(igroup,irap))) > eps) ) THEN
           !            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
           !                      & "-->   ?")') &
           !              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), dim_rap
           ibnd=istart(igroup)
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap
                 ibnd=istart(igroup)+i-1
                 rap_et(ibnd)=0
              ENDDO
           ENDIF
           GOTO 300
        ELSEIF (abs(times(igroup,irap)) > accuracy) THEN
           ibnd=istart(igroup)+shift
           dimen=nint(dble(char_mat(irap,1)))
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dimen*nint(dble(times(igroup,irap)))
                 ibnd=istart(igroup)+shift+i-1
                 rap_et(ibnd)=irap
              ENDDO
              shift=shift+dimen*nint(dble(times(igroup,irap)))
           ENDIF
           !         IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
           !            WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,&
           !                      & 3x,"--> ",a15)') &
           !              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
           !                              dim_rap, name_rap(irap)
           !         ELSE
           !            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
           !                      & "--> ",i3," ",a15)') &
           !              istart(igroup), istart(igroup+1)-1, &
           !              w1(istart(igroup)), dim_rap, NINT(DBLE(times)), name_rap(irap)
           !         END IF
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO
  !WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(w1)
  DEALLOCATE(evcr)
  IF (okvan) CALL deallocate_bec_type (becp)

  RETURN
END SUBROUTINE find_band_sym

SUBROUTINE rotate_psi(evc,evcr,s,ftau,gk,nl,igk,nr1,nr2,nr3, &
     nr1x,nr2x,nr3x,nrxx,ngm,npw)

  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE fft_base,  ONLY : cgather_sym, cscatter_sym, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft

  IMPLICIT NONE

  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, ngm, npw, nbnd
  INTEGER :: s(3,3), ftau(3), gk(3), nl(ngm), igk(npw)

  COMPLEX(DP), ALLOCATABLE :: psic(:), psir(:)
  COMPLEX(DP) :: evc(npw), evcr(npw)
  COMPLEX(DP) :: phase
  REAL(DP) :: arg
  INTEGER :: i, j, k, ri, rj, rk, ir, rir, ipol
  LOGICAL :: zone_border
#if defined  (__MPI)
  !
  COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
  COMPLEX (DP), ALLOCATABLE :: psic_collect(:)
  !
#endif
  !
  ALLOCATE(psic(nrxx))
  ALLOCATE(psir(nrxx))
  !
  zone_border=(gk(1)/=0.or.gk(2)/=0.or.gk(3)/=0)
  !
  psic = ( 0.D0, 0.D0 )
  !
  psic(nl(igk(1:npw))) = evc(1:npw)
  !
  CALL invfft ('Dense', psic, dfftp)
  !
#if defined  (__MPI)
  !
  ALLOCATE (psic_collect(nr1x*nr2x*nr3x))
  ALLOCATE (psir_collect(nr1x*nr2x*nr3x))
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
              ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
              rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
              arg=tpi*( (gk(1)*(i-1))/dble(nr1)+(gk(2)*(j-1))/dble(nr2)+ &
                   (gk(3)*(k-1))/dble(nr3) )
              phase=cmplx(cos(arg),sin(arg),kind=DP)
              psir_collect(ir)=psic_collect(rir)*phase
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
              rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
              psir_collect(ir)=psic_collect(rir)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
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
              ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
              rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
              arg=tpi*( (gk(1)*(i-1))/dble(nr1)+(gk(2)*(j-1))/dble(nr2)+ &
                   (gk(3)*(k-1))/dble(nr3) )
              phase=cmplx(cos(arg),sin(arg),kind=DP)
              psir(ir)=psic(rir)*phase
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
              rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
              psir(ir)=psic(rir)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
#endif
  !
  CALL fwfft ('Dense', psir, dfftp)
  !
  evcr(1:npw) = psir(nl(igk(1:npw)))
  !
  DEALLOCATE(psic)
  DEALLOCATE(psir)
  !
  RETURN
END SUBROUTINE rotate_psi

SUBROUTINE find_band_sym_so (evc,et,at,nbnd,npw,nsym,ngm,s,ftau,d_spin,gk, &
     xk,igk,nl,nr1,nr2,nr3,nr1x,nr2x,nr3x,nrxx,npwx, &
     rap_et,times,ngroup,istart,accuracy)

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
  USE becmod,             ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE mp_global,          ONLY : intra_pool_comm
  USE mp,                 ONLY : mp_sum

  IMPLICIT NONE

  REAL(DP), INTENT(in) :: accuracy

  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, ngm, npw, npwx

  INTEGER ::                  &
       nsym,             &
       nbnd,             &
       ngroup,           &
       istart(nbnd+1),   &
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
       times(nbnd,24), &
       d_spin(2,2,48), &
       evc(npwx*npol, nbnd)

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::         &
       ibnd,      &
       igroup,    &
       dim_rap,   &  ! counters
       irot,      &
       irap,      &
       shift,     &
       iclass,    &
       na, i, j, ig, ipol, jpol, jrap, dimen

  COMPLEX(DP) :: zdotc          ! moltiplication factors

  REAL(DP), ALLOCATABLE ::  w1(:)      ! list of energy eigenvalues in eV
  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), & ! the rotated of each wave function
       trace(:,:)   ! the trace of the symmetry matrix
  ! within a given group
  !
  !    Divide the bands on the basis of the band degeneracy.
  !
  ALLOCATE(w1(nbnd))
  ALLOCATE(evcr(npwx*npol,nbnd))
  ALLOCATE(trace(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )

  rap_et=-1
  w1=et*rytoev
  !
  !  divide the energies in groups of degenerate eigenvalues. Two eigenvalues
  !  are assumed to be degenerate if their difference is less than 0.0001 eV.
  !
  ngroup=1
  istart(ngroup)=1
  DO ibnd=2,nbnd
     IF (abs(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd
     ENDIF
  ENDDO
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
             nl,igk,npol,nr1,nr2,nr3,nr1x,nr2x,nr3x,nrxx,ngm,npw,npwx)
     ENDDO
     !
     !   and apply S in the US case.
     !
     IF ( okvan ) THEN
        CALL calbec( npw, vkb, evcr, becp )
        CALL s_psi( npwx, npw, nbnd, evcr, evcr )
     ENDIF
     !
     !  Compute the trace of the representation for each group of bands
     !
     DO igroup=1,ngroup
        dim_rap=istart(igroup+1)-istart(igroup)
        DO i=1,dim_rap
           ibnd=istart(igroup)+i-1
           trace(iclass,igroup)=trace(iclass,igroup) +            &
                zdotc(2*npwx,evc(1,ibnd),1,evcr(1,ibnd),1)
        ENDDO
        !      write(6,*) igroup, iclass, dim_rap, trace(iclass,igroup)
     ENDDO
  ENDDO
  !
  CALL mp_sum(trace,intra_pool_comm)
  !
  !DO iclass=1,nclass
  !   write(6,'(i5,3(2f11.8,1x))') iclass,trace(iclass,1),trace(iclass,2), &
  !                                       trace(iclass,3)
  !ENDDO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of bands
  !
  !IF (domag) THEN
  !   WRITE(stdout,'(/,5x,"Band symmetry, ",a11," [",a11, &
  !          & "] magnetic double point group,")') gname, gname_is
  !   WRITE(stdout,'(5x,"using ",a11,/)') gname_is
  !ELSE
  !   WRITE(stdout,'(/,5x,"Band symmetry, ",a11," double point group:",/)') gname
  !ENDIF

  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     shift=0
     DO irap=1,nrap
        times(igroup,irap)=(0.d0,0.d0)
        DO iclass=1,nclass
           times(igroup,irap)=times(igroup,irap) &
                +conjg(trace(iclass,igroup))*char_mat_so(irap, &
                which_irr_so(iclass))*dble(nelem_so(iclass))
        ENDDO
        times(igroup,irap)=times(igroup,irap)/2/nsym
        IF ((abs(nint(dble(times(igroup,irap)))-dble(times(igroup,irap)))&
             > accuracy).or. (abs(aimag(times(igroup,irap))) > accuracy) ) THEN
           !            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
           !                      & "-->   ?")') &
           !              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), dim_rap
           ibnd=istart(igroup)
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap
                 ibnd=istart(igroup)+i-1
                 rap_et(ibnd)=0
              ENDDO
           ENDIF
           GOTO 300
        ENDIF
        IF (abs(times(igroup,irap)) > accuracy) THEN
           dimen=nint(dble(char_mat_so(irap,1)))
           ibnd=istart(igroup) + shift
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dimen*nint(dble(times(igroup,irap)))
                 ibnd=istart(igroup)+shift+i-1
                 rap_et(ibnd)=irap
              ENDDO
              shift=shift+dimen*nint(dble(times(igroup,irap)))
           ENDIF
           !         IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
           !            WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,3x,&
           !                      & "--> ",a15)') &
           !              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
           !                             dim_rap, name_rap_so(irap)
           !         ELSE
           !            WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,"eV",3x,i3,&
           !                      & 3x,"--> ",i3," ",a15)') &
           !              istart(igroup), istart(igroup+1)-1, &
           !              w1(istart(igroup)), dim_rap, NINT(DBLE(times)), name_rap_so(irap)
           !         END IF
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO
  !WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(w1)
  DEALLOCATE(evcr)
  IF (okvan) CALL deallocate_bec_type ( becp )
  RETURN
END SUBROUTINE find_band_sym_so

SUBROUTINE rotate_psi_so(evc_nc,evcr,s,ftau,d_spin,has_e,gk,nl,igk,npol, &
     nr1,nr2,nr3,nr1x,nr2x,nr3x,nrxx,ngm,npw,npwx)
  !
  !  This subroutine rotates a spinor wavefunction according to the symmetry
  !  s. d_spin contains the 2x2 rotation matrix in the spin space.
  !  has_e=-1 means that also a 360 degrees rotation is applied in spin space.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE fft_base,  ONLY : cgather_sym, cscatter_sym, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft

  IMPLICIT NONE

  INTEGER :: npol, nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, ngm, npw, nbnd, npwx
  INTEGER :: s(3,3), ftau(3), gk(3), nl(ngm), igk(npw), has_e

  COMPLEX(DP), ALLOCATABLE :: psic(:,:), psir(:,:), evcr_save(:,:)
  COMPLEX(DP) :: evc_nc(npwx,2), evcr(npwx,2), d_spin(2,2)
  COMPLEX(DP) :: phase
  REAL(DP) :: arg, sum
  INTEGER :: i, j, k, ri, rj, rk, ir, rir, ipol, jpol
  LOGICAL :: zone_border
  !
#if defined  (__MPI)
  !
  COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
  COMPLEX (DP), ALLOCATABLE :: psic_collect(:)
  !
  ALLOCATE (psic_collect(nr1x*nr2x*nr3x))
  ALLOCATE (psir_collect(nr1x*nr2x*nr3x))
#endif
  !
  ALLOCATE(psic(nrxx,npol))
  ALLOCATE(psir(nrxx,npol))
  ALLOCATE(evcr_save(npwx,npol))
  !
  zone_border=(gk(1)/=0.or.gk(2)/=0.or.gk(3)/=0)
  !
  psic = ( 0.D0, 0.D0 )
  psir = ( 0.D0, 0.D0 )
  !
  DO ipol=1,npol
     !
     psic(nl(igk(1:npw)),ipol) = evc_nc(1:npw,ipol)
     CALL invfft ('Dense', psic(:,ipol), dfftp)
     !
#if defined  (__MPI)
     !
     !
     CALL cgather_sym( psic(:,ipol), psic_collect )
     !
     psir_collect=(0.d0,0.d0)
     !
     IF (zone_border) THEN
        DO k = 1, nr3
           DO j = 1, nr2
              DO i = 1, nr1
                 CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
                 ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
                 rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
                 arg=tpi*( (gk(1)*(i-1))/dble(nr1)+(gk(2)*(j-1))/dble(nr2)+ &
                      (gk(3)*(k-1))/dble(nr3) )
                 phase=cmplx(cos(arg),sin(arg),kind=DP)
                 psir_collect(ir)=psic_collect(rir)*phase
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO k = 1, nr3
           DO j = 1, nr2
              DO i = 1, nr1
                 CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
                 ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
                 rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
                 psir_collect(ir)=psic_collect(rir)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
     CALL cscatter_sym( psir_collect, psir(:,ipol) )
     !
#else
     IF (zone_border) THEN
        DO k = 1, nr3
           DO j = 1, nr2
              DO i = 1, nr1
                 CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
                 ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
                 rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
                 arg=tpi*( (gk(1)*(i-1))/dble(nr1)+(gk(2)*(j-1))/dble(nr2)+ &
                      (gk(3)*(k-1))/dble(nr3) )
                 phase=cmplx(cos(arg),sin(arg),kind=DP)
                 psir(ir,ipol)=psic(rir,ipol)*phase
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO k = 1, nr3
           DO j = 1, nr2
              DO i = 1, nr1
                 CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
                 ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
                 rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
                 psir(ir,ipol)=psic(rir,ipol)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
#endif
     !
     CALL fwfft ('Dense', psir(:,ipol), dfftp)
     !
     evcr_save(1:npw,ipol) = psir(nl(igk(1:npw)),ipol)
     !
  ENDDO
  evcr=(0.d0,0.d0)
  DO ipol=1,npol
     DO jpol=1,npol
        evcr(:,ipol)=evcr(:,ipol)+conjg(d_spin(jpol,ipol))*evcr_save(:,jpol)
     ENDDO
  ENDDO
  IF (has_e==-1) evcr=-evcr
  !
  DEALLOCATE(evcr_save)
  DEALLOCATE(psic)
  DEALLOCATE(psir)
#if defined (__MPI)
  DEALLOCATE (psic_collect)
  DEALLOCATE (psir_collect)
#endif
  RETURN
END SUBROUTINE rotate_psi_so

SUBROUTINE find_nks1nks2(firstk,lastk,nks1tot,nks1,nks2tot,nks2,spin_component)
  !
  !  This routine selects the first (nks1) and last (nks2) k point calculated
  !  by the  current pool.
  !
  USE lsda_mod, ONLY : nspin
  USE klist,  ONLY : nks, nkstot
  USE mp_global, ONLY : my_pool_id, npool, kunit

  IMPLICIT NONE
  INTEGER, INTENT(out) :: nks1tot,nks1,nks2tot,nks2
  INTEGER, INTENT(in) :: firstk, lastk, spin_component
  INTEGER :: nbase, rest

  IF (nspin==1.or.nspin==4) THEN
     nks1tot=max(1,firstk)
     nks2tot=min(nkstot, lastk)
  ELSEIF (nspin==2) THEN
     IF (spin_component == 1) THEN
        nks1tot=max(1,firstk)
        nks2tot=min(nkstot/2,lastk)
     ELSEIF (spin_component==2) THEN
        nks1tot=nkstot/2 + max(1,firstk)
        nks2tot=nkstot/2 + min(nkstot/2,lastk)
     ENDIF
  ENDIF
  IF (nks1tot>nks2tot) CALL errore('findnks1nks2','wrong nks1tot or nks2tot',1)

#ifdef __MPI
  nks  = kunit * ( nkstot / kunit / npool )
  rest = ( nkstot - nks * npool ) / kunit
  IF ( ( my_pool_id + 1 ) <= rest ) nks = nks + kunit
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit

  nks1=max(1,nks1tot-nbase)
  IF (nks1>nks) nks1=nks+1
  nks2=min(nks,nks2tot-nbase)
  IF (nks2<1) nks2=nks1-1
#else
  nks1=nks1tot
  nks2=nks2tot
#endif

END SUBROUTINE find_nks1nks2

SUBROUTINE find_info_group(nsym,s,t_rev,ftau,d_spink,gk,sname,  &
     s_is,d_spin_is,gk_is, &
     is_symmorphic,search_sym)
  !
  ! This routine receives as input a point group and sets the corresponding
  ! variables for the description of the classes and of the irreducible
  ! representations. It sets also the group name and code.
  ! In the magnetic case it selects the invariat subgroup.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg
  USE noncollin_module,     ONLY : noncolin
  USE spin_orb,             ONLY : domag
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, &
       name_class_so, d_spin, name_class_so1
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ftau_is, gname_is, &
       sname_is, code_group_is

  IMPLICIT NONE

  INTEGER, INTENT(in) :: nsym,        & ! dimension of the group
       s(3,3,48),   & ! rotation matrices
       t_rev(48),   & ! if time reversal is need
       ftau(3,48),  & ! fractionary translation
       gk(3,48)

  INTEGER, INTENT(out) :: s_is(3,3,48),   & ! rotation matrices
       gk_is(3,48)

  COMPLEX(DP),INTENT(out)   :: d_spink(2,2,48),  & ! rotation in spin space
       d_spin_is(2,2,48)   ! rotation in spin space

  LOGICAL, INTENT(out) :: is_symmorphic, &  ! true if the gruop is symmorphic
       search_sym        ! true if gk

  CHARACTER(len=45), INTENT(in) :: sname(48)

  REAL(DP) :: sr(3,3,48)
  INTEGER :: isym

  is_symmorphic=.true.
  search_sym=.true.

  DO isym=1,nsym
     is_symmorphic=( is_symmorphic.and.(ftau(1,isym)==0).and.  &
          (ftau(2,isym)==0).and.  &
          (ftau(3,isym)==0) )
  ENDDO

  IF (.not.is_symmorphic) THEN
     DO isym=1,nsym
        search_sym=( search_sym.and.(gk(1,isym)==0).and.  &
             (gk(2,isym)==0).and.  &
             (gk(3,isym)==0) )
     ENDDO
  ENDIF
  !
  !  Set the group name, divide it in classes and set the irreducible
  !  representations
  !
  nsym_is=0
  DO isym=1,nsym
     CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
     IF (noncolin) THEN
        !
        !  In the noncollinear magnetic case finds the invariant subgroup of the point
        !  group of k. Presently we use only this subgroup to classify the levels.
        !
        IF (domag) THEN
           IF (t_rev(isym)==0) THEN
              nsym_is=nsym_is+1
              CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
              CALL find_u(sr_is(1,1,nsym_is),d_spin_is(1,1,nsym_is))
              s_is(:,:,nsym_is)=s(:,:,isym)
              gk_is(:,nsym_is)=gk(:,isym)
              ftau_is(:,nsym_is)=ftau(:,isym)
              sname_is(nsym_is)=sname(isym)
           ENDIF
        ELSE
           CALL find_u(sr(1,1,isym),d_spink(1,1,isym))
        ENDIF
     ENDIF
  ENDDO
  CALL find_group(nsym,sr,gname,code_group)
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
        CALL divide_class_so(code_group,nsym,sr,d_spink, &
             has_e,nclass,nelem_so,elem_so,which_irr_so)
     ENDIF
  ELSE
     CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
  ENDIF

  RETURN
END SUBROUTINE find_info_group
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE s_axis_to_cart (s, sr, at, bg)
  !----------------------------------------------------------------------
  !
  !     This routine transform a symmetry matrix expressed in the
  !     basis of the crystal axis in the cartesian basis.
  !
  !     last revised 3 may 1995 by A. Dal Corso
  !
  !
  USE kinds
  IMPLICIT NONE
  !
  !     first the input parameters
  !
  INTEGER :: s (3, 3)
  ! input: matrix in crystal axis
  real(DP) :: sr (3, 3), at (3, 3), bg (3, 3)
  ! output: matrix in cartesian axis
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  !
  !     here the local variable
  !

  INTEGER :: apol, bpol, kpol, lpol
  !
  !     counters on polarizations
  !
  DO apol = 1, 3
     DO bpol = 1, 3
        sr (apol, bpol) = 0.d0
        DO kpol = 1, 3
           DO lpol = 1, 3
              sr (apol, bpol) = sr (apol, bpol) + at (apol, kpol) * &
                   dble ( s (lpol, kpol) ) * bg (bpol, lpol)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE s_axis_to_cart
