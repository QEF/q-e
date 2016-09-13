!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#define CI   ( 0.D0, 1.D0 )
!
! ... written by Manu Sharma ( 2001-2005 )
!
!----------------------------------------------------------------------------
SUBROUTINE wf( clwf, c, bec, eigr, eigrb, taub, irb, &
               b1, b2, b3, Uall, what1, wfc, jw, ibrav )
  !----------------------------------------------------------------------------
  !
  ! ... this routine calculates overlap matrices
  !
  ! ... routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,                    ONLY : DP
  USE constants,                ONLY : pi, tpi
  USE ions_base,                ONLY : nsp, na, nax, nat
  USE uspp_param,               ONLY : nvb, ish
  USE cell_base,                ONLY : omega, at, alat, h, ainv
  USE electrons_base,           ONLY : nbspx, nbsp, nupdwn, iupdwn, nspin
  USE smallbox_gvec,                    ONLY : npb, nmb, ngb
  USE gvecw,                    ONLY : ngw
  USE gvect,       ONLY : gstart
  USE control_flags,            ONLY : iverbosity,conv_elec
  USE qgb_mod,                  ONLY : qgb
  USE wannier_base,             ONLY : wfg, nw, weight, indexplus, indexplusz, &
                                       indexminus, indexminusz, tag, tagp,     &
                                       expo, wfsd
  USE uspp_param,               ONLY : nh, nhm
  USE uspp,                     ONLY : nkb
  USE io_global,                ONLY : ionode, stdout
  USE mp,                       ONLY : mp_barrier, mp_sum
  USE mp_wave,                  ONLY : redistwf
  USE mp_global,                ONLY : nproc_bgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  USE fft_interfaces,           ONLY : invfft
  USE fft_base,                 ONLY : dfftp, dfftb
  USE printout_base,            ONLY : printout_base_open, printout_base_unit, &
                                       printout_base_close
  USE cp_main_variables,        ONLY : nfi, iprint_stdout
  USE time_step,                ONLY : tps
  USE input_parameters,         ONLY : tcpbo
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: irb(3,nat), jw, ibrav, clwf
  REAL(DP),    INTENT(INOUT) :: bec(nkb,nbsp)
  REAL(DP),    INTENT(IN)    :: b1(3), b2(3), b3(3), taub(3,nax)
  COMPLEX(DP), INTENT(INOUT) :: c(ngw,nbspx)
  COMPLEX(DP), INTENT(IN)    :: eigr(ngw,nat), eigrb(ngb,nat)
  REAL(DP),    INTENT(INOUT) :: Uall(nbsp,nbsp)
  LOGICAL,     INTENT(IN)    :: what1
  REAL(DP),    INTENT(OUT)   :: wfc(3,nbsp)
  !
  REAL(DP),    ALLOCATABLE :: becwf(:,:), temp3(:,:)
  COMPLEX(DP), ALLOCATABLE :: cwf(:,:), bec2(:), bec3(:), bec2up(:)
  COMPLEX(DP), ALLOCATABLE :: bec2dw(:), bec3up(:), bec3dw(:)
  COMPLEX(DP), ALLOCATABLE :: c_m(:,:), c_p(:,:), c_psp(:,:)
  COMPLEX(DP), ALLOCATABLE :: c_msp(:,:)
  INTEGER,     ALLOCATABLE :: tagz(:)
  REAL(DP),    ALLOCATABLE :: Uspin(:,:)
  COMPLEX(DP), ALLOCATABLE :: X(:,:), Xsp(:,:), X2(:,:), X3(:,:)
  COMPLEX(DP), ALLOCATABLE :: O(:,:,:), Ospin(:,:,:), Oa(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: qv(:)
  REAL(DP),    ALLOCATABLE :: gr(:,:), mt(:), mt0(:), wr(:), W(:,:), EW(:,:)
  INTEGER,     ALLOCATABLE :: f3(:), f4(:)
  COMPLEX(DP), ALLOCATABLE :: U2(:,:)
  !
  INTEGER           :: inl, jnl, iss, isa, is, ia, ijv, i, j, k, l, ig, &
                       ierr, ti, tj, tk, iv, jv, inw, iqv, ibig1, ibig2, &
                       ibig3, ir1, ir2, ir3, ir, m,  &
                       ib, jb, total, nstat, jj, ngpww, irb3
  REAL(DP)    :: t1, t2, t3, taup(3)
  REAL(DP)    :: wrsq, wrsqmin
  COMPLEX(DP) :: qvt
  REAL (DP)   :: temp_vec(3)
  INTEGER           :: adjust,ini, ierr1,nnn, me
  INTEGER           :: igx, igy, igz
  REAL(DP)    :: wfcx, wfcy, wfcz
  REAL(DP)    :: te(6)
  INTEGER     :: iunit
  
  COMPLEX(DP), EXTERNAL :: boxdotgridcplx
  !
#if defined (__MPI)
  !
  INTEGER :: proc, ntot, ncol, mc, ngpwpp(nproc_bgrp)
  INTEGER :: ncol1,nz1, nz_1 
  INTEGER :: nmin(3), nmax(3), n1,n2,nzx,nz,nz_
  INTEGER :: nmin1(3), nmax1(3)
  !
  COMPLEX(DP), ALLOCATABLE :: psitot(:,:), psitot_pl(:,:)
  COMPLEX(DP), ALLOCATABLE :: psitot_mi(:,:)
  INTEGER,     ALLOCATABLE :: ns(:)
  !
#endif
  !
  CALL start_clock('wf_1')
  !
  me = me_bgrp + 1
  !
  ALLOCATE( becwf(nkb,nbsp), temp3(nkb,nbsp), U2(nbsp,nbsp) )
  ALLOCATE( cwf(ngw,nbspx), bec2(nbsp), bec3(nbsp), bec2up(nupdwn(1)) )
  ALLOCATE( bec3up( nupdwn(1) ) )
  IF( nspin == 2 ) THEN
     ALLOCATE( bec2dw( nupdwn(2) ), bec3dw( nupdwn(2) ) )
  ENDIF
  ! 
  te = 0.D0
  !
  ALLOCATE( tagz( nw ))
  !
  tagz(:) = 1
  tagz(3) = 0
  !
  ! ... set up matrix O
  !
  ALLOCATE( O( nw, nbsp, nbsp ), X( nbsp, nbsp ), Oa( nw, nbsp, nbsp ) )
  !
  IF ( nspin == 2 .AND. nvb > 0 ) THEN
     !
     ALLOCATE( X2( nupdwn(1), nupdwn(1) ) )
     ALLOCATE( X3( nupdwn(2), nupdwn(2) ) )
     !
  END IF
  !
#if defined (__MPI)
  !
  ! Compute the number of states to each processor
  !
  ALLOCATE( ns( nproc_bgrp ) )
  ns = nbsp / nproc_bgrp
  DO j = 1, nbsp
     IF( (j-1) < MOD( nbsp, nproc_bgrp ) ) ns( j ) = ns( j ) + 1 
  END DO
  IF( iverbosity > 2) THEN
     DO j=1,nproc_bgrp
        WRITE( stdout, * ) ns(j)
     END DO
  END IF
  !
  nstat = ns( me )

  total = 0   
  DO proc=1,nproc_bgrp
     ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
     total=total+ngpwpp(proc)
     IF( iverbosity > 2) THEN
        WRITE( stdout, * ) "I am proceessor", proc, "and i have ",ns(me)," states."
     END IF
  END DO
  !
  ALLOCATE(psitot(total,nstat))
  ALLOCATE(psitot_pl(total,nstat))
  ALLOCATE(psitot_mi(total,nstat))

  ALLOCATE(c_p(ngw,nbspx))
  ALLOCATE(c_m(ngw,nbspx))
  IF( iverbosity > 2) THEN
     WRITE( stdout, * ) "All allocations done"
  END IF
  !
  ! ... Step 1. Communicate to all Procs so that each proc has all
  ! ... G-vectors and some states instead of all states and some
  ! ... G-vectors. This information is stored in the 1-d array 
  ! ... psitot1.
  !
  !   Step 2. Convert the 1-d array psitot1 into a 2-d array consistent with the
  !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
  !   and nstat states instead of all nbsp states
  !
  !
  CALL redistwf( c, psitot, ngpwpp, ns, intra_bgrp_comm, 1 )
  !
#endif   

  IF( clwf .EQ. 5 ) THEN
     !
     CALL write_psi( c, jw )
     !
  END IF
  !
  !
#if defined (__MPI)
  !
  !   Step 3. do the translation of the 2-d array to get the transtalted
  !   arrays psitot_pl and psittot_mi, corresponding to G+G' and -G+G'
  !   
  DO inw=1,nw   
   !
   !   Intermediate Check. If the translation is only along the z-direction
   !   no interprocessor communication and data rearrangement is required 
   !   because each processor contains all the G- components in the z-dir.
   !
   IF(tagz(inw).EQ.0) THEN
     DO i=1,nbsp
        DO ig=1,ngw
           IF(indexplusz(ig).EQ.-1) THEN
              c_p(ig,i)=(0.D0,0.D0)
           ELSE
              c_p(ig,i)=c(indexplusz(ig),i)
           END IF
           IF(indexminusz(ig).EQ.-1) THEN
              c_m(ig,i)=(0.D0,0.D0)
           ELSE
              c_m(ig,i)=CONJG(c(indexminusz(ig),i))
           END IF
        END DO
     END DO
   ELSE
     DO i=1,ns(me)
        DO ig=1,total
           IF(indexplus(ig,inw).EQ.-1) THEN
              psitot_pl(ig,i)=(0.D0,0.D0)
           ELSE   
              IF(tagp(ig,inw).EQ.1) THEN
                 psitot_pl(ig,i)=CONJG(psitot(indexplus(ig,inw),i))
              ELSE
                 psitot_pl(ig,i)=psitot(indexplus(ig,inw),i)
              END IF
           END IF
           IF(indexminus(ig,inw).EQ.-1) THEN
              psitot_mi(ig,i)=(0.D0,0.D0)
           ELSE
              IF(tag(ig,inw).EQ.1) THEN
                 psitot_mi(ig,i)=CONJG(psitot(indexminus(ig,inw),i))
              ELSE
                 psitot_mi(ig,i)=psitot(indexminus(ig,inw),i)
              END IF
           END IF
        END DO
     END DO
     IF( iverbosity > 2) THEN
        WRITE( stdout, * ) "Step 3. do the translation of the 2-d array...Done, wf"
     END IF
     !
     !   Step 4. Convert the 2-d arrays psitot_p and psitot_m into 1-d
     !   arrays
     !
     !   Step 5. Redistribute among processors. The result is stored in 2-d
     !   arrays c_p and c_m consistent with the notation c(ngw,nbsp), such that
     !   c_p(j,i) contains the coeffiCIent for c(j,i) corresponding to G+G'
     !       and c_m(j,i) contains the coeffiCIent for c(j,i) corresponding to -G+G'
     !
     c_p = 0.D0
     CALL redistwf( c_p, psitot_pl, ngpwpp, ns, intra_bgrp_comm, -1 )
     !
     c_m = 0.D0
     CALL redistwf( c_m, psitot_mi, ngpwpp, ns, intra_bgrp_comm, -1 )
     !
  END IF
    !
#else
    !
  ALLOCATE(c_p(ngw,nbspx))
  ALLOCATE(c_m(ngw,nbspx))
  DO inw=1,nw
     IF(tagz(inw).EQ.0) THEN
        DO i=1,nbsp
           DO ig=1,ngw
              IF(indexplusz(ig).EQ.-1) THEN
                 c_p(ig,i)=(0.D0,0.D0)
              ELSE
                 c_p(ig,i)=c(indexplusz(ig),i)
              END IF
              IF(indexminusz(ig).EQ.-1) THEN
                 c_m(ig,i)=(0.D0,0.D0)
              ELSE
                 c_m(ig,i)=CONJG(c(indexminusz(ig),i))
              END IF
           END DO
        END DO
     ELSE
        DO i=1,nbsp
           DO ig=1,ngw
              IF(indexplus(ig,inw).EQ.-1) THEN
                 c_p(ig,i)=(0.D0,0.D0)
              ELSE
                 IF(tagp(ig,inw).EQ.1) THEN
                    c_p(ig,i)=CONJG(c(indexplus(ig,inw),i))
                 ELSE
                    c_p(ig,i)=c(indexplus(ig,inw),i)
                 END IF
              END IF
              IF(indexminus(ig,inw).EQ.-1) THEN
                 c_m(ig,i)=(0.D0,0.D0)
              ELSE
                 IF(tag(ig,inw).EQ.1) THEN
                    c_m(ig,i)=CONJG(c(indexminus(ig,inw),i))
                 ELSE
                    c_m(ig,i)=c(indexminus(ig,inw),i)
                 END IF
              END IF
           END DO
        END DO
     END IF
    !
#endif
     !
     ! ... Step 6. Calculate Overlaps
     !
     ! ... Augmentation Part first
     !
     ALLOCATE( qv( dfftb%nnr ) )
     !
     X = ZERO
     !
     isa = 1
     DO is = 1, nvb
        DO ia =1, na(is)
           DO iv = 1, nh(is)
              inl = ish(is) + (iv-1)*na(is) + ia
              jv = iv 
              ijv=(jv-1)*jv/2 + iv
              qv( 1 : dfftb%nnr ) = 0.D0 
              DO ig=1,ngb
                 qv(npb(ig))=eigrb(ig,isa)*qgb(ig,ijv,is)
                 qv(nmb(ig))=CONJG(eigrb(ig,isa)*qgb(ig,ijv,is))
              END DO
#if defined(__MPI)
              irb3=irb(3,isa)
#endif
              CALL invfft(qv,dfftb,isa)
              iqv=1
              qvt=(0.D0,0.D0)
              qvt=boxdotgridcplx(irb(1,isa),qv,expo(1,inw))

#if defined(__MPI)
              CALL mp_sum( qvt, intra_bgrp_comm )
#endif
              !
              IF (nspin.EQ.1) THEN
                 bec2(1:nbsp)=(0.D0,0.D0)
                 bec2(1:nbsp)=bec(inl,1:nbsp)*ONE
                 CALL ZSYRK('U','T',nbsp,1,qvt,bec2,1,ONE,X,nbsp)
              ELSE
                 X2=(0.D0,0.D0)
                 X3=(0.D0,0.D0)
                 bec2up(1:nupdwn(1))=(0.D0,0.D0)
                 bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))
                 CALL ZSYRK('U','T',nupdwn(1),1,qvt,bec2up,1,ONE,X2,nupdwn(1))
                 bec2dw(1:nupdwn(2))=(0.D0,0.D0)
                 bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):nbsp)
                 CALL ZSYRK('U','T',nupdwn(2),1,qvt,bec2dw,1,ONE,X3,nupdwn(2))
                 DO i = 1, nupdwn(1)
                    DO j=i, nupdwn(1)
                       X(i,j)=X(i,j)+X2(i,j)
                    END DO
                 END DO
                 DO i = 1,nupdwn(2)
                    DO j=i,nupdwn(2)
                       X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                    END DO
                 END DO
              END IF
              DO jv = iv+1, nh(is)
                 jnl = ish(is) + (jv-1)*na(is) + ia
                 ijv = (jv-1)*jv/2 + iv
                 qv( 1:dfftb%nnr ) = 0.D0
                 DO ig=1,ngb
                    qv(npb(ig))=eigrb(ig,isa)*qgb(ig,ijv,is)
                    qv(nmb(ig))=CONJG(eigrb(ig,isa)*qgb(ig,ijv,is))
                 END DO
                 CALL invfft(qv,dfftb,isa)
                 iqv=1
                 qvt=0.D0
                 qvt=boxdotgridcplx(irb(1,isa),qv,expo(1,inw))
#if defined(__MPI)
                 CALL mp_sum( qvt, intra_bgrp_comm )
#endif
                 !
                 IF (nspin.EQ.1) THEN
                    bec2(1:nbsp)=(0.D0,0.D0)
                    bec3(1:nbsp)=(0.D0,0.D0)
                    bec2(1:nbsp)=bec(inl,1:nbsp)*ONE
                    bec3(1:nbsp)=bec(jnl,1:nbsp)*ONE
                    CALL ZSYR2K('U','T',nbsp,1,qvt,bec2,1,bec3,1,ONE,X,nbsp)
                 ELSE
                    X2=(0.D0,0.D0)
                    X3=(0.D0,0.D0)
                    bec2up(1:nupdwn(1))=(0.D0,0.D0)
                    bec3up(1:nupdwn(1))=(0.D0,0.D0)
                    bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))*ONE
                    bec3up(1:nupdwn(1))=bec(jnl,1:nupdwn(1))*ONE
                    CALL ZSYR2K('U','T',nupdwn(1),1,qvt,bec2up,1,bec3up,1,ONE,X2,nupdwn(1))
                    bec2dw(1:nupdwn(2))=(0.D0,0.D0)
                    bec3dw(1:nupdwn(2))=(0.D0,0.D0)
                    bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):nbsp)*ONE
                    bec3dw(1:nupdwn(2))=bec(jnl,iupdwn(2):nbsp)*ONE
                    CALL ZSYR2K('U','T',nupdwn(2),1,qvt,bec2dw,1,bec3dw,1,ONE,X3,nupdwn(2))
                    DO i = 1, nupdwn(1)
                       DO j=i, nupdwn(1)
                          X(i,j)=X(i,j)+X2(i,j)
                       END DO
                    END DO
                    DO i = 1,nupdwn(2)
                       DO j=i,nupdwn(2)
                          X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                       END DO
                    END DO
                 END IF
              END DO
           END DO
           isa = isa + 1
        END DO
     END DO
     t1=omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     X=X*t1
     DO i=1, nbsp
        DO j=i+1, nbsp
           X(j, i)=X(i, j)
        END DO
     END DO
     Oa(inw, :, :)=X(:, :)
     IF( iverbosity > 2) THEN
        WRITE( stdout, * ) "Augmentation Part Done"
     END IF

     DEALLOCATE( qv )


     !   Then Soft Part
     IF( nspin == 1 ) THEN
        !   Spin Unpolarized calculation
        X=0.D0   
        IF( gstart == 2 ) THEN
           c_m(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:)
        CALL zgemm('C','N',nbsp,nbsp,ngw,ONE,c,ngw,c_p,ngw,ONE,X,nbsp)
        CALL zgemm('T','N',nbsp,nbsp,ngw,ONE,c,ngw,c_m,ngw,ONE,X,nbsp)

        CALL mp_sum ( X, intra_bgrp_comm )

        O(inw,:,:)=Oa(inw,:,:)+X(:,:)

        IF( iverbosity > 2) THEN
           WRITE( stdout, * ) "Soft Part Done"
        END IF

     ELSE
        !   Spin Polarized case
        !   Up Spin First
        ALLOCATE(Xsp(nbsp,nupdwn(1)))
        ALLOCATE(c_psp(ngw,nupdwn(1)))
        ALLOCATE(c_msp(ngw,nupdwn(1)))
        Xsp=0.D0
        c_psp=0.D0 
        c_msp=0.D0
        DO i=1,nupdwn(1)
           c_psp(:,i)=c_p(:,i)
           c_msp(:,i)=c_m(:,i)
        END DO
        IF(gstart.EQ.2) THEN
           c_msp(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:,1,1)
        CALL zgemm('C','N',nbsp,nupdwn(1),ngw,ONE,c,ngw,c_psp,ngw,ONE,Xsp,nbsp)
        CALL zgemm('T','N',nbsp,nupdwn(1),ngw,ONE,c,ngw,c_msp,ngw,ONE,Xsp,nbsp)
#if defined(__MPI)
        CALL mp_sum ( Xsp, intra_bgrp_comm )
#endif
        DO i=1,nupdwn(1)
           DO j=1,nbsp
              X(j,i)=Xsp(j,i)
           END DO
        END DO
        DEALLOCATE(Xsp,c_psp,c_msp)
        !    Then Down Spin
        ALLOCATE(Xsp(nbsp,iupdwn(2):nbsp))
        ALLOCATE(c_psp(ngw,iupdwn(2):nbsp))
        ALLOCATE(c_msp(ngw,iupdwn(2):nbsp))
        Xsp=0.D0
        c_psp=0.D0
        c_msp=0.D0
        DO i=iupdwn(2),nbsp
           c_psp(:,i)=c_p(:,i)
           c_msp(:,i)=c_m(:,i)
        END DO
        IF(gstart.EQ.2) THEN
           c_msp(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:,1,1)
        CALL zgemm('C','N',nbsp,nupdwn(2),ngw,ONE,c,ngw,c_psp,ngw,ONE,Xsp,nbsp)
        CALL zgemm('T','N',nbsp,nupdwn(2),ngw,ONE,c,ngw,c_msp,ngw,ONE,Xsp,nbsp)
#if defined(__MPI)
        CALL mp_sum ( Xsp, intra_bgrp_comm )
#endif
        DO i=iupdwn(2),nbsp
           DO j=1,nbsp
              X(j,i)=Xsp(j,i)
           END DO
        END DO
        DEALLOCATE(Xsp,c_psp,c_msp)
        O(inw,:,:)=Oa(inw,:,:)+X(:,:)
     END IF


  END DO

#if defined(__MPI)
  DEALLOCATE(ns)
#endif

  CALL stop_clock('wf_1')

  DEALLOCATE( X )
  IF ( ALLOCATED( X2 ) )  DEALLOCATE( X2 )
  IF ( ALLOCATED( X3 ) )  DEALLOCATE( X3 )
  !

  CALL start_clock('wf_2')


  IF(clwf.EQ.2) THEN
     !    output the overlap matrix to fort.38
     IF(me.EQ.1) THEN
        REWIND 38
        WRITE(38, '(i5, 2i2, i3, f9.5)') nbsp, nw, nspin, ibrav, alat
        IF (nspin.EQ.2) THEN
           WRITE(38, '(i5)') nupdwn(1)
        END IF
        WRITE(38, *) at(:,1)*alat
        WRITE(38, *) at(:,2)*alat
        WRITE(38, *) at(:,3)*alat
        WRITE(38, *) b1
        WRITE(38, *) b2
        WRITE(38, *) b3
        DO inw=1, nw
           WRITE(38, *) wfg(inw, :), weight(inw)
        END DO
        DO inw=1, nw
           DO i=1, nbsp
              DO j=1, nbsp
                 WRITE(38, *) O(inw, i, j)
              END DO
           END DO
        END DO
        DO i=1, nbsp
           DO j=1, nbsp
              WRITE(38, *) Uall(i, j)
           END DO
        END DO
        CLOSE(38)
     END IF
     CALL stop_run( .TRUE. )
  END IF

  IF(clwf.EQ.3.OR.clwf.EQ.4) THEN
     IF(nspin.EQ.1) THEN
        IF(.NOT.what1) THEN
           IF(wfsd==1) THEN
              CALL ddyn_u(nbsp,O,Uall, 1)
           ELSE IF(wfsd==2) THEN
              CALL wfsteep(nbsp,O,Uall,b1,b2,b3)
           ELSE IF(wfsd==3) THEN
              CALL jacobi_rotation(nbsp,O,Uall,b1,b2,b3)
           END IF
        END IF
        IF( iverbosity > 2) THEN
           WRITE( stdout, * ) "Out from DDYN"
        END IF
     ELSE
        ALLOCATE(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
        DO i=1, nupdwn(1)
           DO j=1, nupdwn(1)
              Uspin(i, j)=Uall(i, j)
              Ospin(:, i, j)=O(:, i, j)
           END DO
        END DO
        IF(.NOT.what1) THEN
           IF(wfsd==1) THEN
             CALL ddyn_u(nupdwn(1), Ospin, Uspin, 1)
           ELSE IF (wfsd==2) THEN
             CALL wfsteep(nupdwn(1), Ospin, Uspin,b1,b2,b3)
           ELSE
            CALL jacobi_rotation(nupdwn(1), Ospin, Uspin,b1,b2,b3)
           END IF
        END IF
        DO i=1, nupdwn(1)
           DO j=1, nupdwn(1)
              Uall(i, j)=Uspin(i, j)
              O(:,i,j)  =Ospin(:,i,j)
           END DO
        END DO
        DEALLOCATE(Uspin, Ospin)
        ALLOCATE(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
        DO i=1, nupdwn(2)
           DO j=1, nupdwn(2)
              Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
              Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
           END DO
        END DO
        IF(.NOT.what1) THEN
           IF(wfsd==1) THEN
              CALL ddyn_u(nupdwn(2), Ospin, Uspin, 2)
           ELSE IF (wfsd==2) THEN
              CALL wfsteep(nupdwn(2), Ospin, Uspin,b1,b2,b3)
           ELSE
              CALL jacobi_rotation(nupdwn(2), Ospin, Uspin,b1,b2,b3)
           END IF
        END IF
        DO i=1, nupdwn(2)
           DO j=1, nupdwn(2)
              Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
              O(:,i+nupdwn(1),j+nupdwn(1))=Ospin(:,i,j)
           END DO
        END DO
        DEALLOCATE(Uspin, Ospin)
     END IF
  END IF

  !       Update C and bec
  cwf=ZERO
  !        cwf(:,:)=c(:,:,1,1)
  becwf=0.0d0
  U2=Uall*ONE
  CALL zgemm('N','N',ngw,nbsp,nbsp,ONE,c,ngw,U2,nbsp,ZERO,cwf,ngw)
  !           call zgemm('nbsp','nbsp',ngw,nbsp,nbsp,ONE,cwf,ngw,U2,nbsp,ZERO,cwf,ngw)
  CALL dgemm('N','N',nkb,nbsp,nbsp,ONE,bec,nkb,Uall,nbsp,ZERO,becwf,nkb)
  U2=ZERO
  IF( iverbosity > 2) THEN
     WRITE( stdout, * ) "Updating Wafefunctions and Bec"
  END IF

  c(:,:)=cwf(:,:)
  bec(:,:)=becwf(:,:)

  IF( iverbosity > 2) THEN
     WRITE( stdout, * ) "Wafefunctions and Bec Updated"
  END IF
  !
  ! calculate wannier-function centers
  !
  ALLOCATE( wr(nw), W(nw,nw), gr(nw,3), EW(nw,nw), f3(nw), f4(nw), mt0(nw), mt(nw) )
  !
  DO inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  END DO
  !
  ! set up a matrix with the element (i,j) is G_i G_j weight(j)
  ! to check the correctness of choices on G vectors
  !
  DO i=1, nw
     DO j=1, nw
        W(i,j)=DOT_PRODUCT(gr(i,:),gr(j,:))*weight(j)
     END DO
  END DO
  !
  EW = W
  DO i=1,nw
     EW(i,i) = EW(i,i)-1.D0
  END DO
  !
  ! ... balance the phase factor if necessary
  !
  ! adjust mt : very inefficient routine added by Young-Su -> must be improved
  DO i=1, nbsp
     mt0(:) = -AIMAG(LOG(O(:,i,i)))/tpi
     wr = MATMUL(EW,mt0)
     wrsq = SUM(wr(:)**2)
     IF ( wrsq .lt. 1.D-6 ) THEN
        mt = mt0
     ELSE
        wrsqmin = 100.D0
COMB:   DO k=3**nw-1,0,-1
           tk=k
           DO j=nw,1,-1
              f3(j)=tk/3**(j-1)
              tk=tk-f3(j)*3**(j-1)
           END DO
           mt(:)=mt0(:)+f3(:)-1
           wr = MATMUL(EW,mt)
           wrsq = SUM(wr(:)**2)
           IF ( wrsq .lt. wrsqmin ) THEN
              wrsqmin = wrsq
              f4(:)=f3(:)-1
           END IF
        END DO COMB
        mt = mt0 + f4
     END IF
     !
     wfc(1, i) = SUM(mt*weight(:)*gr(:,1))*alat
     wfc(2, i) = SUM(mt*weight(:)*gr(:,2))*alat
     wfc(3, i) = SUM(mt*weight(:)*gr(:,3))*alat
     !
  END DO
  !
  DEALLOCATE( wr, W, gr, EW, f3, f4, mt0, mt )
  !
#if defined (__MPI)
  !
  DEALLOCATE( psitot )
  DEALLOCATE( psitot_pl )
  DEALLOCATE( psitot_mi )
  !
#endif
  !
  DEALLOCATE( c_p, c_m )
  !
  DEALLOCATE( O )
  DEALLOCATE( Oa )
  DEALLOCATE( tagz )
  DEALLOCATE( becwf, temp3, U2 )
  DEALLOCATE( cwf, bec2, bec3, bec2up, bec3up )
  IF( ALLOCATED( bec2dw ) ) DEALLOCATE( bec2dw )
  IF( ALLOCATED( bec3dw ) ) DEALLOCATE( bec3dw )

  CALL stop_clock('wf_2')
  !
  RETURN
  !
END SUBROUTINE wf
!
!----------------------------------------------------------------------------
SUBROUTINE ddyn( m, Omat, Umat, b1, b2, b3 )
  !----------------------------------------------------------------------------
  ! ... This part of the subroutine wf has been added by Manu. It performes
  ! ... Damped Dynamics on the A matrix to get the Unitary transformation to
  ! ... obtain the wannier function at time(t+delta). It also updates the
  ! ... quantities bec
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : ionode, stdout
  USE wannier_base,     ONLY : wf_friction, nsteps, tolw, adapt, wf_q, &
                               weight, nw, wfdt
  USE cell_base,        ONLY : alat
  USE constants,        ONLY : tpi, bohr_radius_angs
  USE electrons_base,   ONLY : nbsp
  USE control_flags,    ONLY : iverbosity,conv_elec
  USE mp_global,        ONLY : me_bgrp
  USE printout_base,    ONLY : printout_base_open, printout_base_unit, &
                               printout_base_close
  USE cp_main_variables,        ONLY : nfi, iprint_stdout
  USE time_step,                ONLY : tps
  USE input_parameters,         ONLY : tcpbo
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: f3(nw), f4(nw), i,j,inw
  INTEGER ,INTENT(in) :: m
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  REAL(DP), INTENT(inout) :: Umat(m,m)
  COMPLEX(DP), INTENT(inout) :: Omat(nw,m,m)
  COMPLEX(DP) :: U2(m,m),U3(m,m)
  INTEGER :: adjust,ini, ierr1,nnn
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: W
  REAL(DP) :: t0, fric,U(m,m), t2
  REAL(DP) :: A(m,m),oldt0,Wm(m,m),U1(m,m)
  REAL(DP) :: Aminus(m,m), Aplus(m,m),f2(4*m)
  REAL(DP) :: temp(m,m)
  COMPLEX(DP) :: d(m,m)
  COMPLEX(DP) :: f1(2*m-1), wp(m*(m+1)/2),z(m,m)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:, :) :: X1
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:, :, :) :: Oc
  REAL(DP) , ALLOCATABLE , DIMENSION(:) :: mt
  REAL(DP) :: spread, sp
  REAL(DP) :: wfc(3,nbsp), gr(nw,3)
  INTEGER  :: me, iunit
  !
  me = me_bgrp + 1
  !
  ALLOCATE(mt(nw))
  ALLOCATE(X1(m,m))
  ALLOCATE(Oc(nw,m,m))

  fric=wf_friction
  ALLOCATE (W(m,m),wr(m))

  Umat=0.D0
  DO i=1,m
     Umat(i,i)=1.D0
  END DO

  U2=Umat*ONE

  !
  ! update Oc using the initial guess of Uspin
  !
  DO inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=ZERO
     CALL zgemm ('T', 'N', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
     X1=ZERO
     CALL zgemm ('N','N', m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
     Oc(inw, :, :)=X1(:, :)
  END DO

  U2=ZERO
  U3=ZERO

  oldt0=0.D0
  A=0.D0
  Aminus=A
  temp=Aminus


  !   START ITERATIONS HERE

  DO ini=1, nsteps

     t0=0.D0     !use t0 to store the value of omega
     DO inw=1, nw
        DO i=1, m
           t0=t0+DBLE(CONJG(Oc(inw, i, i))*Oc(inw, i, i))
        END DO
     END DO

     IF(ABS(t0-oldt0).LT.tolw) THEN
        IF(me.EQ.1) THEN
           WRITE(stdout,*) "MLWF Generated at Step",ini,"  Convergence = ",ABS(t0-oldt0)
        END IF
        IF( iverbosity > 2) THEN
           WRITE(stdout,*) "MLWF Generated at Step",ini,"  Convergence = ",ABS(t0-oldt0)
        END IF
        GO TO 241
     END IF

     IF(adapt) THEN
        IF(oldt0.LT.t0) THEN
           fric=fric/2.d0
           A=Aminus
           Aminus=temp
        END IF
     END IF

     !   calculate d(omega)/dA and store result in W
     !   this is the force for the damped dynamics
     !

     W=0.D0
     DO inw=1, nw
        t2=weight(inw)
        DO i=1,m
           DO j=1,m
              W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*CONJG(Oc(inw,i,i)        &
                   -Oc(inw,j,j))+CONJG(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
           END DO
        END DO
     END DO


     !   the verlet scheme to calculate A(t+wfdt)

     Aplus=0.D0

     DO i=1,m
        DO j=i+1,m
           Aplus(i,j)=Aplus(i,j)+(2*wfdt/(2*wfdt+fric))*(2*A(i,j)               &
                -Aminus(i,j)+(wfdt*wfdt/wf_q)*W(i,j)) + (fric/(2*wfdt+fric))*Aminus(i,j)
        ENDDO
     ENDDO

     Aplus=Aplus-TRANSPOSE(Aplus)
     Aplus=(Aplus-A)

     DO i=1, m
        DO j=i,m 
           wp(i + (j-1)*j/2) = CMPLX(0.d0, Aplus(i,j),kind=DP)
        END DO
     END DO

     CALL zhpev('V','U',m,wp,wr,z,m,f1,f2,ierr1)

     IF (ierr1.NE.0) THEN 
        WRITE( stdout, * ) "failed to diagonalize W!"
        STOP
     END IF

     d=0.D0
     DO i=1, m
        d(i, i)=EXP(CI*wr(i)*wfdt)
     END DO      !d=exp(d)

     !   U=z*exp(d)*z+
     !   
     U3=ZERO
     CALL zgemm ('N', 'N', m,m,m,ONE,z,m,d,m,ZERO,U3,m)  
     U2=ZERO
     CALL zgemm ('N','C', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
     U=DBLE(U2)
     U2=ZERO
     U3=ZERO

     temp=Aminus
     Aminus=A
     A=Aplus


     !   update Umat
     !
     U1=ZERO
     CALL dgemm ('N', 'N', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)

     Umat=U1 

     !   update Oc
     !
     U2=Umat*ONE
     U3=ZERO
     DO inw=1, nw
        X1(:, :)=Omat(inw, :, :)
        CALL zgemm ('T', 'N', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
        X1=ZERO
        CALL zgemm ('N','N',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
        Oc(inw, :, :)=X1(:, :)
     END DO
     U2=ZERO
     U3=ZERO

     IF(ABS(t0-oldt0).GE.tolw.AND.ini.GE.nsteps) THEN
        IF(me.EQ.1) THEN
           WRITE(stdout,*) "MLWF Generated at Step",ini,"  Convergence = ",ABS(t0-oldt0)
        END IF
        IF( iverbosity > 2) THEN
           WRITE(stdout,*) "MLWF Generated at Step",ini,"  Convergence = ",ABS(t0-oldt0)
        END IF
        GO TO 241
     END IF

     oldt0=t0

  END DO

241 DEALLOCATE(wr, W)

  !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
  IF( ( MOD( nfi, iprint_stdout ) == 0 ) )  THEN
     IF ( ionode ) THEN
        iunit = printout_base_unit( "spr" )
        CALL printout_base_open( "spr" )
        IF (.NOT.tcpbo) THEN
           WRITE( iunit, '(I8,f16.8)' )nfi,tps
        ELSE   
           IF(conv_elec) WRITE( iunit, '(I8,f16.8)' )nfi,tps
        END IF    
     END IF    
  END IF    

  spread=0.0d0

  DO i=1, m
     !
     mt=1.D0-DBLE(Oc(:,i,i)*CONJG(Oc(:,i,i)))
     sp = (alat*bohr_radius_angs/tpi)**2*SUM(mt*weight)
     !
     !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
     IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
        IF ( ionode ) THEN
           IF (.NOT.tcpbo) THEN
              WRITE( iunit, '(f20.14)' ) sp
           ELSE
             IF(conv_elec) WRITE( iunit, '(f20.14)' ) sp
           END IF   
        END IF   
     END IF
     !
     IF ( sp < 0.D0 ) &
        CALL errore( 'cp-wf', 'Something wrong WF Spread negative', 1 )
     !
     spread=spread+sp
     !
  END DO

  !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
  IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
     IF ( ionode ) CALL printout_base_close( "spr" )
  END IF

  spread=spread/m

  IF(me.EQ.1) THEN
     WRITE(24, '(f10.7)') spread
     WRITE(27,*) "Average spread = ", spread
  END IF
  Omat=Oc
  IF( iverbosity > 2) THEN
     WRITE( stdout, * ) "Average spread = ", spread
  END IF
  !
  DEALLOCATE (mt,X1,Oc)
  !
  IF( iverbosity > 2) THEN
     WRITE( stdout, * ) "Leaving DDYN"
  END IF
  RETURN
END SUBROUTINE ddyn
!
!----------------------------------------------------------------------------
SUBROUTINE wfunc_init( clwf, b1, b2, b3, ibrav )
  !----------------------------------------------------------------------------
  !
  USE io_global,          ONLY : stdout
  USE kinds,              ONLY : DP
  USE gvect, ONLY : g, mill, gstart
  USE gvecw,              ONLY : ngw
  USE electrons_base,     ONLY : nbsp
  USE wannier_base,       ONLY : gnx, gnn, indexplus, indexminus, &
                                 indexplusz, indexminusz, tag, tagp, &
                                 wfg, weight, nw
  USE uspp_param,         ONLY : nvb
  USE mp,                 ONLY : mp_barrier, mp_bcast, mp_gather, mp_set_displs, mp_sum
  USE mp_global,          ONLY : nproc_bgrp, me_bgrp, intra_bgrp_comm, root_bgrp
  USE fft_base,           ONLY : dfftp
  USE parallel_include     
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  INTEGER,  INTENT(in) :: clwf, ibrav
#if defined(__MPI)
  INTEGER :: ntot, proc, ierr, i,j,inw,ngppp(nproc_bgrp)
  INTEGER :: ii,ig,displs(nproc_bgrp)
#else
  INTEGER :: ierr, i,j,inw, ntot
  INTEGER :: ii,ig
#endif
  REAL (DP), ALLOCATABLE:: bigg(:,:)
  INTEGER, ALLOCATABLE :: bign(:,:)
  INTEGER :: igcount,nw1,jj,nw2, in, kk
  INTEGER, ALLOCATABLE :: i_1(:), j_1(:), k_1(:)
  INTEGER :: ti, tj, tk
  REAL(DP) ::t1, vt
  INTEGER :: me
  INTEGER :: ng_me          !HK: the number of g vector in proc me
  INTEGER :: g_start, g_end !HK: ig loop range for proc me
  !
  CALL start_clock('wf_init')

  me = me_bgrp + 1
  !
! IF ( nbsp < nproc_bgrp ) &
!    CALL errore( 'cp-wf', &
!               & 'Number of Processors is greater than the number of states', 1 )
  !
  ALLOCATE(gnx(3,ngw))
  ALLOCATE(gnn(3,ngw))
  vt=1.0d-4
  j=0
  DO i=1,ngw
     gnx(1,i)=g(1,i)
     gnx(2,i)=g(2,i)
     gnx(3,i)=g(3,i)
     gnn(1,i)=mill(1,i)
     gnn(2,i)=mill(2,i)
     gnn(3,i)=mill(3,i)
  END DO

#if defined(__MPI)

  ntot=0
  DO i=1,nproc_bgrp
     ngppp(i)=(dfftp%nwl(i)+1)/2
  END DO

  CALL mp_set_displs( ngppp, displs, ntot, nproc_bgrp )

  ALLOCATE(bigg(3,ntot)); bigg = 0._DP
  ALLOCATE(bign(3,ntot)); bign= 0

#else
  ntot=ngw
  ALLOCATE(bigg(3,ntot))
  ALLOCATE(bign(3,ntot))
  bigg(1:3,1:ntot)=gnx(1:3,1:ntot)
  bign(1:3,1:ntot)=gnn(1:3,1:ntot)
#endif
  !
  CALL setwfg( ibrav, b1, b2, b3 )
  !
  nw1 = nw

  WRITE( stdout, * ) "WANNIER SETUP : check G vectors and weights"
  DO i=1,nw1 
     WRITE( stdout,'("inw = ",I1,":",3I4,F11.6)') i,wfg(i,:), weight(i)
  END DO
  
  WRITE( stdout, * ) "Translations to be done", nw1
  ALLOCATE(indexplus(ntot,nw1))
  ALLOCATE(indexminus(ntot,nw1))
  ALLOCATE(tag(ntot,nw1))
  ALLOCATE(tagp(ntot,nw1))
  ALLOCATE(indexplusz(ngw))
  ALLOCATE(indexminusz(ngw))
  ALLOCATE(i_1(nw1))
  ALLOCATE(j_1(nw1))
  ALLOCATE(k_1(nw1))

  indexplus=0
  indexminus=0
  tag=0
  tagp=0
  indexplusz=0
  indexminusz=0
  i_1(:)=wfg(:,1)
  j_1(:)=wfg(:,2)
  k_1(:)=wfg(:,3)


  WRITE( stdout, * ) "ibrav selected:", ibrav
  !
  IF(nvb.GT.0) CALL small_box_wf(i_1, j_1, k_1, nw1)

#if defined(__MPI)
  !
  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_gather( gnx, bigg, ngppp, displs, root_bgrp, intra_bgrp_comm )
  !
  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_gather( gnn, bign, ngppp, displs, root_bgrp, intra_bgrp_comm )
  !
  CALL mp_bcast( bigg,  root_bgrp, intra_bgrp_comm ) !parallelize wf_init with MPI (serial could be very slow for large system or high cutoff)
  CALL mp_bcast( bign,  root_bgrp, intra_bgrp_comm ) !parallelization could be done better with MPI_SCATHER and MPI_GATHER (todo)
  ng_me = ceiling(dble(ntot)/dble(nproc_bgrp))
  g_start = (me-1)*ng_me + 1
  g_end   = me*ng_me
  if (g_end.gt.ntot)   g_end = ntot !HK: truncate upper bound to ntot
  !
#else
  g_start = 1
  g_end   = ntot
#endif

  IF(me.EQ.1) THEN
     IF(clwf.EQ.5) THEN
#if defined(__MPI)
        DO ii=1,ntot
           WRITE(21,*) bigg(:,ii)
        END DO
#else
        DO ii=1,ngw
           WRITE(21,*) g(1,ii), g(2,ii), g(3,ii)
        END DO
#endif
        CLOSE(21)
     END IF
  END IF

  DO inw=1,nw1
     IF(i_1(inw).EQ.0.AND.j_1(inw).EQ.0) THEN
        IF(gstart.EQ.2) THEN
           indexminusz(1)=-1
        END IF
!$omp parallel do private(ti,tj,tk)
        DO ig=1,ngw
           ti=(gnn(1,ig)+i_1(inw))
           tj=(gnn(2,ig)+j_1(inw))
           tk=(gnn(3,ig)+k_1(inw))
           DO ii=1,ngw
              IF(gnn(1,ii).EQ.ti.AND.gnn(2,ii).EQ.tj.AND.gnn(3,ii).EQ.tk) THEN
                 indexplusz(ig)=ii
                 GO TO 224
              END IF
           END DO
           indexplusz(ig)=-1
224        ti=(-gnn(1,ig)+i_1(inw))
           tj=(-gnn(2,ig)+j_1(inw))
           tk=(-gnn(3,ig)+k_1(inw))
           IF(ti.LT.0.OR.(ti.EQ.0.AND.(tj.LT.0.OR.(tj.EQ.0.AND.tk.LT.0)))) THEN
              DO ii=1,ngw
                 IF(gnn(1,ii).EQ.-ti.AND.gnn(2,ii).EQ.-tj.AND.gnn(3,ii).EQ.-tk) THEN
                    indexminusz(ig)=ii
                    GO TO 223
                 END IF
              END DO
              indexminusz(ig)=-1
           ELSE
              DO ii=1,ngw
                 IF(gnn(1,ii).EQ.ti.AND.gnn(2,ii).EQ.tj.AND.gnn(3,ii).EQ.tk) THEN
                    indexminusz(ig)=ii
                    GO TO 223
                 END IF
              END DO
              indexminusz(ig)=-1
           END IF
223        CONTINUE
        END DO
!$omp end parallel do
        WRITE( stdout, * ) "Translation", inw, "for", ngw, "G vectors"
     ELSE
#if defined(__MPI)
        IF(me.EQ.1) THEN   
#endif
           IF(gstart.EQ.2) THEN
              indexminus(1,inw)=-1
           END IF
#if defined(__MPI)
        END IF
#endif
!$omp parallel do private(ti,tj,tk)
        DO ig=g_start,g_end
           ti=(bign(1,ig)+i_1(inw))
           tj=(bign(2,ig)+j_1(inw))
           tk=(bign(3,ig)+k_1(inw))
           IF(ti.LT.0.OR.(ti.EQ.0.AND.(tj.LT.0.OR.(tj.EQ.0.AND.tk.LT.0)))) THEN
              DO ii=1,ntot
                 IF(bign(1,ii).EQ.-ti.AND.bign(2,ii).EQ.-tj.AND.bign(3,ii).EQ.-tk) THEN
                    indexplus(ig,inw)=ii
                    tagp(ig,inw)=1
                    GO TO 214
                 END IF
              END DO
              indexplus(ig,inw)=-1
              tagp(ig,inw)=1
           ELSE
              DO ii=1,ntot
                 IF(bign(1,ii).EQ.ti.AND.bign(2,ii).EQ.tj.AND.bign(3,ii).EQ.tk) THEN
                    indexplus(ig,inw)=ii
                    tagp(ig,inw)=-1
                    GO TO 214
                 END IF
              END DO
              indexplus(ig,inw)=-1
              tagp(ig,inw)=-1
           END IF
214        ti=-bign(1,ig)+i_1(inw)
           tj=-bign(2,ig)+j_1(inw)
           tk=-bign(3,ig)+k_1(inw)
           IF(ti.LT.0.OR.(ti.EQ.0.AND.(tj.LT.0.OR.(tj.EQ.0.AND.tk.LT.0)))) THEN
              DO ii=1,ntot
                 IF(bign(1,ii).EQ.-ti.AND.bign(2,ii).EQ.-tj.AND.bign(3,ii).EQ.-tk) THEN
                    indexminus(ig,inw)=ii
                    tag(ig,inw)=1
                    GO TO 213
                 END IF
              END DO
              indexminus(ig,inw)=-1
              tag(ig,inw)=1
           ELSE 
              DO ii=1,ntot
                 IF(bign(1,ii).EQ.ti.AND.bign(2,ii).EQ.tj.AND.bign(3,ii).EQ.tk) THEN
                    indexminus(ig,inw)=ii
                    tag(ig,inw)=-1
                    GO TO 213
                 END IF
              END DO
              indexminus(ig,inw)=-1
              tag(ig,inw)=-1
           END IF
213        CONTINUE
        END DO
!$omp end parallel do 
        WRITE( stdout, * ) "Translation", inw, "for", ntot, "G vectors"
     END IF
  END DO

#if defined(__MPI)

  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_sum( indexplus,  intra_bgrp_comm )
  CALL mp_sum( indexminus, intra_bgrp_comm )
  CALL mp_sum( tag,        intra_bgrp_comm )
  CALL mp_sum( tagp,       intra_bgrp_comm )
  !
  CALL mp_bcast( indexplus,  root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( indexminus, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( tag,        root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( tagp,       root_bgrp, intra_bgrp_comm )

#endif
  DEALLOCATE(bigg)
  DEALLOCATE(bign)
  DEALLOCATE(i_1,j_1,k_1)

  CALL stop_clock('wf_init')

  RETURN
END SUBROUTINE wfunc_init
!
!----------------------------------------------------------------------------
SUBROUTINE grid_map()
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE efcalc,                 ONLY : xdist, ydist, zdist
  USE fft_base,               ONLY : dffts
  USE mp_global,              ONLY : me_bgrp
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: ir1, ir2, ir3, ibig3, me, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx
  !
  me = me_bgrp + 1
  !
  ALLOCATE(xdist(dffts%nnr))
  ALLOCATE(ydist(dffts%nnr))
  ALLOCATE(zdist(dffts%nnr))
  !
  nr1s = dffts%nr1
  nr2s = dffts%nr2
  nr1sx = dffts%nr1x
  nr2sx = dffts%nr2x
  nr3sx = dffts%nr3x
  nr3s = dffts%nr3
  DO ir3=1,nr3s
#if defined(__MPI)
     ibig3 = ir3 - dffts%ipp( me )
     IF(ibig3.GT.0.AND.ibig3.LE.dffts%npp(me)) THEN
#else
        ibig3=ir3
#endif
        DO ir2=1,nr2s
           DO ir1=1,nr1s
              xdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
                   &                  ((ir1-1)/DBLE(nr1sx))
              ydist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                   &
                   &                  ((ir2-1)/DBLE(nr2sx))
              zdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
                   &                  ((ir3-1)/DBLE(nr3sx))
              !         
           END DO
        END DO
#if defined(__MPI)
     END IF
#endif
  END DO
  RETURN
END SUBROUTINE grid_map
!
!----------------------------------------------------------------------------
SUBROUTINE setwfg( ibrav, b1, b2, b3 )
  !----------------------------------------------------------------------------
  !
  ! ... added by Young-Su Lee ( Nov 2006 )
  ! Find G vectors for a given ibrav and celldms
  !
  USE kinds,              ONLY : DP
  USE cell_base,          ONLY : tpiba, celldm
  USE wannier_base,       ONLY : wfg, nw, weight
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3)
  INTEGER,  INTENT(IN) :: ibrav
  REAL(DP) :: tweight(6), t0, t1, t2, t3, t4, t5, t6
  INTEGER  :: twfg(6,3), kk
  
  twfg(:,:) = 0

  twfg(1,1)=1   
  twfg(1,2)=0    
  twfg(1,3)=0   

  twfg(2,1)=0    
  twfg(2,2)=1   
  twfg(2,3)=0   

  twfg(3,1)=0   
  twfg(3,2)=0   
  twfg(3,3)=1   

  SELECT CASE(ibrav)

  CASE(1)
     !
     !  Cubic P [sc]
     !
     nw = 3
     !
     !
  CASE(2)
     !
     !  Cubic F [fcc]
     !
     nw = 4

     twfg(4,1)=-1    
     twfg(4,2)=-1   
     twfg(4,3)=-1  
     !
     !
  CASE(3)
     !
     !  Cubic I [bcc]
     !
     nw = 6

     twfg(4,1)=1 
     twfg(4,2)=1 
     twfg(4,3)=0   

     twfg(5,1)=0 
     twfg(5,2)=1 
     twfg(5,3)=1 

     twfg(6,1)=-1
     twfg(6,2)=0 
     twfg(6,3)=1 
     !
     !
  CASE(4)
     !
     !  Hexagonal and Trigonal P
     !
     nw = 4

     twfg(4,1)=1     
     twfg(4,2)=-1   
     twfg(4,3)=0   
     !
     !
  CASE(5)
     !
     !  Trigonal R
     !
     t0 = 1.D0/3.D0
     !
     IF ( celldm(4) .ge. t0 ) THEN    
        !
        nw = 4
        !
        twfg(4,1)=1
        twfg(4,2)=1
        twfg(4,3)=1
        !
     ELSE
        !
        IF ( celldm(4) .gt. 0 ) THEN
           !
           nw = 6
           ! 
           twfg(4,1)=1
           twfg(4,2)=1
           twfg(4,3)=0
         
           twfg(5,1)=0
           twfg(5,2)=1
           twfg(5,3)=1
         
           twfg(6,1)=1
           twfg(6,2)=0
           twfg(6,3)=1
           !
        ELSE IF ( celldm(4) .eq. 0 ) THEN
           !
           nw = 3
           !
        ELSE
           !
           nw = 6
           ! 
           twfg(4,1)=1
           twfg(4,2)=-1
           twfg(4,3)=0
         
           twfg(5,1)=0
           twfg(5,2)=1
           twfg(5,3)=-1
         
           twfg(6,1)=-1
           twfg(6,2)=0
           twfg(6,3)=1
           !
        END IF
        !
     END IF

  CASE(6)
     !
     !  Tetragonal P [st]
     !
     nw = 3
     !
     !
  CASE(7)
     !
     !  Tetragonal I [bct]
     !
     nw = 6

     twfg(4,1)=1  
     twfg(4,2)=0  
     twfg(4,3)=1  

     twfg(5,1)=0  
     twfg(5,2)=1  
     twfg(5,3)=-1 

     twfg(6,1)=1  
     twfg(6,2)=1  
     twfg(6,3)=0  
     !
     !
  CASE(8)
     !
     !  Orthorhombic P
     !
     nw = 3
     !
     !
  CASE(9)
     !
     !  Orthorhombic C
     !
     IF (celldm(2).EQ.1) THEN   ! Tetragonal P 
        !
        nw=3
        !
     ELSE
        !
        nw = 4
        !
        IF ( celldm(2) < 1 ) THEN
           !
           twfg(4,1)=1
           twfg(4,2)=-1
           twfg(4,3)=0
           !
        ELSE
           !
           twfg(4,1)=1
           twfg(4,2)=1
           twfg(4,3)=0
           !
        END IF
        !
     END IF
     !
     !
  CASE(10)
     !
     !  Orthorhombic F
     !
     twfg(4,1)=1    
     twfg(4,2)=1   
     twfg(4,3)=1   
     !
     IF ( celldm(2) .eq. 1 .AND. celldm(3) .eq. 1 )  THEN ! Cubic F
        !
        nw = 4
        !
     ELSE
        !
        nw = 6
        !
        IF ( celldm(2) .eq. 1 .AND. celldm(3) .ne. 1) THEN ! Tetragonal I
           !
           twfg(5,1)=1  
           twfg(5,2)=1  
           twfg(5,3)=0  
           twfg(6,1)=0  
           twfg(6,2)=1  
           twfg(6,3)=1  
           !
        ELSE IF ( celldm(2) .ne. 1 .AND. celldm(3) .eq. 1) THEN ! Tetragonal I
           !
           twfg(5,1)=1   
           twfg(5,2)=1   
           twfg(5,3)=0   
           twfg(6,1)=1   
           twfg(6,2)=0   
           twfg(6,3)=1   
           !
        ELSE IF ( celldm(2) .eq. celldm(3) ) THEN ! Tetragonal I
           !
           twfg(5,1)=0   
           twfg(5,2)=1   
           twfg(5,3)=1   
           twfg(6,1)=1   
           twfg(6,2)=0   
           twfg(6,3)=1   
           !
        ELSE IF ( celldm(2) .gt. 1 .and. celldm(3) .gt. 1 ) THEN
           !
           twfg(5,1)=0   
           twfg(5,2)=1   
           twfg(5,3)=1   
           twfg(6,1)=1   
           twfg(6,2)=0   
           twfg(6,3)=1   
           !
        ELSE IF ( celldm(2) .lt. celldm(3) ) THEN
           !
           twfg(5,1)=1   
           twfg(5,2)=1   
           twfg(5,3)=0   
           twfg(6,1)=1   
           twfg(6,2)=0   
           twfg(6,3)=1   
           !
        ELSE           
           !
           twfg(5,1)=1   
           twfg(5,2)=1   
           twfg(5,3)=0   
           twfg(6,1)=0   
           twfg(6,2)=1   
           twfg(6,3)=1   
           !
        END IF
        !
     END IF
     !
     !
  CASE(11) 
     !
     !  Orthorhombic I
     !
     nw = 6
     !
     twfg(4,1)=1     
     twfg(4,2)=1     
     twfg(4,3)=0     

     twfg(5,1)=0     
     twfg(5,2)=1     
     twfg(5,3)=1     

     twfg(6,1)=-1   
     twfg(6,2)=0     
     twfg(6,3)=1     
     !
     !
  CASE(12)
     !
     !  Monoclinic P
     !
     IF ( celldm(4) .eq. 0 ) THEN ! Orthorhombic P
        !
        nw = 3
        !
     ELSE
        !
        nw = 4
        !
        t1 = SQRT(DOT_PRODUCT(b1,b1))
        t2 = SQRT(DOT_PRODUCT(b2,b2))
        t4 = DOT_PRODUCT(b1,b2)/t1/t2
        !
        t0 = - t4 * t1 / t2
        kk = NINT(t0)
        !
        IF((kk.EQ.0).AND.(t0.GT.0)) kk=1
        IF((kk.EQ.0).AND.(t0.LT.0)) kk=-1
      
        twfg(4,1)=1     
        twfg(4,2)=kk    
        twfg(4,3)=0     
        !
     END IF
     !
     !
  CASE(0,13,14)
     !
     !  Monoclinic C, Triclinic P, Free Cell
     !
     nw = 6 
     !
     t1 = SQRT(DOT_PRODUCT(b1,b1))
     t2 = SQRT(DOT_PRODUCT(b2,b2))
     t3 = SQRT(DOT_PRODUCT(b3,b3))
     t4 = DOT_PRODUCT(b1,b2)/t1/t2
     t5 = DOT_PRODUCT(b2,b3)/t2/t3
     t6 = DOT_PRODUCT(b3,b1)/t3/t1
     !
     t0 = - t4 * t1 / t2
     kk = NINT(t0)
     !
     IF((kk.EQ.0).AND.(t0.GE.0)) kk=1
     IF((kk.EQ.0).AND.(t0.LT.0)) kk=-1
     
     twfg(4,1)=1 
     twfg(4,2)=kk
     twfg(4,3)=0 
     !
     t0 = - t5 * t2 / t3
     kk = NINT(t0)
     !
     IF((kk.EQ.0).AND.(t0.GE.0)) kk=1
     IF((kk.EQ.0).AND.(t0.LT.0)) kk=-1
     !
     twfg(5,1)=0  
     twfg(5,2)=1   
     twfg(5,3)=kk 
     !
     t0 = - t6 * t3 / t1
     kk = NINT(t0)
     !
     IF((kk.EQ.0).AND.(t0.GE.0)) kk=1
     IF((kk.EQ.0).AND.(t0.LT.0)) kk=-1
     ! 
     twfg(6,1)=kk 
     twfg(6,2)=0 
     twfg(6,3)=1
     !   
  CASE DEFAULT
     ! alternate ibrav values -5, -9, -12: not sure what to do 
     !
     CALL errore ('setwfc','ibrav not implemented', ABS(ibrav))
     !
  END SELECT
  !
  CALL tric_wts2( b1, b2, b3, nw, twfg, tweight ) 
  !
  ALLOCATE(wfg(nw,3), weight(nw))
  !
  wfg(:,:) = twfg(1:nw,:)
  weight(:) = tweight(1:nw)
  !
  RETURN
  !
END SUBROUTINE setwfg
!
!----------------------------------------------------------------------------
SUBROUTINE tric_wts( rp1, rp2, rp3, alat, wts )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine computes the weights to be used for
  ! ... R.P. translations in the WF calculation in the case
  ! ... of ibrav=0 or ibrav=14
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  USE cell_base, ONLY : tpiba, tpiba2
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: rp1(3), rp2(3), rp3(3)
  REAL(DP), INTENT(IN)  :: alat
  REAL(DP), INTENT(OUT) :: wts(6)
  ! 
  REAL(DP) :: b1x, b2x, b3x, b1y, b2y, b3y, b1z, b2z, b3z
  INTEGER        :: i
  !
  !
  b1x = rp1(1)*tpiba
  b2x = rp2(1)*tpiba
  b3x = rp3(1)*tpiba
  b1y = rp1(2)*tpiba
  b2y = rp2(2)*tpiba
  b3y = rp3(2)*tpiba
  b1z = rp1(3)*tpiba
  b2z = rp2(3)*tpiba
  b3z = rp3(3)*tpiba
  !        WRITE( stdout, * ) 'COMPUTING WEIGHTS NOW ...'

  wts(1) = tpiba2*(-b1z*b2x*b2z*b3x + b2y**2*b3x**2 + b1z*b2z*b3x**2 + & 
       b2z**2*b3x**2 - b1z*b2y*b2z*b3y - 2.D0*b2x*b2y*b3x*b3y + & 
       b2x**2*b3y**2 + b1z*b2z*b3y**2 + b2z**2*b3y**2 + & 
       b1z*b2x**2*b3z + b1z*b2y**2*b3z - b1z*b2x*b3x*b3z - & 
       2.D0*b2x*b2z*b3x*b3z - b1z*b2y*b3y*b3z - &
       2.D0*b2y*b2z*b3y*b3z + b2x**2*b3z**2 + b2y**2*b3z**2 + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2y*(b2x + b3x)*b3y - &
       b2z*(b2x + b3x)*b3z +  b2x*(b3y**2 + b3z**2)) + &
       b1y*(b2x**2*b3y - b2x*b3x*(b2y + b3y) + &
       b2z*b3y*(b2z - b3z) + b2y*(b3x**2 - b2z*b3z +  b3z**2)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y &
       + b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(2) = tpiba2*(b1z**2*(b2x*b3x + b3x**2 + b3y*(b2y + b3y)) + &
       b1y**2*(b2x*b3x + b3x**2 + b3z*(b2z + b3z)) - &
       b1z*(-b2z*(b3x**2 + b3y**2) + (b2x*b3x + b2y*b3y)*b3z + &
       b1x*(b2z*b3x + (b2x + 2.D0* b3x)*b3z)) - &
       b1y*(b1x*(b2y*b3x + (b2x + 2.D0*b3x)*b3y) + &
       b3y*(b1z*b2z + b2x*b3x + 2.D0*b1z*b3z + b2z*b3z) - &
       b2y*(b3x**2 - b1z* b3z + b3z**2)) + &
       b1x*(-b2y*b3x*b3y + b2x*b3y**2 - b2z*b3x*b3z + b2x*b3z**2 + &
       b1x*(b2y*b3y + b3y**2 + b3z*(b2z + b3z))))/ &                        
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(3) = tpiba2*(b1z**2*(b2x**2 + b2x*b3x + b2y*(b2y + b3y)) - &
       b1y*(2.D0*b1z*b2y*b2z + b2x*b2y*b3x - b2x**2*b3y + &
       b1z*b2z*b3y - b2z**2*b3y + b1x*(2.D0*b2x*b2y + b2y*b3x + b2x*b3y) + &
       b1z*b2y*b3z + b2y*b2z*b3z) + b1y**2*(b2x**2 + b2x*b3x + b2z*(b2z + b3z)) - &
       b1z*(b2x*b2z*b3x + b2y*b2z*b3y - b2x**2*b3z - b2y**2*b3z + &
       b1x*(2.D0*b2x*b2z + b2z*b3x + b2x*b3z)) + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2x*b2y*b3y - b2x*b2z*b3z + &
       b1x*(b2y**2 + b2y*b3y + b2z*(b2z +   b3z))))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + b1y*b2x*b3z - &
       b1x*b2y*b3z)**2) 

  wts(4) = tpiba2*(b1z*(-b2z*(b3x**2 + b3y**2) + (b2x*b3x + b2y*b3y)*b3z) + & 
       b1y*(b3y*(b2x*b3x + b2z*b3z) - b2y*(b3x**2 + b3z**2)) + &
       b1x*(b2y*b3x*b3y + b2z*b3x*b3z - b2x*(b3y**2 +  b3z**2)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(5) =  -tpiba2*(b1z**2*(b2x*b3x + b2y*b3y) - b1x*b1z*(b2z*b3x + b2x*b3z) - &
       b1y*(b1x*b2y*b3x + b1x*b2x*b3y + b1z*b2z*b3y + b1z*b2y*b3z) + &
       b1y**2*(b2x*b3x + b2z*b3z) + b1x**2*(b2y*b3y + b2z*b3z))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(6) = -tpiba2*(b1z*(-b2x*b2z*b3x - b2y*b2z*b3y + b2x**2*b3z + b2y**2*b3z) + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2x*b2y*b3y -  b2x*b2z*b3z) + &
       b1y*(-b2x*b2y*b3x + b2x**2*b3y + b2z*(b2z*b3y -  b2y*b3z)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)
  !
  RETURN
  !
END SUBROUTINE tric_wts
!
!----------------------------------------------------------------------------
SUBROUTINE tric_wts2( rp1, rp2, rp3, nw, wfg, weight )
  !----------------------------------------------------------------------------
  !
  ! ... added by Young-Su Lee ( Nov 2006 )
  !
  ! Find the least square solutions of weights for G vectors
  ! If the set of G vectors and calculated weights do not conform to the condition,
  !  SUM_i weight_i G_ia G_ib = delta_ab
  ! the code stops.
  !
  USE kinds,              ONLY : DP
  USE io_global,          ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: rp1(3), rp2(3), rp3(3)
  INTEGER, INTENT(IN)   :: wfg(6,3), nw
  REAL(DP), INTENT(OUT) :: weight(6)
  !
  REAL(DP) :: gp(6,nw), A(6,nw), gr(nw,3), S(6), R(6), WORK(1000), t
  INTEGER  :: i, LWORK, INFO
  !
  DO i=1, nw
     gr(i,:) = wfg(i,1)*rp1(:)+wfg(i,2)*rp2(:)+wfg(i,3)*rp3(:)
  END DO
                                                                                                                       
  DO i=1, nw
     gp(1,i)=gr(i,1)*gr(i,1)
     gp(2,i)=gr(i,2)*gr(i,2)
     gp(3,i)=gr(i,3)*gr(i,3)
     gp(4,i)=gr(i,1)*gr(i,2)
     gp(5,i)=gr(i,2)*gr(i,3)
     gp(6,i)=gr(i,3)*gr(i,1)
  END DO
  !
  R = 0.D0
  R(1:3) = 1.D0
  !
  LWORK=1000
  A = gp
  S = R
  !
  CALL DGELS( 'N', 6, nw, 1, A, 6, S, 6, WORK, LWORK, INFO )
  !
  IF (INFO .ne. 0) THEN
     WRITE( stdout, * ) "failed to get a weight factor for ",INFO,"th vector"
     STOP
  END IF
  !
  weight(1:nw) = S(1:nw)
  S=matmul(gp,weight(1:nw))
  !
  DO i=1, nw
     IF ( weight(i) .lt. 0.D0 ) THEN
        WRITE( stdout, * ) "WARNING: weight factor less than zero"
     END IF
  END DO
  !
  DO i=1,6
     t = abs(S(i)-R(i))
     IF ( t .gt. 1.D-8 ) THEN
        WRITE( stdout, * ) "G vectors do not satisfy the completeness condition",i,t
        !HK
        WRITE( stdout, * ) "NOTE: If you are in a variable cell calculation, "
        WRITE( stdout, * ) " please take a look if your ibrav is consistent with the restarted cell"
        WRITE( stdout, * ) " use ibrav=0 and CELL_PARAMETERS card."
        STOP
     END IF
  END DO
  !
  RETURN
  !
END SUBROUTINE tric_wts2
!
!----------------------------------------------------------------------------
SUBROUTINE small_box_wf( i_1, j_1, k_1, nw1 )
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE constants,              ONLY : fpi
  USE wannier_base,           ONLY : expo
  USE fft_base,               ONLY : dfftp
  USE mp_global,              ONLY : me_bgrp
  USE parallel_include
  !
  IMPLICIT NONE

  INTEGER ir1, ir2, ir3, ibig3 , inw
  REAL(DP) x
  INTEGER , INTENT(in) :: nw1, i_1(nw1), j_1(nw1), k_1(nw1)
  INTEGER :: me

  me = me_bgrp + 1

  ALLOCATE(expo(dfftp%nnr,nw1))

  DO inw=1,nw1

     WRITE( stdout, * ) inw ,":", i_1(inw), j_1(inw), k_1(inw)

     DO ir3=1,dfftp%nr3
#if defined(__MPI)
        ibig3 = ir3 - dfftp%ipp( me )
        IF(ibig3.GT.0.AND.ibig3.LE.dfftp%npp(me)) THEN
#else
           ibig3=ir3
#endif
           DO ir2=1,dfftp%nr2
              DO ir1=1,dfftp%nr1
                 x =  (((ir1-1)/DBLE(dfftp%nr1x))*i_1(inw) +                          &
                      &                  ((ir2-1)/DBLE(dfftp%nr2x))*j_1(inw) +             &
                      &                  ((ir3-1)/DBLE(dfftp%nr3x))*k_1(inw))*0.5d0*fpi
                 expo(ir1+(ir2-1)*dfftp%nr1x+(ibig3-1)*dfftp%nr1x*dfftp%nr2x,inw) = CMPLX(COS(x), -SIN(x),kind=DP)
              END DO
           END DO
#if defined(__MPI)
        END IF
#endif
     END DO
  END DO
  RETURN
END SUBROUTINE small_box_wf
!
!-----------------------------------------------------------------------
FUNCTION boxdotgridcplx(irb,qv,vr)
  !-----------------------------------------------------------------------
  !
  ! Calculate \sum_i qv(r_i)*vr(r_i)  with r_i on box grid
  ! array qv(r) is defined on box grid, array vr(r)on dense grid
  ! irb   : position of the box in the dense grid
  ! Parallel execution: remember to sum the contributions from other nodes
  !
  !      use ion_parameters
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp, dfftb
  USE mp_global,       ONLY : me_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER,           INTENT(IN):: irb(3)
  COMPLEX(DP), INTENT(IN):: qv(dfftb%nnr), vr(dfftp%nnr)
  COMPLEX(DP)            :: boxdotgridcplx
  !
  INTEGER :: ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig, me
  !
  me = me_bgrp + 1
  !
  boxdotgridcplx = ZERO

  DO ir3=1,dfftb%nr3
     ibig3=irb(3)+ir3-1
     ibig3=1+MOD(ibig3-1,dfftp%nr3)
#if defined(__MPI)
     ibig3 = ibig3 - dfftp%ipp( me )
     IF (ibig3.GT.0.AND.ibig3.LE.dfftp%npp(me)) THEN
#endif
        DO ir2=1,dfftb%nr2
           ibig2=irb(2)+ir2-1
           ibig2=1+MOD(ibig2-1,dfftp%nr2)
           DO ir1=1,dfftb%nr1
              ibig1=irb(1)+ir1-1
              ibig1=1+MOD(ibig1-1,dfftp%nr1)
              ibig=ibig1 + (ibig2-1)*dfftp%nr1x + (ibig3-1)*dfftp%nr1x*dfftp%nr2x
              ir  =ir1 + (ir2-1)*dfftb%nr1x + (ir3-1)*dfftb%nr1x*dfftb%nr2x
              boxdotgridcplx = boxdotgridcplx + qv(ir)*vr(ibig)
           END DO
        END DO
#if defined(__MPI)
     ENDIF
#endif
  END DO
  !
  RETURN
  !
END FUNCTION boxdotgridcplx
!
!----------------------------------------------------------------------------
SUBROUTINE write_rho_g( rhog )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE io_global,          ONLY : stdout
  USE gvect,              ONLY : ngm
  USE gvect, ONLY : g
  USE electrons_base,     ONLY : nspin
  USE fft_base,           ONLY : dfftp
  USE mp_global,          ONLY : nproc_bgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  USE mp,                 ONLY : mp_barrier, mp_gather, mp_set_displs
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) ,INTENT(IN) :: rhog(ngm,nspin) 
  REAL(DP),   ALLOCATABLE:: gnx(:,:), bigg(:,:)
  COMPLEX(DP),ALLOCATABLE :: bigrho(:)
  COMPLEX(DP) :: rhotmp_g(ngm)
  INTEGER           :: ntot, i, j, me
#if defined(__MPI)
  INTEGER proc, ierr, ngdens(nproc_bgrp), displs(nproc_bgrp)
#endif
  CHARACTER (LEN=6)  :: name
  CHARACTER (LEN=15) :: name2

  me = me_bgrp + 1

  ALLOCATE(gnx(3,ngm))

  DO i=1,ngm
     gnx(1,i)=g(1,i)
     gnx(2,i)=g(2,i)
     gnx(3,i)=g(3,i)
  END DO

#if defined(__MPI)

  DO i=1,nproc_bgrp
     ngdens(i)=(dfftp%ngl(i)+1)/2
  END DO

  CALL mp_set_displs( ngdens, displs, ntot, nproc_bgrp )

  IF(me.EQ.1) THEN
     ALLOCATE(bigg(3,ntot))
  END IF

  CALL mp_barrier(intra_bgrp_comm)
 
  CALL mp_gather( gnx, bigg, ngdens, displs,  root_bgrp, intra_bgrp_comm )

  DO i=1,nspin

     rhotmp_g(1:ngm)=rhog(1:ngm,i)

     IF(me.EQ.1) THEN 
        ALLOCATE (bigrho(ntot))
     END IF

     CALL mp_barrier(intra_bgrp_comm)

     CALL mp_gather( rhotmp_g, bigrho, ngdens, displs,  root_bgrp, intra_bgrp_comm )

     IF(me.EQ.1) THEN
        IF(i.EQ.1) name2="CH_DEN_G_PARA.1"
        IF(i.EQ.2) name2="CH_DEN_G_PARA.2"
        OPEN(unit=57, file=name2) 
        DO j=1,ntot
           WRITE(57,*) bigrho(j)
        END DO
        CLOSE(57)
        DEALLOCATE(bigrho)
     END IF

     WRITE( stdout, * ) "Charge density written to ", name2

  END DO

  IF(me.EQ.1) THEN
     name="G_PARA"
     OPEN(unit=56, file=name) 
     DO i=1,ntot
        WRITE(56,*) bigg(:,i)
     END DO
     CLOSE(56)
     DEALLOCATE(bigg)
  END IF
  WRITE( stdout, * ) "G-vectors written to G_PARA"
#else
  ntot=ngm
  ALLOCATE(bigg(3,ntot))
  bigg(1:3,1:ntot)=gnx(1:3,1:ngm)
  DO i=1,nspin
     ALLOCATE(bigrho(ntot))
     bigrho(1:ngm)=rhog(1:ngm,i)

     IF(i.EQ.1) name2="CH_DEN_G_SERL.1"
     IF(i.EQ.2) name2="CH_DEN_G_SERL.2"

     OPEN(unit=57, file=name2) 
     DO j=1,ntot
        WRITE(57,*) bigrho(j)
     END DO
     CLOSE(57)
     DEALLOCATE(bigrho)

     WRITE( stdout, * ) "Charge density written to", name2

  END DO

  name="G_SERL"
  OPEN(unit=56, file=name) 
  DO i=1,ntot
     WRITE(56,*) bigg(:,i)
  END DO
  CLOSE(56)
  DEALLOCATE(bigg)
  WRITE( stdout, * ) "G-vectors written to G_SERL"
#endif
  !
  DEALLOCATE(gnx)
  !
  RETURN
  !
END SUBROUTINE write_rho_g
!
!----------------------------------------------------------------------------
SUBROUTINE macroscopic_average( rhog, tau0, e_tuned )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE gvect, ONLY : g
  USE gvect,              ONLY : ngm
  USE electrons_base,     ONLY : nspin
  USE tune,               ONLY : npts, xdir, ydir, zdir, B, &
                                 shift, start, av0, av1
  USE cell_base,          ONLY : at, alat, tpiba, omega
  USE ions_base,          ONLY : nsp, na, zv, nax
  USE constants,          ONLY : pi, tpi
  USE mp,                 ONLY : mp_barrier, mp_bcast,  mp_gather, mp_set_displs
  USE fft_base,           ONLY : dfftp
  USE mp_global,          ONLY : nproc_bgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE:: gnx(:,:), bigg(:,:)
  COMPLEX(DP) ,INTENT(in) :: rhog(ngm,nspin)
  COMPLEX(DP),ALLOCATABLE :: bigrho(:)
  COMPLEX(DP), ALLOCATABLE :: rhotmp_g(:)
  INTEGER ntot, i, j, ngz, l, isa
  INTEGER ,ALLOCATABLE :: g_red(:,:)
#if defined(__MPI)
  INTEGER proc, ierr, ngdens(nproc_bgrp), displs( nproc_bgrp )
#endif
  REAL(DP) zlen,vtot, pos(3,nax,nsp), a_direct(3,3),a_trans(3,3), e_slp, e_int
  REAL(DP), INTENT(out) :: e_tuned(3)
  REAL(DP), INTENT(in) :: tau0(3,nax)
  REAL(DP),ALLOCATABLE :: v_mr(:), dz(:), gz(:), g_1(:,:), vbar(:), cd(:), v_final(:)
  REAL(DP), ALLOCATABLE:: cdion(:), cdel(:), v_line(:), dist(:)
  COMPLEX(DP),ALLOCATABLE :: rho_ion(:),v_1(:),vmac(:),rho_tot(:),rhogz(:), bigrhog(:)
  INTEGER :: me

  me = me_bgrp + 1

  ALLOCATE(gnx(3,ngm))

  DO i=1,ngm
     gnx(1,i)=g(1,i)
     gnx(2,i)=g(2,i)
     gnx(3,i)=g(3,i)
  END DO

#if defined(__MPI)

  DO i=1,nproc_bgrp
     ngdens(i)=(dfftp%ngl(i)+1)/2
  END DO

  CALL mp_set_displs( ngdens, displs, ntot, nproc_bgrp )

#else

  ntot=ngm

#endif

  ALLOCATE(bigg(3,ntot))
  ALLOCATE (bigrho(ntot))
  ALLOCATE (bigrhog(2*ntot-1))

#if defined(__MPI)
  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_gather( gnx, bigg, ngdens, displs, root_bgrp,intra_bgrp_comm )
  !
  CALL mp_bcast( bigg, root_bgrp, intra_bgrp_comm )
  !
  ALLOCATE( rhotmp_g( ngm ) )

  rhotmp_g(1:ngm)=rhog(1:ngm,1)

  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_gather( rhotmp_g, bigrho, ngdens, displs, root_bgrp,intra_bgrp_comm )
  !
  DEALLOCATE( rhotmp_g )
  !
  CALL mp_bcast( bigrho, root_bgrp, intra_bgrp_comm )
  !
#else
  !
  bigg(1:3,1:ntot)=gnx(1:3,1:ngm)
  bigrho(1:ngm)=rhog(1:ngm,1)
  !
#endif

  ALLOCATE(g_1(3,2*ntot-1))
  ALLOCATE(g_red(3,2*ntot-1))

  ALLOCATE(v_mr(npts))
  ALLOCATE(v_final(npts))
  ALLOCATE(dz(npts))
  ALLOCATE(vbar(npts))
  ALLOCATE(cd(npts))
  ALLOCATE(cdel(npts))
  ALLOCATE(cdion(npts))

  !-- needed for non-orthogonal cells

  a_direct(1,1:3)=at(1:3,1)*alat
  a_direct(2,1:3)=at(1:3,2)*alat
  a_direct(3,1:3)=at(1:3,3)*alat

  a_trans=TRANSPOSE(a_direct)

  !--- Construct rho(-g) from rho(g). rgo(-g)=rho*(g)

  bigrhog(1:ntot)=bigrho(1:ntot)
  g_1(:,1:ntot)=bigg(:,1:ntot)
  DO i=2,ntot
     bigrhog(ntot+i-1)=CONJG(bigrho(i))
     g_1(:,ntot+i-1)=-bigg(:,i)
  END DO

  !--- needed fot non-orthogonal cells

  DO i=1,2*ntot-1
     g_red(:,i)=NINT(MATMUL(a_trans(:,:),g_1(:,i))*tpiba/tpi)
  END DO

  !--- define the direction of the line

  xdir=1
  ydir=2

  IF ((zdir).EQ.1) xdir=3
  IF ((zdir).EQ.2) ydir=3

  IF(zdir.EQ.1) zlen=DSQRT(at(1,1)**2+at(2,1)**2+at(3,1)**2)
  IF(zdir.EQ.2) zlen=DSQRT(at(1,2)**2+at(2,2)**2+at(3,2)**2)
  IF(zdir.EQ.3) zlen=DSQRT(at(1,3)**2+at(2,3)**2+at(3,3)**2)
  zlen = zlen*alat

  !--- We need the potentiail only along zdir, so pick the appropriate G-vectors with Gxdir=Gydir=0

  ngz=0
  DO i=1,2*ntot-1
     IF((g_red(xdir,i).EQ.0).AND.(g_red(ydir,i).EQ.0)) ngz=ngz+1
  END DO

  ALLOCATE(gz(ngz))
  ALLOCATE(rhogz(ngz))
  ALLOCATE(rho_ion(ngz))
  ALLOCATE(rho_tot(ngz))
  ALLOCATE(vmac(ngz))
  ALLOCATE(v_1(ngz))

  !--- The G-vectors are output in units of 2*pi/a, so convert them to the correct values

  j=0
  DO i=1,2*ntot-1
     IF((g_red(xdir,i).EQ.0).AND.(g_red(ydir,i).EQ.0)) THEN
        j=j+1
        gz(j)=g_1(zdir,i)*tpiba
        rhogz(j)=bigrhog(i)
     END IF
  END DO

  isa = 0
  DO i=1,nsp
     DO j=1,na(i)
        isa = isa + 1
        pos(:,j,i)=tau0(:,isa)
     END DO
  END DO

  !--- Construct the ionic Charge density in G-space

  rho_ion = ZERO
  !
  DO j=1,ngz
     DO i=1,nsp
        DO l=1,na(i)
           rho_ion(j)=rho_ion(j)+zv(i)*EXP(-CI*gz(j)*pos(zdir,l,i))*EXP(-gz(j)**2/(4.D0*ONE))
        END DO
     END DO
  END DO

  rho_ion=rho_ion/omega

  !--- Construct the total Charge density in G-space

  rho_tot=rho_ion-rhogz

  !--- Construct the electrostatic potential and macroscopic average in G-space

  v_1(1)=ZERO
  vmac(1)=ZERO
  v_1(2:ngz)=4*pi*rho_tot(2:ngz)/gz(2:ngz)**2
  vmac(2:)=v_1(2:)*SIN(gz(2:)*b)/(gz(2:)*b)


  !--- Calculate planewise average in R-space and FFT V(Gz) ---> V(z) ... well not really FFT but FT

  vbar=0.D0
  v_mr=0.D0
  cdel=0.D0
  cdion=0.D0
  cd=0.D0
  DO j=1,npts
     dz(j)=(j-1)*zlen/(npts*1.D0)
     DO i=1,ngz
        vbar(j)=vbar(j)-DBLE(EXP(CI*gz(i)*dz(j))*v_1(i))
        v_mr(j)=v_mr(j)-DBLE(EXP(CI*gz(i)*dz(j))*vmac(i))
        cdel(j)=cdel(j)-DBLE(EXP(CI*gz(i)*dz(j))*rhogz(i))
        cdion(j)=cdion(j)+DBLE(EXP(CI*gz(i)*dz(j))*rho_ion(i))
        cd(j)=cd(j)+DBLE(EXP(CI*gz(i)*dz(j))*rho_tot(i))
     END DO
     !           WRITE( stdout, * ) vbar(j), v_mr(j), cdel(j), cdion(j)
  END DO
  IF (shift) THEN
     vtot=(v_mr(start)+v_mr(start-1))/2.D0
     v_final(1:npts-start+1)=v_mr(start:npts)-vtot
     v_final(npts-start+2:npts)=v_mr(1:start-1)-vtot
  ELSE
     vtot=(v_mr(1)+v_mr(npts))/2.D0
     v_final(1:npts)=v_mr(1:npts)-vtot
  END IF

  e_tuned=0.D0

  ALLOCATE(v_line(1:av1-av0+1))
  ALLOCATE(dist(1:av1-av0+1))


  v_line(1:av1-av0+1)=v_final(av0:av1)
  dist(1:av1-av0+1) =dz(av0:av1)

  e_tuned(zdir)=-(v_final(av1)-v_final(av0))/((av1-av0)*zlen/(npts*1.D0))


  DEALLOCATE(bigg,g_1,bigrho,bigrhog,g_red)     
  DEALLOCATE(gnx,v_mr,v_final,dz,vbar,cd,cdel,cdion)
  DEALLOCATE(v_line, dist)
  DEALLOCATE(gz,rhogz,rho_ion,rho_tot,vmac,v_1)

  RETURN
END SUBROUTINE macroscopic_average
!
!----------------------------------------------------------------------------
SUBROUTINE least_square( npts, x, y, slope, intercept )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN) :: npts
  REAL(DP), INTENT(IN) :: x(npts), y(npts)
  REAL(DP), INTENT(OUT):: slope, intercept
  !
  INTEGER        :: i
  REAL(DP) :: sumx,sumy,sumx2,sumxy,sumsqx
  REAL(DP) :: xav,yav

  sumxy=0.D0
  sumx =0.D0
  sumy =0.D0
  sumx2=0.D0
  DO i=1,npts
     sumxy=sumxy+x(i)*y(i)
     sumx =sumx +x(i)
     sumy =sumy +y(i)
     sumx2=sumx2+x(i)*x(i)
  END DO
  sumsqx=sumx**2
  xav=sumx/DBLE(npts)
  yav=sumy/DBLE(npts)

  slope=(npts*sumxy - sumx*sumy)/(npts*sumx2 - sumsqx)

  intercept=yav-slope*xav

  RETURN

END SUBROUTINE least_square
!
!----------------------------------------------------------------------------
SUBROUTINE wfsteep( m, Omat, Umat, b1, b2, b3 )
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout, ionode
  USE wannier_base,           ONLY : nw, weight, nit, tolw, wfdt, maxwfdt, nsd
  USE control_flags,          ONLY : iverbosity,conv_elec
  USE cell_base,              ONLY : alat
  USE constants,              ONLY : tpi, bohr_radius_angs
  USE mp_global,              ONLY : me_bgrp
  USE printout_base,          ONLY : printout_base_open, printout_base_unit, &
                                     printout_base_close
  USE cp_main_variables,      ONLY : nfi, iprint_stdout
  USE time_step,              ONLY : tps
  USE input_parameters,         ONLY : tcpbo
  USE parallel_include
  !
  IMPLICIT NONE

  !    (m,m) is the size of the matrix Ospin.
  !    Ospin is input overlap matrix.
  !    Uspin is the output unitary transformation.
  !             Rough guess for Uspin can be carried in.
  !
  !     conjugated gradient to search maximization
  !
  INTEGER, INTENT(in) :: m
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  COMPLEX(DP), INTENT(inout) :: Omat(nw, m, m)
  REAL(DP), INTENT(inout) :: Umat(m,m)
  !
  INTEGER :: i, j, k, l, ig, ierr, ti, tj, tk, inw, ir, adjust
  INTEGER :: f3(nw), f4(nw), ierr1
  REAL(DP) :: slope, slope2, t1, t2, t3, mt(nw),t21,temp1,maxdt
  REAL(DP) :: U(m,m), wfc(3, m), Wm(m,m), schd(m,m), f2(4*m), gr(nw, 3)
  REAL(DP) :: Uspin2(m,m),temp2,wfdtold,oldt1,t01, d3(m,m), d4(m,m), U1(m,m)
  REAL(DP) :: spread, sp
  REAL(DP), ALLOCATABLE  :: wr(:)
  REAL(DP), ALLOCATABLE  :: W(:,:)
  COMPLEX(DP) :: ct1, ct2, ct3, z(m, m), X(m, m), d(m,m), d2(m,m)
  COMPLEX(DP) :: f1(2*m-1), wp(m*(m+1)/2), Oc(nw, m, m)
  COMPLEX(DP) ::  Oc2(nw, m, m),wp1(m*(m+1)/2), X1(m,m), U2(m,m), U3(m,m)
  INTEGER :: me, iunit
  !
  me = me_bgrp + 1
  !
  ALLOCATE(W(m,m), wr(m))
  !
  Umat=0.D0
  DO i=1,m
     Umat(i,i)=1.D0
  END DO
  Oc=ZERO
  Oc2=ZERO
  X1=ZERO
  U2=Umat*ONE
  !
  ! update Oc using the initial guess of Uspin
  !
  DO inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=ZERO
     CALL zgemm ('T', 'N', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
     X1=ZERO
     CALL zgemm ('N','N', m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
     Oc(inw, :, :)=X1(:, :)
  END DO

  U2=ZERO
  U3=ZERO

  W=0.D0
  schd=0.D0
  oldt1=0.D0
  wfdtold=0.D0

  DO k=1, nit
     t01=0.D0     !use t1 to store the value of omiga
     DO inw=1, nw
        DO i=1, m
           t01=t01+DBLE(CONJG(Oc(inw, i, i))*Oc(inw, i, i))
        END DO
     END DO

     !    WRITE( stdout, * ) t01

     IF(ABS(oldt1-t01).LT.tolw) THEN 
        IF(me.EQ.1) THEN
           WRITE(stdout,*) "MLWF Generated at Step",k
        END IF
        IF( iverbosity > 2 ) THEN
           WRITE( stdout, * ) "MLWF Generated at Step",k
        END IF
        GO TO 40
     END IF

     !    oldt1=t01

     !   calculate d(omiga)/dW and store result in W
     !   W should be a real symmetric matrix for gamma-point calculation
     !
     Wm=W
     W=0.D0
     DO inw=1, nw
        t2=weight(inw)
        DO i=1,m
           DO j=i+1,m
              W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*CONJG(Oc(inw,i,i)        &
                   -Oc(inw,j,j))+CONJG(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
           END DO
        END DO
     END DO
     W=W-TRANSPOSE(W)

     !   calculate slope=d(omiga)/d(lamda)
     slope=SUM(W**2)

     !   calculate slope2=d2(omiga)/d(lamda)2
     slope2=0.D0
     DO ti=1, m
        DO tj=1, m
           DO tk=1, m
              t2=0.D0
              DO inw=1, nw
                 t2=t2+DBLE(Oc(inw,tj,tk)*CONJG(Oc(inw,tj,tj)+Oc(inw,tk,tk) &
                      -2.D0*Oc(inw,ti,ti))-4.D0*Oc(inw,ti,tk)          &
                      *CONJG(Oc(inw,ti,tj)))*weight(inw)
              END DO
              slope2=slope2+W(tk,ti)*W(ti,tj)*2.D0*t2
           END DO
        END DO
     END DO
     slope2=2.D0*slope2

     !   use parabola approximation. Defined by 1 point and 2 slopes
     IF (slope2.LT.0) wfdt=-slope/2.D0/slope2
     IF (maxwfdt.GT.0.AND.wfdt.GT.maxwfdt) wfdt=maxwfdt

     IF (k.LT.nsd) THEN
        schd=W    !use steepest-descent technique

        !   calculate slope=d(omiga)/d(lamda)
        slope=SUM(schd**2)

        !       schd=schd*maxwfdt
        DO i=1, m
           DO j=i, m
              wp1(i + (j-1)*j/2) = CMPLX(0.d0, schd(i,j),kind=DP)
           END DO
        END DO

        CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
        IF (ierr.NE.0) STOP 'failed to diagonalize W!'

     ELSE
        !
        CALL dgemm ('T','N', m,m,m,ONE,Wm,m,Wm,m,ZERO,d3,m)

        t1=0.D0
        DO i=1, m
           t1=t1+d3(i, i)
        END DO
        IF (t1.NE.0) THEN
           d4=(W-Wm)
           CALL dgemm ('T','N', m,m,m,ONE,W,m,d4,m,ZERO,d3,m)
           t2=0.D0
           DO i=1, m
              t2=t2+d3(i, i)
           END DO
           t3=t2/t1
           schd=W+schd*t3
        ELSE
           schd=W
        END IF
        !
        !   calculate the new d(Lambda) for the new Search Direction
        !   added by Manu. September 19, 2001
        !
        !   calculate slope=d(omiga)/d(lamda)
        slope=SUM(schd**2)
        !------------------------------------------------------------------------
        !   schd=schd*maxwfdt
        DO i=1, m
           DO j=i, m
              wp1(i + (j-1)*j/2) = CMPLX(0.d0, schd(i,j),kind=DP)
           END DO
        END DO

        CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
        IF (ierr.NE.0) STOP 'failed to diagonalize W!'

        maxdt=maxwfdt

11      d=0.D0
        DO i=1, m
           d(i, i)=EXP(CI*(maxwfdt)*wr(i))
        END DO

        U3=ZERO
        CALL zgemm ('N', 'N', m,m,m,ONE,z,m,d,m,ZERO,U3,m)
        U2=ZERO
        CALL zgemm ('N','C', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
        U=DBLE(U2)
        U2=ZERO
        U3=ZERO
        !
        !   update Uspin
        U1=ZERO
        CALL dgemm ('N', 'N', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)
        Umat=U1

        !
        !   update Oc
        !
        U2=Umat*ONE
        U3=ZERO
        DO inw=1, nw
           X1(:,:)=Omat(inw,:,:)
           CALL zgemm ('T', 'N', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
           X1=ZERO
           CALL zgemm ('N','N',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
           Oc2(inw, :,:)=X(:,:)
        END DO
        U2=ZERO 
        U3=ZERO
        !
        t21=0.D0     !use t21 to store the value of omiga
        DO inw=1, nw
           DO i=1, m
              t21=t21+DBLE(CONJG(Oc2(inw, i, i))*Oc2(inw, i, i))
           END DO
        END DO

        temp1=-((t01-t21)+slope*maxwfdt)/(maxwfdt**2)
        temp2=slope
        wfdt=-temp2/(2*temp1)

        IF (wfdt.GT.maxwfdt.OR.wfdt.LT.0.D0) THEN
           maxwfdt=2*maxwfdt
           GO TO 11
        END IF

        maxwfdt=maxdt
        !
        !
        !   use parabola approximation. Defined by 2 point and 1 slopes
        !    if (slope2.lt.0) wfdt=-slope/2.D0/slope2
        !    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
        !
        !    write(6, '(e12.5E2,1x,e11.5E2,1x,f6.2)') slope2, slope, wfdt
        !-------------------------------------------------------------------------
        !
        !      schd is the new searching direction
        !
     END IF

     d=0.D0
     DO i=1, m
        d(i, i)=EXP(CI*wfdt*wr(i))
     END DO          !d=exp(d)


     !   U=z*exp(d)*z+
     !
     U3=ZERO
     CALL zgemm ('N', 'N', m,m,m,ONE,z,m,d,m,ZERO,U3,m)
     U2=ZERO
     CALL zgemm ('N','C', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
     U=DBLE(U2)
     U2=ZERO
     U3=ZERO

     !   update Uspin
     !
     U1=ZERO
     CALL dgemm ('N', 'N', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)
     Umat=U1

     !   update Oc
     !
     U2=Umat*ONE
     U3=ZERO
     DO inw=1, nw
        X1(:, :)=Omat(inw, :, :)
        CALL zgemm ('T', 'N', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
        X1=ZERO
        CALL zgemm ('N','N',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
        Oc(inw, :, :)=X1(:, :)
     END DO
     U2=ZERO
     U3=ZERO
     IF(ABS(t01-oldt1).GE.tolw.AND.k.GE.nit) THEN
        IF(me.EQ.1) THEN
           WRITE(*,*) "MLWF Not generated after",k,"Steps."
        END IF
        IF( iverbosity > 2 ) THEN
           WRITE( stdout, * ) "MLWF Not generated after",k,"Steps."
        END IF
        GO TO 40
     END IF
     oldt1=t01
  END DO

40 DEALLOCATE(W, wr)

  !
  ! calculate the spread
  !
  !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
  IF( ( MOD( nfi, iprint_stdout ) == 0 ) )  THEN
     IF ( ionode ) THEN
        iunit = printout_base_unit( "spr" )
        CALL printout_base_open( "spr" )
        IF (.NOT.tcpbo) THEN
           WRITE( iunit, '(I8,f16.8)' )nfi,tps
        ELSE   
           IF(conv_elec) WRITE( iunit, '(I8,f16.8)' )nfi,tps
        END IF    
     END IF    
  END IF    

  spread=0.0_DP

  DO i=1, m
     !
     mt=1.D0-DBLE(Oc(:,i,i)*CONJG(Oc(:,i,i)))
     sp = (alat*bohr_radius_angs/tpi)**2*SUM(mt*weight)
     !
     !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
     IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
        IF ( ionode ) THEN
           IF (.NOT.tcpbo) THEN
              WRITE( iunit, '(f20.14)' ) sp
           ELSE
             IF(conv_elec) WRITE( iunit, '(f20.14)' ) sp
           END IF   
        END IF   
     END IF
     !
     IF( sp < 0.D0 ) &
        CALL errore( 'cp-wf', 'Something wrong WF Spread negative', 1 )
     !
     spread=spread+sp
     !
  END DO

  !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
  IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
     IF ( ionode ) CALL printout_base_close( "spr" )
  END IF

  spread=spread/DBLE(m)

  IF(me.EQ.1) THEN
     WRITE(24, '(f10.7)') spread
     WRITE(27,*) "Average spread = ", spread
  END IF
  !
  Omat=Oc
  !
  RETURN
END SUBROUTINE wfsteep
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE write_psi( c, jw )
  !----------------------------------------------------------------------------
  ! ... for calwf 5             - M.S
  ! ... collect wavefunctions on first node and write to file
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout, ionode
  USE gvecw ,                 ONLY : ngw
  USE electrons_base,         ONLY : nbspx
  USE mp,                     ONLY : mp_barrier, mp_set_displs, mp_gather
  USE fft_base,               ONLY : dfftp
  USE mp_global,              ONLY : nproc_bgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: jw
  COMPLEX(DP) :: c(ngw,nbspx)
  !
  INTEGER ::i, ig, proc, ntot, ngpwpp(nproc_bgrp)
  INTEGER ::displs(nproc_bgrp)
  COMPLEX(DP), ALLOCATABLE:: psitot(:)

#if defined (__MPI)
  !
  DO proc=1,nproc_bgrp
     ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
  END DO
  !
  CALL mp_set_displs( ngpwpp, displs, ntot, nproc_bgrp )
  !
  ! allocate the needed work spaces
  !
  IF ( me_bgrp == root_bgrp ) THEN
     ALLOCATE(psitot(ntot))
  ELSE
     ALLOCATE(psitot(1))
  END IF
  !
  ! ... gather all psis arrays on the first node, in psitot
  !
  CALL mp_barrier( intra_bgrp_comm )
  !
  CALL mp_gather( c(:,jw), psitot, ngpwpp, displs, root_bgrp, intra_bgrp_comm )
  !
  ! write the node-number-independent array
  !
  IF( me_bgrp == root_bgrp ) THEN
     DO i=1,ntot
        WRITE(22,*) psitot(i)
     END DO
  END IF
  !
  DEALLOCATE(psitot)

#else
  !
  DO i=1,ngw
     WRITE(22,*) c(i,jw)
  END DO
  !
#endif

  IF( ionode ) WRITE( stdout, * ) "State Written", jw
  !
  CALL stop_run( .TRUE. )
  !
  RETURN
  !
END SUBROUTINE write_psi
!
!----------------------------------------------------------------------------
SUBROUTINE jacobi_rotation( m, Omat, Umat, b1, b2, b3 )
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE wannier_base,           ONLY : nw, weight, nit, tolw, wfdt, maxwfdt, nsd
  USE cell_base,              ONLY : alat
  USE constants,              ONLY : tpi
  USE mp_global,              ONLY : me_bgrp
  USE printout_base,          ONLY : printout_base_open, printout_base_unit, &
                                     printout_base_close
  USE parallel_include
  !
  IMPLICIT NONE
  !    (m,m) is the size of the matrix Ospin.
  !    Ospin is input overlap matrix.
  !    Uspin is the output unitary transformation.
  !
  !    Jacobi rotations method is used to minimize the spread.
  !    (F. Gygi, J.-L. Fatterbert and E. Schwegler, Comput. Phys. Commun. 155, 1 (2003))
  !
  !    This subroutine has been written by Sylvie Stucki and Audrius Alkauskas
  !    in the Chair of Atomic Scale Simulation in Lausanne (Switzerland)
  !    under the direction of Prof. Alfredo Pasquarello.
  !
  INTEGER, intent(in) :: m
  COMPLEX(DP), DIMENSION(nw, m, m), intent(inout) :: Omat
  REAL(DP), DIMENSION(m, m), intent(inout) :: Umat
  LOGICAL :: stopCriteria
  INTEGER :: iterationNumber, lig, col, i, nbMat
  INTEGER, PARAMETER :: dimG=2
  REAL(DP), DIMENSION(2*nw, m, m):: OmatReal
  REAL(DP), DIMENSION(dimG, dimG):: matrixG
  REAL(DP), DIMENSION(dimG)      :: eigenVec
  REAL(DP) :: a1, a2 ! are the components aii-ajj and aij-aji of the matrixes used to build matrixG
  REAL(DP) :: r, c, s ! For a single rotation
  REAL(DP) :: bMinusa, outDiag ! To compute the eigenvector linked to the largest eigenvalue of matrixG
  REAL(DP) :: newMat_ll, newMat_cc, newMat_lc, presentSpread, saveSpread, mt(nw)
  REAL(DP), DIMENSION (m,2) :: newMat_cols
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  INTEGER :: me, iunit
  REAL(DP), DIMENSION(m, m) :: matriceTest
  !
  me = me_bgrp + 1
  nbMat=2*nw
  !
  WRITE(24,    *) 'Spreads before optimization'
  DO i=1, m
     !
     mt=1.D0-DBLE(Omat(:,i,i)*CONJG(Omat(:,i,i)))
     presentSpread = SUM(mt*weight)
     presentSpread = (alat/tpi)*DSQRT(presentSpread)
     WRITE(24,    *) 'Spread of the ', i, '-th wannier function is ' , presentSpread
     IF( presentSpread < 0.D0 ) &
          CALL errore( 'cp-wf', 'Something is wrong, WannierF spread negative', 1 )
     !
  ENDDO
  !
  Umat=0.D0
  DO i=1,m
     Umat(i,i)=1.D0
  END DO
  do i=1,m
     write (*, *) Umat(i, :)
  end do
  do i = 1, nw
     OmatReal((2*i-1), :, :) = real(Omat(i, :, :), DP)
     OmatReal(2*i, :, :) = aimag(Omat(i, :, :))
  end do
  !
  iterationNumber = 0 
  stopCriteria = .false.
  !
  ! Calculation of the spread
  presentSpread = 0.
  do i=1, nbMat
     do lig=1, m-1
        do col = lig+1, m
           presentSpread = presentSpread + OmatReal(i, lig, col)*OmatReal(i, lig, col)
        end do
     end do
  end do
  print *, "Initial spread : ", presentSpread
  saveSpread=presentSpread
  !
  ! ATTENTION! limite d'iteration = nit !!!!
  do while ((.not. stopCriteria) .and. (iterationNumber<nit))
     iterationNumber=iterationNumber + 1
     print *, "Tournus numero : ", iterationNumber
     do lig = 1, m-1
        do col =lig+1, m
           matrixG(1,1) = 0.d0
           matrixG(1,2) = 0.d0
           matrixG(2,1) = 0.d0
           matrixG(2,2) = 0.d0
           !
           !Building matrix G
           !         
           do i= 1, nbMat
              a1 = OmatReal(i, lig, lig)-OmatReal(i, col, col)
              a2 = OmatReal(i, lig, col)+OmatReal(i, col, lig)
              !print *, "a1 and a2 :", a1, a2
              matrixG(1,1) = matrixG(1,1) + a1*a1
              matrixG(1,2) = matrixG(1,2) + a1*a2
              matrixG(2,1) = matrixG(2,1) + a2*a1
              matrixG(2,2) = matrixG(2,2) + a2*a2
           end do
           !
           ! Computation of the eigenvector associated with the largest eigenvalue of matrixG
           bMinusa = matrixG(2,2)-matrixG(1,1)
           outDiag = matrixG(1,2)
           !
           if (abs(outDiag) .gt. 1d-10) then
              eigenVec(1) = 1.d0
              eigenVec(2) = (bMinusa + sqrt(bMinusa*bMinusa + &
                                    4.d0*outDiag*outDiag))/(2.d0*outDiag)
           else
              if (bMinusa .lt. 0) then
                 eigenVec(1)=1.d0
                 eigenVec(2)=0.d0
              else
                 eigenVec(1)=0.d0
                 eigenVec(2)=1.d0
              end if
           end if
           !  
           r = sqrt(eigenVec(1)*eigenVec(1) + eigenVec(2)*eigenVec(2))
           c = sqrt((eigenVec(1)+r)/(2.d0*r))
           s = eigenVec(2)/sqrt(2.d0*r*(eigenVec(1)+r))
           !
           !Update of the matrixes
           !
           do i= 1, nbMat
              ! Computation of the new components
              newMat_ll = c*c*OmatReal(i, lig, lig) +s*s*OmatReal(i, col, col) &
                   + 2*s*c*OmatReal(i, lig, col)
              newMat_cc = s*s*OmatReal(i, lig, lig) +c*c*OmatReal(i, col, col) &
                   - 2*s*c*OmatReal(i, lig, col)
              !
              newMat_lc = s*c* (-OmatReal(i, lig, lig) + OmatReal(i, col, col))&
                   + (c*c-s*s) * OmatReal(i, lig, col)
              newMat_cols(:, 1) = c*OmatReal(i, :, lig) + s*OmatReal(i, :, col)
              newMat_cols(:, 2) = -s*OmatReal(i, :, lig) + c*OmatReal(i, :, col)
              !
              ! copy of the components
              !  
              OmatReal(i, :, lig) = newMat_cols(:,1)
              OmatReal(i, :, col) = newMat_cols(:,2)
              OmatReal(i, lig, :) = newMat_cols(:,1)
              OmatReal(i, col, :) = newMat_cols(:,2)
              !
              OmatReal(i, lig, col) = newMat_lc
              OmatReal(i, col, lig) = newMat_lc
              OmatReal(i, lig, lig) = newMat_ll
              OmatReal(i, col, col) = newMat_cc
           end do
           !
           ! Update of the matrix of base transformation
           ! We reuse newMat_cols to keep the new values before the copy
           newMat_cols(:, 1) = c*Umat(:, lig) + s*Umat(:, col)
           newMat_cols(:, 2) = -s*Umat(:, lig) + c*Umat(:, col)
           !
           Umat(:, lig) = newMat_cols(:, 1)
           Umat(:, col) = newMat_cols(:, 2)
        end do
     end do
     !matriceTest = matmul(Umat, transpose(Umat))
     !write (*, *) matriceTest
     !
     ! Calculation of the new spread
     !
     presentSpread = 0.d0
     !
     do i=1, nbMat
        do lig=1, m-1
           do col = lig+1, m
              presentSpread = presentSpread + OmatReal(i, lig, col)*OmatReal(i, lig, col)
           end do
        end do
     end do
     !
     print *, iterationNumber, presentSpread
     !
     if (saveSpread-presentSpread < tolw) then
        print *, "Arret : ", iterationNumber, presentSpread
        stopCriteria = .true.
     end if
     saveSpread = presentSpread
  enddo
  !
  do i=1, nw
     !
     Omat(i, :, :) = OmatReal((2*i-1), :, :) + CI*OmatReal(2*i, :, :)
     !
  end do
  ! 
  WRITE(24,    *) 'Spreads after optimization'
  DO i=1, m
     !
     mt=1.D0-DBLE(Omat(:,i,i)*CONJG(Omat(:,i,i)))
     presentSpread = SUM(mt*weight)
     presentSpread = (alat/tpi)*DSQRT(presentSpread)
     WRITE(24,    *) 'Spread of the ', i, '-th wannier function is ' , presentSpread
     IF( presentSpread < 0.D0 ) &
          CALL errore( 'cp-wf', 'Something is wrong, WannierF spread negative', 1 )
     !
  ENDDO
  
  RETURN
END SUBROUTINE jacobi_rotation
!==============================================================================

    SUBROUTINE ddyn_u(nbsp, O, U, iss)

! input: the overlap matrix O
! ouput: the unitary transformation matrix U

       USE kinds,            ONLY : DP
       USE wannier_base,     ONLY : wf_friction, nsteps, tolw, adapt, wf_q, weight, nw, wfdt
       USE cell_base,        ONLY : alat
       USE constants,        ONLY : tpi, autoaf => BOHR_RADIUS_ANGS
       USE mp_global,        ONLY : nproc_image, me_image, intra_image_comm
       USE cp_main_variables, ONLY: descla
       USE cp_interfaces,     ONLY: distribute_lambda, collect_lambda
       USE printout_base,     ONLY : printout_base_open, printout_base_unit, printout_base_close
       USE cp_main_variables, ONLY : nfi, iprint_stdout
       USE time_step,         ONLY : tps
       USE input_parameters,  ONLY : tcpbo
       USE control_flags,     ONLY : conv_elec 
       USE parallel_include
       USE io_global,         ONLY : ionode, stdout

       IMPLICIT NONE
   
       INTEGER ,      INTENT(in)    :: nbsp
       REAL(DP),      INTENT(out)   :: U(nbsp,nbsp)
       COMPLEX(DP),   INTENT(inout) :: O(nw,nbsp,nbsp)
       INTEGER ,      INTENT(in)    :: iss
   
       INTEGER                      :: ista(0:nproc_image-1),iend(0:nproc_image-1)
       REAL(DP),    ALLOCATABLE,  DIMENSION(:,:) :: identy,Um,Up,U0,Ul,W,X2,X3,tmpr2,tmpi2, tmpr,tmpi
       COMPLEX(DP), ALLOCATABLE,  DIMENSION(:, :, :) :: Oc, Ocold, Ol
   
       INTEGER  ::  i, j, inw, nmax, iter, iunit,ini, nr, nc, ir, ic, ierr, nlam
       REAL(DP) :: t0, myt0, fric, t2(nw), mt(nw), oldt0,fric1,spread,sp, eps, wfdt2, fricp, fricm

       nlam = 1
       !IF( SIZE( descla ) < 2 ) THEN
       !   IF( descla(1)%active_node > 0 ) &
       !      nlam = descla(1)%nrcx
       !ELSE
       !   IF( ( descla(1)%active_node > 0 ) .OR. ( descla(2)%active_node > 0 ) ) &
       !      nlam = MAX( descla(1)%nrcx, descla(2)%nrcx )
       !END IF
       IF ( descla(iss)%active_node > 0) &
           nlam = descla(iss)%nrcx
   
       ALLOCATE(Oc(nbsp,nbsp, nw), Ocold(nbsp,nbsp,nw), Ol(nlam,nlam,nw))
       ALLOCATE(Up(nlam,nlam), U0(nlam,nlam), Um(nlam,nlam), Ul(nlam,nlam),  X2(nbsp,nbsp), X3(nbsp,nbsp))
       ALLOCATE(W(nlam,nlam), identy(nlam,nlam))
       ALLOCATE(tmpr(nlam,nlam), tmpi(nlam,nlam), tmpr2(nlam,nlam), tmpi2(nlam,nlam))

       eps=1.0E-13_DP
       nmax=50
       fric=wf_friction
       oldt0=0.D0
      
       nr = descla(iss)%nr
       nc = descla(iss)%nc
       ir = descla(iss)%ir
       ic = descla(iss)%ic

       do inw = 1, nw
          X2(:,:) =  REAL(O(inw, :, :))
          X3(:,:) = AIMAG(O(inw, :, :))
     
          call distribute_lambda(X2, tmpr, descla(iss))
          call distribute_lambda(X3, tmpi, descla(iss))
      
          Oc(:,:,inw) = DCMPLX(X2,X3)
          Ol(:,:,inw) = DCMPLX(tmpr,tmpi)
       enddo
      
       Ocold = Oc
      
       X2=0.D0
       DO i=1,nbsp
          X2(i,i)=1.D0
       END DO
       call distribute_lambda(X2, identy, descla(iss))
    
       Ul = identy
    
       call para_range(1, nbsp, nproc_image, ista, iend)

       DO inw = 1,nw
          t2(inw) = 2.D0 * weight(inw)
       ENDDO
    
       U0 = identy            
       Up = identy            

       DO ini=1, nsteps
    
          Um = identy 
          U0 = identy 
   
          fric1 = fric/(2*wfdt)  
          wfdt2 = wfdt*wfdt/(wf_q*(1+fric1))
          fricp = 2.d0/(1+fric1)
          fricm = (1-fric1)/(1+fric1)
    
          DO j = 1, nc
             DO i = 1, nr
                W(i,j) = 0.D0
                DO inw = 1, nw
                   W(i,j) = W(i,j) + &
                   t2(inw)*REAL((Ocold(i+ir-1,j+ic-1,inw)+Ocold(j+ic-1,i+ir-1,inw))*CONJG(Oc(j+ic-1,j+ic-1,inw)) )
                ENDDO
                Up(i,j)=fricp*U0(i,j)-fricm*Um(i,j)+wfdt2*W(i,j)   
             END DO
          END DO
    
          CALL ortho_u(Up,U0,nlam,identy,eps,nmax,nbsp, iss)
          CALL sqr_mm_cannon( 'N', 'N', nbsp, 1.0d0, Ul, nlam, Up, nlam, 0.0d0, tmpr, nlam, descla(iss))  
          Ul = tmpr

          Ocold = Oc
          DO inw = 1, nw
             tmpr(:,:)=REAL(Ol(:,:,inw))
             tmpi(:,:)=AIMAG(Ol(:,:,inw))
   
             CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, Ul, nlam, tmpr, nlam, 0.0d0, tmpr2, nlam, descla(iss)) 
             CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, Ul, nlam, tmpi, nlam, 0.0d0, tmpi2, nlam, descla(iss)) 
             CALL sqr_mm_cannon( 'N', 'N', nbsp, 1.0d0, tmpr2, nlam, Ul, nlam, 0.0d0, tmpr, nlam, descla(iss)) 
             CALL sqr_mm_cannon( 'N', 'N', nbsp, 1.0d0, tmpi2, nlam, Ul, nlam, 0.0d0, tmpi, nlam, descla(iss)) 
   
             call collect_lambda(X2, tmpr, descla(iss) )
             call collect_lambda(X3, tmpi, descla(iss))
   
             Oc(:,:,inw)=CMPLX(X2,X3)
          ENDDO
   
!======================================================================

       myt0=0.D0
       DO i = ista(me_image), iend(me_image)
!      DO i = 1, nbsp
          DO inw=1, nw
             myt0=myt0+DBLE(CONJG(Oc( i, i, inw))*Oc( i, i, inw))
          END DO
       END DO
#if defined(__MPI)
       CALL mpi_allreduce (myt0, t0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           intra_image_comm, ierr)
#else
       t0 = myt0
#endif
!       if(mod(ini,10) == 0)WRITE(stdout,*), 'spread at ', ini, ' = ', t0
       IF(ABS(t0-oldt0).LT.tolw) THEN
         !WRITE(stdout,*) "MLWF Generated at Step",ini,"  Convergence = ",ABS(t0-oldt0)
         WRITE(stdout,'(3X,"MLWF step",I5,"  Convergence = ",ES10.3,"  Generated")')ini,ABS(t0-oldt0) ! BS
         exit
       ELSEIF(ini.GE.nsteps) THEN
         !WRITE(stdout,*) "MLWF Not generated after",ini,"Steps."
          WRITE(stdout,'(3X,"MLWF step",I5,"  Convergence = ",ES10.3,"  Not Generated")')ini,ABS(t0-oldt0) ! BS
       END IF

       IF(oldt0 .GT. t0) fric=fric/2.D0

       oldt0=t0

     END DO !cycl for nsteps

     !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
     IF( ( MOD( nfi, iprint_stdout ) == 0 ) )  THEN
        IF ( ionode ) THEN
           iunit = printout_base_unit( "spr" )
           CALL printout_base_open( "spr" )
           IF (.NOT.tcpbo) THEN
              WRITE( iunit, '(I8,f16.8)' )nfi,tps
           ELSE   
              IF(conv_elec) WRITE( iunit, '(I8,f16.8)' )nfi,tps
           END IF    
        END IF    
     END IF    

     spread=0.0d0
     DO i=1, nbsp
        mt=1.D0-DBLE(Oc(i,i,:)*CONJG(Oc(i,i,:)))
        sp = (alat*autoaf/tpi)**2*SUM(mt*weight)
        !
        !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
        IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
           IF ( ionode ) THEN
              IF (.NOT.tcpbo) THEN
                 WRITE( iunit, '(f20.14)' ) sp
              ELSE
                IF(conv_elec) WRITE( iunit, '(f20.14)' ) sp
              END IF   
           END IF   
        END IF
        !
        IF ( sp < 0.D0 ) &
           CALL errore( 'cp-wf', 'Something wrong WF Spread negative', 1 )  
        spread=spread+sp
     END DO

     !BS .. to print WANNIER spreads at every iprint steps in .spr file ..    
     IF( ( MOD( nfi, iprint_stdout ) == 0 ) ) THEN
        IF ( ionode ) CALL printout_base_close( "spr" )
     END IF

     spread=spread/nbsp
    !IF(me_image.EQ.0) write(*,'(3X,"Average spread =",ES10.3)') spread ! BS

     call collect_lambda(U,Ul,descla(iss))

     do inw = 1, nw
        O(inw,:,:) = Oc(:,:,inw)
     enddo

     DEALLOCATE(tmpr,tmpi,tmpr2,tmpi2,identy,Oc,Ocold,Up,U0,Um,W,Ol, Ul, X2, X3)

     RETURN
     END SUBROUTINE ddyn_u

!-------------------------------------------------------------------------

      SUBROUTINE  ortho_u(up,u0,nlam,identy,eps,nmax,nbsp, iss)
!-----------------------------------------------------------------------
!     input = up (non-unitary), u0 (must be unitary)
!     output = up (unitary),
!     the method used is similar (same) to the ortho (les houches 1988)
!      x^t x + b x + x^t b^t + a = 1
!     where b = up^t u0   a = up^t up
!
      USE kinds,            ONLY : DP
      USE mp_global,        ONLY : me_image, intra_image_comm
      USE mp,               ONLY : mp_max
      USE cp_main_variables,       ONLY  : descla
 
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: nlam, nmax, nbsp, iss
      REAL(DP),INTENT(INOUT) :: up(nlam,nlam)
      REAL(DP),INTENT(IN)    :: u0(nlam,nlam), identy(nlam,nlam),eps
      REAL(DP) :: delta
      INTEGER  :: i,j,iter, nc, nr

      REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: xloc,tmp, tmp2, tmp2t, amat, bmat
      ALLOCATE( xloc(nlam,nlam))
      ALLOCATE( tmp(nlam,nlam),tmp2(nlam,nlam),tmp2t(nlam,nlam) )
      ALLOCATE( amat(nlam,nlam),bmat(nlam,nlam))

      nr = descla(iss)%nr
      nc = descla(iss)%nc

      CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, up, nlam, up, nlam, 0.0d0, amat, nlam, descla(iss))
      CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, up, nlam, u0, nlam, 0.0d0, bmat, nlam, descla(iss))

      amat = identy-amat
      bmat = identy-bmat
      xloc = 0.5d0*amat

      delta = 1.0E10_DP
      DO iter = 1,nmax

         CALL sqr_mm_cannon( 'N', 'N', nbsp, 1.0d0, bmat, nlam, xloc, nlam, 0.0d0, tmp2, nlam, descla(iss))
         CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, xloc, nlam, xloc, nlam, 0.0d0, tmp, nlam, descla(iss))

         CALL sqr_tr_cannon( nbsp, tmp2, nlam, tmp2t, nlam, descla(iss) )

         do j=1,nc
            do i=1,nr
               xloc(i,j)=0.5d0*(amat(i,j) + tmp2(i,j) + tmp2t(i,j) - tmp(i,j) )
            end do
         end do

         IF(iter .GE. 3) THEN

            tmp = up         ! upnew = up + u0*xloc
            CALL sqr_mm_cannon( 'N', 'N', nbsp, 1.0d0, u0, nlam, xloc, nlam, 1.0d0, tmp, nlam, descla(iss))

            tmp2 = identy
            CALL sqr_mm_cannon( 'T', 'N', nbsp, 1.0d0, tmp, nlam, tmp, nlam, -1.0d0, tmp2, nlam, descla(iss))
            delta = 0.d0
            do j=1,nc
            do i=1,nr
               delta=max(delta,abs(tmp2(i,j)) )
            end do
            end do
            CALL mp_max( delta, intra_image_comm )
            IF( delta .le. eps ) exit
         ENDIF

      ENDDO

      up = tmp
      DEALLOCATE(xloc, amat, bmat, tmp, tmp2, tmp2t)
      return
    
      END SUBROUTINE ortho_u

      SUBROUTINE para_range(n1, n2, nprocs, ista, iend)
      INTEGER n1,n2,nprocs,ista(0:nprocs-1),iend(0:nprocs-1)
      INTEGER iwork1,iwork2

      iwork1 = (n2 - n1 + 1) / nprocs
      iwork2 = MOD(n2 - n1 + 1, nprocs)
      do irank = 0,nprocs-1
         ista(irank) = irank * iwork1 + n1 + MIN(irank, iwork2)
         iend(irank) = ista(irank) + iwork1 - 1
         IF (iwork2 > irank) iend(irank) = iend(irank) + 1
      enddo

      return
      end
