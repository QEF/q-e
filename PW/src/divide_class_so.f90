!
! Copyright (C) 2006-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------------
SUBROUTINE divide_class_so(code_group,nrot,smat,d_spin,has_e,nclass, &
                           nelem,elem, which_irr)
!-----------------------------------------------------------------------------
!
! This subroutine receives as input a set of nrot 3x3 matrices smat, 
! and nrot complex 2x2 matrices d_spin, which  are assumed to be the 
! operations of the point group given by code_group. Only the operations
! that do not contain the 2\pi rotation (-E) are given in input.
! smat are in cartesian coordinates.
! This routine divides the double group in classes and find:
!
! nclass         the number of classes of the double group
! nelem(iclass)  for each class, the number of elements of the class
! elem(i,iclass) 1<i<nelem(iclass) for each class tells which matrices 
!                smat belong to that class
! has_e(i,iclass) =-1 if the operation is multiplied by -E, 1 otherwise
! which_irr(iclass) for each class gives the position of that class in the
!                character table associated with the group and provided
!                by the routine set_irr_rap_so. NB: changing the order of
!                the elements in the character table must reflect in 
!                a change to which_irr. Presently the character tables 
!                are those given by G.F. Koster, Space Group and their
!                representations.
!                Several equivalent names for the irreducible representation
!                are given. D, G, L, S are used for Delta, Gamma, Lambda 
!                and Sigma.
!

USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: & 
          code_group,   &  ! The code of the point group
          nrot,         &  ! The number of symmetry operation
          nclass,       &  ! The number of classes
          nelem(24),    &  ! The elements of each class 
          elem(12,24),  &  ! Which elements in the smat list for each class
          has_e(12,24), &  ! if -1 the element is multiplied by -E
          which_irr(24)    ! See above 

REAL(DP) :: smat(3,3,nrot), cmat(3,3), ax(3), bx(3), cx(3)
REAL(DP) :: smate(3,3,2*nrot)
COMPLEX(DP) :: d_spin(2,2,48), d_spine(2,2,96), c_spin(2,2)

INTEGER :: done(96), irot, jrot, krot, iclass, i, other, other1
INTEGER :: tipo_sym, set_e, ipol, axis, axis1, axis2, ts, nused, iaxis(4), &
           iax, ibx, icx, aclass, bclass, cclass,  &
           imax, imbx, imcx, amclass, bmclass, cmclass, ind2(3)  
REAL(DP), PARAMETER :: eps = 1.d-7
REAL(DP) :: angle_rot, angle_rot_s, ars, ax_save(3,3:5), angle_vectors
LOGICAL :: compare_mat_so, is_axis, isok, isok1, is_parallel
LOGICAL :: done_ax(6)
!
! Divide the group in classes.
!
DO irot=1,nrot
   smate(:,:,irot)=smat(:,:,irot)
   smate(:,:,irot+nrot)=smat(:,:,irot)
   d_spine(:,:,irot)=d_spin(:,:,irot)
   d_spine(:,:,irot+nrot)=-d_spin(:,:,irot)
END DO
!
!  If there are doubts that the input matrices are not a point group uncomment
!  the call to this routine.
!
!CALL check_tgroup(2*nrot,d_spine,smate)
!

nclass=0
nelem=0
done=0
DO irot=1,2*nrot
   IF (done(irot)==0) THEN
      nclass=nclass+1
      DO jrot=1,2*nrot
         CALL coniug_mat_so(smate(1,1,jrot),d_spine(1,1,jrot), &
                            smate(1,1,irot),d_spine(1,1,irot), &
                            cmat,c_spin)
         DO krot=1,2*nrot
            IF (compare_mat_so(cmat,c_spin,smate(1,1,krot),d_spine(1,1,krot)) &
                          .AND.done(krot)==0) THEN
               nelem(nclass)=nelem(nclass)+1
               IF (krot.le.nrot) THEN
                  elem(nelem(nclass),nclass)=krot
                  has_e(nelem(nclass),nclass)=1
               ELSE
                  elem(nelem(nclass),nclass)=krot-nrot
                  has_e(nelem(nclass),nclass)=-1
               ENDIF
               done(krot)=1 
            ENDIF 
         ENDDO
      ENDDO
   ENDIF
ENDDO
!
!  For each class we should now decide which_irr. This depends on the group
!  and on the tables of characters of the irreducible representation.
!
which_irr(1)=1
IF (code_group==1) THEN
!
!  C_1 
!
   IF (nclass /= 2) CALL errore('divide_class_so','Wrong classes for C_1',1)
   which_irr(2)=2
ELSEIF (code_group==2.OR.code_group==3.OR.code_group==4) THEN
!
!  C_i, C_s, C_2 
!
   IF (nclass /= 4) &
            CALL errore('divide_class_so','Wrong classes for C_i, C_s or C_2',1)

   DO iclass=2,nclass
      IF (tipo_sym(smat(1,1,elem(1,iclass)))==1) THEN
         which_irr(iclass)=2
      ELSE
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      END IF
   END DO
ELSEIF (code_group==5) THEN
!
!  C_3  
!
! The function angle_rot(smat) provides the rotation angle of the matrix smat
!
   IF (nclass /= 6) CALL errore('divide_class_so','Wrong classes for C_3',1)
   DO iclass=2,nclass
      IF (tipo_sym(smat(1,1,elem(1,iclass)))==1) THEN
         which_irr(iclass)=2
      ELSE
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         ENDIF
      ENDIF
   ENDDO
ELSEIF (code_group==6) THEN
!
!  C_4 
!
   IF (nclass /= 8) CALL errore('divide_class_so','Wrong classes for C_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-90.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSEIF (ABS(ars-270.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),7)
         ELSE
            CALL errore('divide_class_so','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==7) THEN
!
!  C_6 
!
   IF (nclass /= 12) CALL errore('divide_class_so','Wrong classes for C_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),9)
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         ELSE
            CALL errore('divide_class_so','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==8) THEN
!
!  D_2  
!
   IF (nclass /= 5) CALL errore('divide_class_so','Wrong classes for D_2',1)
!
!  first search -E
!
   nused=1
   DO iclass=2, nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE
         iaxis(nused)=iclass
         nused=nused+1
      ENDIF
   ENDDO 

   CALL versor(smat(1,1,elem(1,iaxis(1))),ax)
   CALL which_c2(ax,iax)
   CALL versor(smat(1,1,elem(1,iaxis(2))),bx)
   CALL which_c2(bx,ibx)
   CALL versor(smat(1,1,elem(1,iaxis(3))),cx)
   CALL which_c2(cx,icx)

   CALL is_d2(iax, ibx, icx, ind2)
 
   which_irr(iaxis(1))=ind2(1)+2
   which_irr(iaxis(2))=ind2(2)+2
   which_irr(iaxis(3))=ind2(3)+2

ELSEIF (code_group==9) THEN
!
!  D_3 
!
   IF (nclass /= 6) CALL errore('divide_class_so','Wrong classes for D_3',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSEIF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==10) THEN
!
!  D_4 
!
   IF (nclass /= 7) CALL errore('divide_class_so','Wrong classes for D_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         axis1=MOD(ipol,3)+1
         axis2=MOD(ipol+1,3)+1
         IF (axis==0) call errore('divide_class_so','unknown D_4 axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=5
         ELSEIF (is_axis(ax,axis1).or.is_axis(ax,axis2)) THEN
            which_irr(iclass)=6
         ELSE
            which_irr(iclass)=7
         END IF
      ELSEIF (ts.ne.3) THEN
         CALL errore('divide_class_so','wrong sym_type for D_4',1)
      END IF
   END DO
ELSEIF (code_group==11) THEN
!
!  D_6 
!
   IF (nclass /= 9) CALL errore('divide_class_so','Wrong classes for D_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps) ) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         ENDIF
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,3)) THEN
            which_irr(iclass)=7
         ELSE 
            CALL which_c2(ax, iax)
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=8
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=9
            ELSE
               CALL errore('divide_sym_so','D_6, C_2 axis not recognized',1)
            END IF
         END IF
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      END IF
   END DO
ELSEIF (code_group==12) THEN
!
!  C_2v 
!
   IF (nclass /= 5) CALL errore('divide_class_so','Wrong classes for C_2v',1)
   iax=0
   ibx=0
   icx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2( ax, iax) 
         which_irr(iclass)=3
      ELSEIF (ts==5) THEN
         IF (ibx==0) THEN
            CALL mirror_axis(smat(1,1,elem(1,iclass)), bx)
            CALL which_c2( bx, ibx) 
            bclass=iclass
         ELSE
            CALL mirror_axis(smat(1,1,elem(1,iclass)), bx)
            CALL which_c2( bx, icx) 
            cclass=iclass
         ENDIF
      ENDIF
   ENDDO
   CALL is_c2v(iax, ibx, icx, isok)
   IF (isok) THEN
      which_irr(bclass)=4
      which_irr(cclass)=5
   ELSE
      CALL is_c2v(iax, icx, ibx, isok1)
      IF (.NOT.isok1) CALL errore('divide_class_so','problem with C_2v',1)
      which_irr(bclass)=5
      which_irr(cclass)=4
   ENDIF

ELSEIF (code_group==13) THEN
!
!  C_3v 
!
   IF (nclass /= 6) CALL errore('divide_class_so','Wrong classes for C_3v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSEIF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==14) THEN
!
!  C_4v 
!
   IF (nclass /= 7) CALL errore('divide_class_so','Wrong classes for C_4v',1)

   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSEIF (ts==4) THEN
         which_irr(iclass)=5
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax, iax)
         IF (iax < 4) THEN 
            which_irr(iclass)=6
         ELSE
            which_irr(iclass)=7
         ENDIF
      ENDIF
   ENDDO

ELSEIF (code_group==15) THEN
!
!  C_6v 
!
   IF (nclass /= 9) CALL errore('divide_class_so','Wrong classes for C_6v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         ENDIF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=7
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax,iax)
         IF (iax==2 .OR. iax==12 .OR. iax==13) THEN
            which_irr(iclass)=9
         ELSEIF (iax==1 .OR. iax==10 .OR. iax==11) THEN
            which_irr(iclass)=8
         ELSE
            CALL errore('divide_class_so','C_6v mirror symmetry not recognized',1)
         ENDIF
      ENDIF
   ENDDO
ELSEIF (code_group==16) THEN
!
!  C_2h 
!
   IF (nclass /= 8) CALL errore('divide_class_so','Wrong classes for C_2h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSEIF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSEIF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==17) THEN
!
!  C_3h 
!
   IF (nclass /= 12) CALL errore('divide_class_so','Wrong classes for C_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         END IF
      ELSEIF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSEIF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),9)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==18) THEN
!
!  C_4h 
!
   IF (nclass /= 16) CALL errore('divide_class_so','Wrong classes for C_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         IF (angle_rot(smat(1,1,elem(1,iclass)))-90.d0<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),7)
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSEIF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),9)
      ELSEIF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),13)
      ELSEIF (ts==6) THEN
         IF (angle_rot_s(smat(1,1,elem(1,iclass)))-90.d0<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),15)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==19) THEN
!
!  C_6h 
!
   IF (nclass /= 24) CALL errore('divide_class_so','Wrong classes for C_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),9)
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSEIF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),13)
      ELSEIF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),19)
      ELSEIF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),21)
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),23)
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),15)
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),17)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==20) THEN
!
!  D_2h 
!
!  mirror_axis gives the normal to the mirror plane
!
   IF (nclass /= 10) CALL errore('divide_class_so','Wrong classes for D_2h',1)
!
!  First check if the axis are parallel to x, y or z
!
   iax=0
   ibx=0
   icx=0
   imax=0
   imbx=0
   imcx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (iax==0) THEN
            CALL which_c2(ax, iax)
            aclass=iclass
         ELSEIF (ibx==0) THEN
            CALL which_c2(ax, ibx)
            bclass=iclass
         ELSEIF (icx==0) THEN
            CALL which_c2(ax, icx)
            cclass=iclass
         ELSE
            CALL errore('divide_class_so','D_2h too many C_2 axis',1)
         ENDIF 
      ELSEIF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),6)
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
         IF (imax==0) THEN
            CALL which_c2(ax, imax)
            amclass=iclass
         ELSEIF (imbx==0) THEN
            CALL which_c2(ax, imbx)
            bmclass=iclass
         ELSEIF (imcx==0) THEN
            CALL which_c2(ax, imcx)
            cmclass=iclass
         ELSE
            CALL errore('divide_class_so','D_2h too many mirrors',1)
         ENDIF 
      ELSE
         CALL errore('divide_class_so','D_2h operation not recognized',1)
      ENDIF
   ENDDO

   CALL is_d2( iax, ibx, icx, ind2)

   which_irr(aclass)=ind2(1)+2
   which_irr(bclass)=ind2(2)+2
   which_irr(cclass)=ind2(3)+2
 
   IF (imax==iax) which_irr(amclass) = which_irr(aclass) + 5
   IF (imax==ibx) which_irr(amclass) = which_irr(bclass) + 5
   IF (imax==icx) which_irr(amclass) = which_irr(cclass) + 5
   IF (imbx==iax) which_irr(bmclass) = which_irr(aclass) + 5
   IF (imbx==ibx) which_irr(bmclass) = which_irr(bclass) + 5
   IF (imbx==icx) which_irr(bmclass) = which_irr(cclass) + 5
   IF (imcx==iax) which_irr(cmclass) = which_irr(aclass) + 5
   IF (imcx==ibx) which_irr(cmclass) = which_irr(bclass) + 5
   IF (imcx==icx) which_irr(cmclass) = which_irr(cclass) + 5

ELSEIF (code_group==21) THEN
!
!  D_3h 
!
   IF (nclass /= 9) CALL errore('divide_class_so','Wrong classes for D_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSE IF (ts==4) THEN
         which_irr(iclass)=5
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)>1) THEN
            which_irr(iclass)=9
         ELSE 
            which_irr(iclass)=6
         END IF
      ELSE IF (ts==6) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      END IF
   END DO
ELSEIF (code_group==22) THEN
!
!  D_4h 
!
!  First search the order 4 axis
!
   IF (nclass /= 14) CALL errore('divide_class_so','Wrong classes for D_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         IF (axis==0) call errore('divide_class_so','unknown D_4h axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=0
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=5
         ELSE
           DO ipol=1,3
              IF (is_axis(ax,ipol)) which_irr(iclass)=6
           ENDDO
           IF (which_irr(iclass)==0) which_irr(iclass)=7
         END IF
      ELSEIF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),8)
      ELSEIF (ts==5) THEN
         which_irr(iclass)=0
         CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=12
         ELSE 
            DO ipol=1,3
               IF (is_axis(ax,ipol)) which_irr(iclass)=13
            ENDDO
            IF (which_irr(iclass)==0) which_irr(iclass)=14
         END IF
      ELSEIF (ts==6) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),10)
      END IF
   END DO
ELSEIF (code_group==23) THEN
!
!  D_6h 
!
   IF (nclass /= 18) CALL errore('divide_class_so','Wrong classes for D_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         END IF
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==2) THEN
            which_irr(iclass)=7
         ELSE
            CALL versor(smat(1,1,elem(1,iclass)),ax)
            CALL which_c2(ax,iax) 
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=8
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=9
            ELSE
               CALL errore('divide_class_so','Problem with C_2 of D_6h',1)
            ENDIF
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),10)
      ELSE IF (ts==5) THEN
          IF (nelem(iclass)==2) THEN
             which_irr(iclass)=16
          ELSE 
             CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
             CALL which_c2(ax, iax)
             IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
!
!   sigma_d
!
                which_irr(iclass)=17
             ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
!
!   sigma_v
!
                which_irr(iclass)=18
             ELSE
               CALL errore('divide_class_so','Problem with mirror of D_6h',1)
             ENDIF
          END IF
      ELSE IF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),14)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),12)
         END IF
      END IF
   END DO
ELSEIF (code_group==24) THEN
!
!  D_2d 
!
   IF (nclass /= 7) CALL errore('divide_class_so','Wrong classes for D_2d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
          which_irr(iclass)=2
      ELSE IF (ts==6) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==2) THEN
            which_irr(iclass)=5
         ELSE
            which_irr(iclass)=6
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=7
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==25) THEN
!
!  D_3d 
!
   IF (nclass /= 12) CALL errore('divide_class_so','Wrong classes for D_3d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSEIF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSE IF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSE IF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSE IF (ts==6) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),9)
      ELSE IF (ts==5) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),11)
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==26) THEN
!
!  S_4 
!
   IF (nclass /= 8) CALL errore('divide_class_so','Wrong classes for S_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),5)
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-90.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),7)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSE IF (code_group==27) THEN
!
!  S_6 
!
   IF (nclass /= 12) CALL errore('divide_class_so','Wrong classes for S_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),5)
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),7)
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),9)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==28) THEN
!
!  T 
!
   IF (nclass /= 7) CALL errore('divide_class_so','Wrong classes for T',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),4)
         ELSE IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-240.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),6)
         ENDIF
      ELSE
         CALL errore('divide_class_so','wrong sym type',1)
      END IF
   END DO
ELSE IF (code_group==29) THEN
!
!  T_h 
!
   IF (nclass /= 14) CALL errore('divide_class_so','Wrong classes for T_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),4)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),6)
         END IF
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),8)
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),13)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=10
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==30) THEN
!
!  T_d 
!
   IF (nclass /= 8) CALL errore('divide_class_so','Wrong classes for T_d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==3) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),3)
      ELSE IF (ts==4) THEN
         which_irr(iclass)=5
      ELSE IF (ts==6) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),6)
      ELSE IF (ts==5) THEN
         which_irr(iclass)=8
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==31) THEN
!
!  O 
!
   IF (nclass /= 8) CALL errore('divide_class_so','Wrong classes for O',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==1) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==3) THEN
            which_irr(iclass)=5
         ELSE
            which_irr(iclass)=8
         ENDIF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),6)
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==32) THEN
!
!  O_h
!
   IF (nclass /= 16) CALL errore('divide_class_so','Wrong classes for O_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         IF (nelem(iclass)==6) THEN
            which_irr(iclass)=5
         ELSE
            which_irr(iclass)=8
         END IF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),3)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),6)
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=set_e(has_e(1,iclass),9)
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)==12) THEN
            which_irr(iclass)=16
         ELSE
            which_irr(iclass)=13
         END IF
      ELSE IF (ts==6) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=set_e(has_e(1,iclass),11)
         ELSE
            which_irr(iclass)=set_e(has_e(1,iclass),14)
         END IF
      ELSE IF (ts==1) THEN
           which_irr(iclass)=2
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   ENDDO
ELSE
 CALL errore('divide_class_so','code_group not correct',1)
ENDIF

RETURN
END SUBROUTINE divide_class_so

!-----------------------------------------------------------------------------
SUBROUTINE coniug_mat_so(a,a_spin,b,b_spin,c,c_spin)
!-----------------------------------------------------------------------------
USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP) :: a(3,3), b(3,3), c(3,3)
COMPLEX(DP) :: a_spin(2,2), b_spin(2,2), c_spin(2,2)

c=MATMUL(a,MATMUL(b,TRANSPOSE(a)))
c_spin=MATMUL(a_spin,MATMUL(b_spin,TRANSPOSE(CONJG(a_spin))))

RETURN
END SUBROUTINE coniug_mat_so
!
!-----------------------------------------------------------------------------
FUNCTION compare_mat_so(a,a_spin,b,b_spin)
!-----------------------------------------------------------------------------
!
!  This function compare two 3x3 matrices and two 2x2 matrices 
!  and returns .true. if they coincide.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: a(3,3), b(3,3), csum
REAL(DP), PARAMETER :: eps=1.d-7
COMPLEX(DP) :: a_spin(2,2), b_spin(2,2)
LOGICAL :: compare_mat_so
INTEGER :: i, j

csum=0.d0
DO i=1,2
   DO j=1,2
      csum=csum+ABS(a_spin(i,j)-b_spin(i,j))
   END DO
END DO
compare_mat_so=((ABS(MAXVAL(a-b))<eps).AND.             &
                (ABS(MINVAL(a-b))<eps).AND. (csum < eps ))

RETURN
END FUNCTION compare_mat_so

!-----------------------------------------------------------------------------
SUBROUTINE set_irr_rap_so(code_group,nclass_ref,nrap_ref,char_mat,& 
                          name_rap,name_class,name_class1)
!-----------------------------------------------------------------------------
!
!  This subroutine collects the character tables of the 32 crystallographic
!  double point groups. Various names have been used in the litterature 
!  to identify D, G, L, S are used for Delta, Gamma, Lambda and Sigma.
!   
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: nclass_ref, &    ! Output: number of classes
           nrap_ref,   &    ! Output: number of irreducible representation
           code_group       ! Input: code of the group 

CHARACTER(LEN=15) :: name_rap(24)   ! Output: name of the representations
CHARACTER(LEN=5) :: name_class(24), & ! Output: name of the classes
                    name_class1(24) ! Output: name of the classes

COMPLEX(DP) :: char_mat(12,24) ! Output: character matrix

REAL(DP) :: sqr3d2, sqrt2, sqrt3, dsq2

sqrt2 =SQRT(2.d0)
sqrt3 =SQRT(3.d0)
sqr3d2=sqrt3*0.5d0
dsq2  =sqrt2*0.5d0

char_mat=(1.d0,0.d0)
char_mat(:,2)=(-1.d0,0.d0)

name_class1="     "
name_class(1)="E   "
name_class(2)="-E  "


IF (code_group==1) THEN
!
! C_1
!
   nclass_ref=2
   nrap_ref=1

   name_rap(1)="G_2  "

ELSEIF (code_group==2) THEN
!
! C_i
!
   nclass_ref=4

   name_class(3)="i    "
   name_class(4)="-i   "

   nrap_ref=2

   name_rap(1)="G_2+"
   char_mat(1,4)=(-1.d0,0.d0)

   name_rap(2)="G_2-"
   char_mat(2,3)=(-1.d0,0.d0)

ELSEIF (code_group==3) THEN
!
! C_s
!
   nclass_ref=4
   name_class(3)="s    "
   name_class(4)="-s   "

   nrap_ref=2

   name_rap(1)="G_3   "
   char_mat(1,3)=(0.d0, 1.d0)
   char_mat(1,4)=(0.d0,-1.d0)

   name_rap(2)="G_4   "
   char_mat(2,3)=(0.d0,-1.d0)
   char_mat(2,4)=(0.d0, 1.d0)

ELSEIF (code_group==4) THEN
!
! C_2
!
   nclass_ref=4
   name_class(3)="C2   "
   name_class(4)="-C2  "

   nrap_ref=2

   name_rap(1)="G_3   "
   char_mat(1,3)=(0.d0, 1.d0)
   char_mat(1,4)=(0.d0,-1.d0)

   name_rap(2)="G_4   "
   char_mat(2,3)=(0.d0,-1.d0)
   char_mat(2,4)=(0.d0, 1.d0)


ELSEIF (code_group==5) THEN
!
! C_3    NB: The signs of the characters of the classes C3^2 -C3^2  
!            are changed with respect to Koster, Space groups and 
!            their representations. They match the table in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups.
!
   nclass_ref=6
   name_class(3)="C3   "
   name_class(4)="-C3  "
   name_class(5)="C3^2 "
   name_class(6)="-C3^2"

   nrap_ref=3

   name_rap(1)="G_4  "
   char_mat(1,3)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(1,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(1,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(1,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)

   name_rap(2)="G_5  "
   char_mat(2,3)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(2,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   
   name_rap(3)="G_6  "
   char_mat(3,3)=(-1.0d0,0.d0)
!   char_mat(3,6)=(-1.0d0,0.d0)
   char_mat(3,5)=(-1.0d0,0.d0)


ELSEIF (code_group==6) THEN
!
! C_4   NB: The signs of the characters of the class C4^3 -C4^3 
!            are changed with respect to Koster, Space groups and
!            their representation. 
!
   nclass_ref=8
   name_class(3)="C4   "
   name_class(4)="-C4  "
   name_class(5)="C2   "
   name_class(6)="-C2  "
   name_class(7)="C4^3 "
   name_class(8)="-C4^3"

   nrap_ref=4

   name_rap(1)="G_5  "
   char_mat(1,3)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(1,4)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(1,5)=( 0.d0,1.d0)
   char_mat(1,6)=( 0.d0,-1.d0)
!   char_mat(1,7)=CMPLX(-dsq2, dsq2,kind=DP)
!   char_mat(1,8)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,7)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,8)=CMPLX(-dsq2, dsq2,kind=DP)

   name_rap(2)="G_6  "
   char_mat(2,3)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(2,4)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(2,5)=( 0.d0,-1.d0)
   char_mat(2,6)=( 0.d0, 1.d0)
!   char_mat(2,7)=CMPLX(-dsq2,-dsq2,kind=DP)
!   char_mat(2,8)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,7)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,8)=CMPLX(-dsq2,-dsq2,kind=DP)

   name_rap(3)="G_7  "
   char_mat(3,3)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(3,4)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(3,5)=( 0.d0,1.d0)
   char_mat(3,6)=( 0.d0,-1.d0)
!   char_mat(3,7)=CMPLX( dsq2,-dsq2,kind=DP)
!   char_mat(3,8)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,7)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,8)=CMPLX( dsq2,-dsq2,kind=DP)

   name_rap(4)="G_8  "
   char_mat(4,3)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(4,4)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(4,5)=( 0.d0,-1.d0)
   char_mat(4,6)=( 0.d0, 1.d0)
!   char_mat(4,7)=CMPLX( dsq2, dsq2,kind=DP)
!   char_mat(4,8)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,7)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,8)=CMPLX( dsq2, dsq2,kind=DP)


ELSEIF (code_group==7) THEN
!
! C_6   NB: The signs of characters of the C3^2 -C3^2 C6^5 -C6^5 classes
!           are changed with respect to Koster, Space groups and their
!           representation. They match the table in Koster, Dimmock, 
!           Wheeler, Statz, Properties of the 32 point groups.
!
   nclass_ref=12
   name_class(3)="C6  "
   name_class(4)="-C6 "
   name_class(5)="C3 "
   name_class(6)="-C3 "
   name_class(7)="C2 "
   name_class(8)="-C2 "
   name_class(9)="C3^2"
   name_class(10)="-C3^2"
   name_class(11)="C6^5"
   name_class(12)="-C6^5"

   nrap_ref=6

   name_rap(1)="G_7  "
   char_mat(1,2)=(-1.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(1,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(1,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(1,7)=( 0.d0, 1.d0)
   char_mat(1,8)=( 0.d0,-1.d0)
!   char_mat(1,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(1,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(1,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(1,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(1,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)

   name_rap(2)="G_8  "
   char_mat(2,2)=(-1.d0, 0.d0)
   char_mat(2,3)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(2,4)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,7)=( 0.d0,-1.d0)
   char_mat(2,8)=( 0.d0, 1.d0)
!   char_mat(2,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(2,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)

   name_rap(3)="G_9   "
   char_mat(3,2)=(-1.d0, 0.d0)
   char_mat(3,3)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(3,4)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(3,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,7)=( 0.d0,-1.d0)
   char_mat(3,8)=( 0.d0, 1.d0)
!   char_mat(5,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(5,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(5,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(5,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)

   name_rap(4)="G_10  "
   char_mat(4,2)=(-1.d0, 0.d0)
   char_mat(4,3)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(4,4)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(4,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(4,7)=( 0.d0, 1.d0)
   char_mat(4,8)=( 0.d0,-1.d0)
!   char_mat(4,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(4,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(4,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(4,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)

   name_rap(5)="G_11 "
   char_mat(5,2)=(-1.d0, 0.d0)
   char_mat(5,3)=( 0.d0, 1.d0)
   char_mat(5,4)=( 0.d0,-1.d0)
   char_mat(5,5)=(-1.d0, 0.d0)
   char_mat(5,7)=( 0.d0,-1.d0)
   char_mat(5,8)=( 0.d0, 1.d0)
!   char_mat(5,10)=(-1.d0, 0.d0)
!   char_mat(5,11)=( 0.d0, 1.d0)
!   char_mat(5,12)=( 0.d0,-1.d0)
   char_mat(5,9)=(-1.d0, 0.d0)
   char_mat(5,11)=( 0.d0,-1.d0)
   char_mat(5,12)=( 0.d0, 1.d0)

   name_rap(6)="G_12  "
   char_mat(6,2)=(-1.d0, 0.d0)
   char_mat(6,3)=( 0.d0,-1.d0)
   char_mat(6,4)=( 0.d0, 1.d0)
   char_mat(6,5)=(-1.d0, 0.d0)
   char_mat(6,7)=( 0.d0, 1.d0)
   char_mat(6,8)=( 0.d0,-1.d0)
!   char_mat(6,10)=(-1.d0, 0.d0)
!   char_mat(6,11)=( 0.d0,-1.d0)
!   char_mat(6,12)=( 0.d0, 1.d0)
   char_mat(6,9)=(-1.d0, 0.d0)
   char_mat(6,11)=( 0.d0, 1.d0)
   char_mat(6,12)=( 0.d0,-1.d0)

ELSEIF (code_group==8) THEN
!
! D_2
!
   nclass_ref=5
   name_class(3)="C2 "
   name_class1(3)="-C2 "
   name_class(4)="C2' "
   name_class1(4)="-C2' "
   name_class(5)="C2''"
   name_class1(5)="-C2''"

   nrap_ref=1

   name_rap(1)="G_5  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=( 0.d0, 0.d0)
   char_mat(1,4)=( 0.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)

ELSEIF (code_group==9) THEN
!
! D_3
!
   nclass_ref=6
   name_class(3)="2C3  "
   name_class(4)="-2C3  "
   name_class(5)=" 3C2'"
   name_class(6)="-3C2'"

   nrap_ref=3

   name_rap(1)="G_4  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)

   name_rap(2)="G_5  "
   char_mat(2,2)=(-1.d0, 0.d0)
   char_mat(2,3)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 1.d0)
   char_mat(2,6)=( 0.d0,-1.d0)

   name_rap(3)="G_6  "
   char_mat(3,2)=(-1.d0, 0.d0)
   char_mat(3,3)=(-1.d0, 0.d0)
   char_mat(3,5)=( 0.d0,-1.d0)
   char_mat(3,6)=( 0.d0, 1.d0)


ELSEIF (code_group==10) THEN
!
! D_4
!
   nclass_ref=7
   name_class(3)="2C4  "
   name_class(4)="-2C4  "
   name_class(5)=" C2 "
   name_class1(5)="-C2 "
   name_class(6)="2C2'"
   name_class1(6)="-2C2'"
   name_class(7)="2C2'' "
   name_class1(7)="-2C2''"

   nrap_ref=2

   name_rap(1)="G_6  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)

   name_rap(2)="G_7  "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=( 0.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)


ELSEIF (code_group==12) THEN
!
! C_2v
!
   nclass_ref=5
   name_class(3)=" C2"
   name_class1(3)="-C2"
   name_class(4)=" s_v"
   name_class1(4)="-s_v"
   name_class(5)=" s_v'"
   name_class1(5)="-s_v'"

   nrap_ref=1

   name_rap(1)="G_5  D_5"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=( 0.d0, 0.d0)
   char_mat(1,4)=( 0.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)

ELSEIF (code_group==13) THEN
!
! C_3v
!
   nclass_ref=6
   name_class(3)="2C3  "
   name_class(4)="-2C3 "
   name_class(5)="3s_v"
   name_class(6)="-3s_v"

   nrap_ref=3

   name_rap(1)="G_4   L_6 "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)

   name_rap(2)="G_5   L_4 "
   char_mat(2,3)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 1.d0)
   char_mat(2,6)=( 0.d0,-1.d0)

   name_rap(3)="G_6   L_5 "
   char_mat(3,3)=(-1.d0, 0.d0)
   char_mat(3,5)=( 0.d0,-1.d0)
   char_mat(3,6)=( 0.d0, 1.d0)

ELSEIF (code_group==14) THEN
!
! C_4v
!
   nclass_ref=7
   name_class(3)="2C4  "
   name_class(4)="-2C4 "
   name_class(5)=" C2"
   name_class1(5)="-C2"
   name_class(6)=" 2s_v"
   name_class1(6)="-2s_v"
   name_class(7)=" 2s_d"
   name_class1(7)="-2s_d"

   nrap_ref=2

   name_rap(1)="G_6   D_6"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)

   name_rap(2)="G_7   D_7"
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=( 0.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)

ELSEIF (code_group==15.OR.code_group==11) THEN
!
! C_6v, D_6
!
   nclass_ref=9
   name_class(3)=" 2C6"
   name_class(4)="-2C6"
   name_class(5)=" 2C3 "
   name_class(6)="-2C3"
   name_class(7)=" C2 "
   name_class1(7)="-C2 "
   IF (code_group==15) THEN
      name_class(8)=" 3s_v"
      name_class1(8)="-3s_v"
      name_class(9)=" 3s_d"
      name_class1(9)="-3s_d"
   ELSE
      name_class(8)=" 3C2'"
      name_class1(8)="-3C2'"
      name_class(9)="3C2''"
      name_class1(9)="-3C2''"
   END IF

   nrap_ref=3

   name_rap(1)="G_7  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(1,6)=(-1.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)
   char_mat(1,8)=( 0.d0, 0.d0)
   char_mat(1,9)=( 0.d0, 0.d0)

   name_rap(2)="G_8  "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(2,6)=(-1.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)
   char_mat(2,8)=( 0.d0, 0.d0)
   char_mat(2,9)=( 0.d0, 0.d0)

   name_rap(3)="G_9  "
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=( 0.d0, 0.d0)
   char_mat(3,4)=( 0.d0, 0.d0)
   char_mat(3,5)=(-2.d0, 0.d0)
   char_mat(3,6)=( 2.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=( 0.d0, 0.d0)
   char_mat(3,9)=( 0.d0, 0.d0)

ELSEIF (code_group==16) THEN
!
! C_2h
!
   nclass_ref=8
   name_class(3)="C2  "
   name_class(4)="-C2 "
   name_class(5)="i   "
   name_class(6)="-i  "
   name_class(7)="s_h "
   name_class(8)="-s_h"

   nrap_ref=4

   name_rap(1)="G_3+"
   char_mat(1,3)=( 0.d0, 1.d0)
   char_mat(1,4)=( 0.d0,-1.d0)
   char_mat(1,6)=(-1.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 1.d0)
   char_mat(1,8)=( 0.d0,-1.d0)

   name_rap(2)="G_4+"
   char_mat(2,3)=( 0.d0,-1.d0)
   char_mat(2,4)=( 0.d0, 1.d0)
   char_mat(2,6)=(-1.d0, 0.d0)
   char_mat(2,7)=( 0.d0,-1.d0)
   char_mat(2,8)=( 0.d0, 1.d0)

   name_rap(3)="G_3-"
   char_mat(3,3)=( 0.d0, 1.d0)
   char_mat(3,4)=( 0.d0,-1.d0)
   char_mat(3,5)=(-1.d0, 0.d0)
   char_mat(3,7)=( 0.d0,-1.d0)
   char_mat(3,8)=( 0.d0, 1.d0)

   name_rap(4)="G_4- "
   char_mat(4,3)=( 0.d0,-1.d0)
   char_mat(4,4)=( 0.d0, 1.d0)
   char_mat(4,5)=(-1.d0, 0.d0)
   char_mat(4,7)=( 0.d0, 1.d0)
   char_mat(4,8)=( 0.d0,-1.d0)

ELSEIF (code_group==17) THEN
!
!     Changed several signs. Now the character table is equal to the table
!     of C_6 and the signs match the table in Koster, Dimmock, 
!     Wheeler, Statz, Properties of the 32 point groups.
!
! C_3h
!
   nclass_ref=12
   name_class(3)="C3  "
   name_class(4)="-C3 "
   name_class(5)="C3^2"
   name_class(6)="-C3^2"
   name_class(7)=" s_h "
   name_class(8)="-s_h "
   name_class(9)=" S3 "
   name_class(10)="-S3 "
   name_class(11)=" S3^5"
   name_class(12)="-S3^5"

   nrap_ref=6

   name_rap(1)="G_7  "
   char_mat(1,2)=(-1.d0, 0.d0)
   char_mat(1,3)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(1,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(1,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(1,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(1,7)=( 0.d0, 1.d0)
   char_mat(1,8)=( 0.d0,-1.d0)
!   char_mat(1,9)= CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(1,10)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,9)= CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,10)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(1,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(1,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)

   name_rap(2)="G_8  "
   char_mat(2,2)=(-1.d0, 0.d0)
   char_mat(2,3)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(2,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,7)=( 0.d0,-1.d0)
   char_mat(2,8)=( 0.d0, 1.d0)
!   char_mat(2,9)= CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(2,10)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,9)= CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,10)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(2,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(2,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)

   name_rap(3)="G_9  "
   char_mat(3,2)=(-1.d0, 0.d0)
   char_mat(3,3)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,7)=( 0.d0,-1.d0)
   char_mat(3,8)=( 0.d0, 1.d0)
!   char_mat(3,9)= CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(3,10)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,9)= CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,10)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(3,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(3,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)

   name_rap(4)="G_10  "
   char_mat(4,2)=(-1.d0, 0.d0)
   char_mat(4,3)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(4,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(4,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,7)=( 0.d0, 1.d0)
   char_mat(4,8)=( 0.d0,-1.d0)
!   char_mat(4,9)= CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(4,10)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,9)= CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,10)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(4,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(4,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)

   name_rap(5)="G_11 "
   char_mat(5,2)=(-1.d0, 0.d0)
   char_mat(5,3)=(-1.d0, 0.d0)
!   char_mat(5,5)=( 1.d0, 0.d0)
!   char_mat(5,6)=(-1.d0, 0.d0)
!   char_mat(5,7)=( 0.d0, 1.d0)
!   char_mat(5,8)=( 0.d0,-1.d0)
   char_mat(5,5)=(-1.d0, 0.d0)
   char_mat(5,6)=( 1.d0, 0.d0)
   char_mat(5,7)=( 0.d0,-1.d0)
   char_mat(5,8)=( 0.d0, 1.d0)
   char_mat(5,9)=( 0.d0,-1.d0)
   char_mat(5,10)=( 0.d0, 1.d0)
!   char_mat(5,11)=( 0.d0,-1.d0)
!   char_mat(5,12)=( 0.d0, 1.d0)
   char_mat(5,11)=( 0.d0, 1.d0)
   char_mat(5,12)=( 0.d0,-1.d0)

   name_rap(6)="G_12 "
   char_mat(6,2)=(-1.d0, 0.d0)
   char_mat(6,3)=(-1.d0, 0.d0)
!   char_mat(6,5)=( 1.d0, 0.d0)
!   char_mat(6,6)=(-1.d0, 0.d0)
!   char_mat(6,7)=( 0.d0,-1.d0)
!   char_mat(6,8)=( 0.d0, 1.d0)
   char_mat(6,5)=(-1.d0, 0.d0)
   char_mat(6,6)=( 1.d0, 0.d0)
   char_mat(6,7)=( 0.d0, 1.d0)
   char_mat(6,8)=( 0.d0,-1.d0)
   char_mat(6,9)=( 0.d0, 1.d0)
   char_mat(6,10)=( 0.d0,-1.d0)
!   char_mat(6,11)=( 0.d0, 1.d0)
!   char_mat(6,12)=( 0.d0,-1.d0)
   char_mat(6,11)=( 0.d0,-1.d0)
   char_mat(6,12)=( 0.d0, 1.d0)

ELSEIF (code_group==18) THEN
!
! C_4h   NB: The signs of the characters of the classes C4^3 -C4^2 S4 -S4 
!            are changed with respect to Koster, Space groups and 
!            their representations. They match the table in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups.
!
   nclass_ref=16
   name_class(3)="C4 "
   name_class(4)="-C4 "
   name_class(5)="C2  "
   name_class(6)="-C2 "
   name_class(7)="C4^3"
   name_class(8)="-C4^3"
   name_class(9)="i "
   name_class(10)="-i "
   name_class(11)="S4^3"
   name_class(12)="-S4^3"
   name_class(13)="s_h "
   name_class(14)="-s_h"
   name_class(15)="S4 "
   name_class(16)="-S4"

   nrap_ref=8

   name_rap(1)="G_5+ "
   char_mat(1,3)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(1,4)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(1,5)=( 0.d0, 1.d0)
   char_mat(1,6)=( 0.d0,-1.d0)
   char_mat(1,7)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,8)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(1,10)=(-1.d0, 0.d0)
   char_mat(1,11)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(1,12)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(1,13)=( 0.d0, 1.d0)
   char_mat(1,14)=( 0.d0,-1.d0)
   char_mat(1,15)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,16)=CMPLX(-dsq2, dsq2,kind=DP)

   name_rap(2)="G_6+ "
   char_mat(2,3)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(2,4)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(2,5)=( 0.d0,-1.d0)
   char_mat(2,6)=( 0.d0, 1.d0)
   char_mat(2,7)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,8)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(2,10)=(-1.d0, 0.d0)
   char_mat(2,11)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(2,12)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(2,13)=( 0.d0,-1.d0)
   char_mat(2,14)=( 0.d0, 1.d0)
   char_mat(2,15)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,16)=CMPLX(-dsq2,-dsq2,kind=DP)

   name_rap(3)="G_7+ "
   char_mat(3,3)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(3,4)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(3,5)=( 0.d0, 1.d0)
   char_mat(3,6)=( 0.d0,-1.d0)
   char_mat(3,7)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,8)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(3,10)=(-1.d0, 0.d0)
   char_mat(3,11)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(3,12)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(3,13)=( 0.d0, 1.d0)
   char_mat(3,14)=( 0.d0,-1.d0)
   char_mat(3,15)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,16)=CMPLX( dsq2,-dsq2,kind=DP)

   name_rap(4)="G_8+ "
   char_mat(4,3)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(4,4)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(4,5)=( 0.d0,-1.d0)
   char_mat(4,6)=( 0.d0, 1.d0)
   char_mat(4,7)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,8)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(4,10)=(-1.d0, 0.d0)
   char_mat(4,11)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(4,12)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(4,13)=( 0.d0,-1.d0)
   char_mat(4,14)=( 0.d0, 1.d0)
   char_mat(4,15)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,16)=CMPLX( dsq2, dsq2,kind=DP)

   name_rap(5)="G_5- "
   char_mat(5,3)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(5,4)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(5,5)=( 0.d0, 1.d0)
   char_mat(5,6)=( 0.d0,-1.d0)
   char_mat(5,7)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(5,8)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(5,9)=(-1.d0, 0.d0)
   char_mat(5,11)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(5,12)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(5,13)=( 0.d0,-1.d0)
   char_mat(5,14)=( 0.d0, 1.d0)
   char_mat(5,15)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(5,16)=CMPLX( dsq2,-dsq2,kind=DP)

   name_rap(6)="G_6- "
   char_mat(6,3)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(6,4)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(6,5)=( 0.d0,-1.d0)
   char_mat(6,6)=( 0.d0, 1.d0)
   char_mat(6,7)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(6,8)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(6,9)=(-1.d0, 0.d0)
   char_mat(6,11)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(6,12)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(6,13)=( 0.d0, 1.d0)
   char_mat(6,14)=( 0.d0,-1.d0)
   char_mat(6,15)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(6,16)=CMPLX( dsq2, dsq2,kind=DP)

   name_rap(7)="G_7- "
   char_mat(7,3)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(7,4)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(7,5)=( 0.d0, 1.d0)
   char_mat(7,6)=( 0.d0,-1.d0)
   char_mat(7,7)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(7,8)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(7,9)=(-1.d0, 0.d0)
   char_mat(7,11)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(7,12)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(7,13)=( 0.d0,-1.d0)
   char_mat(7,14)=( 0.d0, 1.d0)
   char_mat(7,15)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(7,16)=CMPLX(-dsq2, dsq2,kind=DP)

   name_rap(8)="G_8- "
   char_mat(8,3)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(8,4)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(8,5)=( 0.d0,-1.d0)
   char_mat(8,6)=( 0.d0, 1.d0)
   char_mat(8,7)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(8,8)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(8,9)=(-1.d0, 0.d0)
   char_mat(8,11)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(8,12)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(8,13)=( 0.d0, 1.d0)
   char_mat(8,14)=( 0.d0,-1.d0)
   char_mat(8,15)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(8,16)=CMPLX(-dsq2,-dsq2,kind=DP)

ELSEIF (code_group==19) THEN
!
! C_6h  NB: The signs of the characters of C3^2 -C3^2 C6^5 -C6^5 S6 -S6 S3 -S3
!           are changed with respect to Koster, Space groups and their 
!           representation. They match the table in Koster, Dimmock, 
!           Wheeler, Statz, Properties of the 32 point groups.
!           They are also consistent with the fact that C_6h=C_6xC_i
!
   nclass_ref=24
   name_class(3)="C6  "
   name_class(4)="-C6 "
   name_class(5)="C3  "
   name_class(6)="-C3 "
   name_class(7)="C2  "
   name_class(8)="-C2 "
   name_class(9)="C3^2"
   name_class(10)="-C3^2"
   name_class(11)="C6^5"
   name_class(12)="-C6^5"
   name_class(13)="i  "
   name_class(14)="-i "
   name_class(15)="S3^5"
   name_class(16)="-S3^5"
   name_class(17)="S6^5"
   name_class(18)="-S6^5"
   name_class(19)="s_h "
   name_class(20)="-s_h"
   name_class(21)="S6 "
   name_class(22)="-S6 "
   name_class(23)="S3 "
   name_class(24)="-S3 "

   nrap_ref=12

   name_rap(1)="G_7+ "
   char_mat(1,3)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(1,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(1,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(1,7)=( 0.d0, 1.0d0)
   char_mat(1,8)=( 0.d0,-1.0d0)
!   char_mat(1,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(1,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(1,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(1,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(1,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(1,14)=(-1.0d0, 0.d0)
   char_mat(1,15)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(1,16)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(1,17)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(1,18)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(1,19)=( 0.d0, 1.0d0)
   char_mat(1,20)=( 0.d0,-1.0d0)
!   char_mat(1,21)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(1,22)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(1,23)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(1,24)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,21)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(1,22)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(1,23)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(1,24)=CMPLX(-sqr3d2, 0.5d0,kind=DP)

   name_rap(2)="G_8+ "
   char_mat(2,3)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(2,4)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,7)=( 0.d0,-1.0d0)
   char_mat(2,8)=( 0.d0, 1.0d0)
!   char_mat(2,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(2,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(2,14)=(-1.0d0, 0.d0)
   char_mat(2,15)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(2,16)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(2,17)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,18)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,19)=( 0.d0,-1.0d0)
   char_mat(2,20)=( 0.d0, 1.0d0)
!   char_mat(2,21)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,22)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,23)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(2,24)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,21)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,22)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,23)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(2,24)=CMPLX(-sqr3d2,-0.5d0,kind=DP)


   name_rap(3)="G_9+ "
   char_mat(3,3)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(3,4)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(3,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,7)=( 0.d0,-1.0d0)
   char_mat(3,8)=( 0.d0, 1.0d0)
!   char_mat(3,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(3,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(3,14)=(-1.0d0, 0.d0)
   char_mat(3,15)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(3,16)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(3,17)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,18)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,19)=( 0.d0,-1.0d0)
   char_mat(3,20)=( 0.d0, 1.0d0)
!   char_mat(3,21)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,22)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,23)=CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(3,24)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,21)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,22)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,23)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(3,24)=CMPLX( sqr3d2,-0.5d0,kind=DP)

   name_rap(4)="G_10+ "
   char_mat(4,3)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(4,4)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(4,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(4,7)=( 0.d0, 1.0d0)
   char_mat(4,8)=( 0.d0,-1.0d0)
!   char_mat(4,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(4,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(4,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(4,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(4,14)=(-1.0d0, 0.d0)
   char_mat(4,15)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(4,16)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(4,17)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,18)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(4,19)=( 0.d0, 1.0d0)
   char_mat(4,20)=( 0.d0,-1.0d0)
!   char_mat(4,21)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(4,22)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(4,23)=CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(4,24)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,21)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,22)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,23)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(4,24)=CMPLX( sqr3d2, 0.5d0,kind=DP)

   name_rap(5)="G_11+ "
   char_mat(5,3)=( 0.d0, 1.0d0)
   char_mat(5,4)=( 0.d0,-1.0d0)
   char_mat(5,5)=(-1.d0, 0.d0 )
   char_mat(5,7)=( 0.d0,-1.0d0)
   char_mat(5,8)=( 0.d0, 1.0d0)
!   char_mat(5,10)=(-1.d0, 0.d0 )
!   char_mat(5,11)=( 0.d0, 1.0d0)
!   char_mat(5,12)=( 0.d0,-1.0d0)
   char_mat(5,9)=(-1.d0, 0.d0 )
   char_mat(5,11)=( 0.d0,-1.0d0)
   char_mat(5,12)=( 0.d0, 1.0d0)
   char_mat(5,14)=(-1.0d0, 0.d0)
   char_mat(5,15)=( 0.d0, 1.0d0)
   char_mat(5,16)=( 0.d0,-1.0d0)
   char_mat(5,17)=(-1.0d0, 0.d0)
   char_mat(5,19)=( 0.d0,-1.0d0)
   char_mat(5,20)=( 0.d0, 1.0d0)
!   char_mat(5,22)=(-1.0d0, 0.d0)
!   char_mat(5,23)=( 0.d0, 1.0d0)
!   char_mat(5,24)=( 0.d0,-1.0d0)
   char_mat(5,21)=(-1.0d0, 0.d0)
   char_mat(5,23)=( 0.d0,-1.0d0)
   char_mat(5,24)=( 0.d0, 1.0d0)

   name_rap(6)="G_12+ "
   char_mat(6,3)=( 0.d0,-1.0d0)
   char_mat(6,4)=( 0.d0, 1.0d0)
   char_mat(6,5)=(-1.d0, 0.d0 )
   char_mat(6,7)=( 0.d0, 1.0d0)
   char_mat(6,8)=( 0.d0,-1.0d0)
!   char_mat(6,10)=(-1.d0, 0.d0 )
!   char_mat(6,11)=( 0.d0,-1.0d0)
!   char_mat(6,12)=( 0.d0, 1.0d0)
   char_mat(6,9)=(-1.d0, 0.d0 )
   char_mat(6,11)=( 0.d0, 1.0d0)
   char_mat(6,12)=( 0.d0,-1.0d0)
   char_mat(6,14)=(-1.0d0, 0.d0)
   char_mat(6,15)=( 0.d0,-1.0d0)
   char_mat(6,16)=( 0.d0, 1.0d0)
   char_mat(6,17)=(-1.0d0, 0.d0)
   char_mat(6,19)=( 0.d0, 1.0d0)
   char_mat(6,20)=( 0.d0,-1.0d0)
!   char_mat(6,22)=(-1.0d0, 0.d0)
!   char_mat(6,23)=( 0.d0,-1.0d0)
!   char_mat(6,24)=( 0.d0, 1.0d0)
   char_mat(6,21)=(-1.0d0, 0.d0)
   char_mat(6,23)=( 0.d0, 1.0d0)
   char_mat(6,24)=( 0.d0,-1.0d0)


   name_rap(7)="G_7- "
   char_mat(7,3)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(7,4)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(7,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(7,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(7,7)=( 0.d0, 1.0d0)
   char_mat(7,8)=( 0.d0,-1.0d0)
!   char_mat(7,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(7,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(7,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(7,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(7,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(7,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(7,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(7,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(7,13)=(-1.0d0, 0.d0)
   char_mat(7,15)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(7,16)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(7,17)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(7,18)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(7,19)=( 0.d0,-1.0d0)
   char_mat(7,20)=( 0.d0, 1.0d0)
!   char_mat(7,21)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(7,22)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(7,23)=CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(7,24)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(7,21)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(7,22)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(7,23)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(7,24)=CMPLX( sqr3d2,-0.5d0,kind=DP)

   name_rap(8)="G_8- "
   char_mat(8,3)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(8,4)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(8,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(8,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(8,7)=( 0.d0,-1.0d0)
   char_mat(8,8)=( 0.d0, 1.0d0)
!   char_mat(8,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(8,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(8,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(8,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(8,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(8,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(8,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(8,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(8,13)=(-1.0d0, 0.d0)
   char_mat(8,15)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(8,16)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(8,17)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(8,18)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(8,19)=( 0.d0, 1.0d0)
   char_mat(8,20)=( 0.d0,-1.0d0)
!   char_mat(8,21)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(8,22)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(8,23)=CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(8,24)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(8,21)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(8,22)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(8,23)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(8,24)=CMPLX( sqr3d2, 0.5d0,kind=DP)


   name_rap(9)="G_9- "
   char_mat(9,3)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(9,4)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(9,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(9,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(9,7)=( 0.d0,-1.0d0)
   char_mat(9,8)=( 0.d0, 1.0d0)
!   char_mat(9,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(9,10)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(9,11)=CMPLX( sqr3d2,-0.5d0,kind=DP)
!   char_mat(9,12)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(9,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(9,10)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(9,11)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(9,12)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(9,13)=(-1.0d0, 0.d0)
   char_mat(9,15)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(9,16)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(9,17)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(9,18)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(9,19)=( 0.d0, 1.0d0)
   char_mat(9,20)=( 0.d0,-1.0d0)
!   char_mat(9,21)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(9,22)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(9,23)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
!   char_mat(9,24)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(9,21)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(9,22)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(9,23)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(9,24)=CMPLX(-sqr3d2, 0.5d0,kind=DP)

   name_rap(10)="G_10- "
   char_mat(10,3)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(10,4)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(10,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(10,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(10,7)=( 0.d0, 1.0d0)
   char_mat(10,8)=( 0.d0,-1.0d0)
!   char_mat(10,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(10,10)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(10,11)=CMPLX( sqr3d2, 0.5d0,kind=DP)
!   char_mat(10,12)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(10,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(10,10)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(10,11)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
   char_mat(10,12)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(10,13)=(-1.0d0, 0.d0)
   char_mat(10,15)=CMPLX( sqr3d2,-0.5d0,kind=DP)
   char_mat(10,16)=CMPLX(-sqr3d2, 0.5d0,kind=DP)
   char_mat(10,17)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(10,18)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(10,19)=( 0.d0,-1.0d0)
   char_mat(10,20)=( 0.d0, 1.0d0)
!   char_mat(10,21)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(10,22)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(10,23)=CMPLX(-sqr3d2,-0.5d0,kind=DP)
!   char_mat(10,24)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(10,21)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(10,22)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(10,23)=CMPLX( sqr3d2, 0.5d0,kind=DP)
   char_mat(10,24)=CMPLX(-sqr3d2,-0.5d0,kind=DP)

   name_rap(11)="G_11-"
   char_mat(11,3)=( 0.d0, 1.0d0)
   char_mat(11,4)=( 0.d0,-1.0d0)
   char_mat(11,5)=(-1.d0, 0.d0 )
   char_mat(11,7)=( 0.d0,-1.0d0)
   char_mat(11,8)=( 0.d0, 1.0d0)
!   char_mat(11,10)=(-1.d0, 0.d0 )
!   char_mat(11,11)=( 0.d0, 1.0d0)
!   char_mat(11,12)=( 0.d0,-1.0d0)
   char_mat(11,9)=(-1.d0, 0.d0 )
   char_mat(11,11)=( 0.d0,-1.0d0)
   char_mat(11,12)=( 0.d0, 1.0d0)
   char_mat(11,13)=(-1.0d0, 0.d0)
   char_mat(11,15)=( 0.d0,-1.0d0)
   char_mat(11,16)=( 0.d0, 1.0d0)
   char_mat(11,18)=(-1.0d0, 0.d0)
   char_mat(11,19)=( 0.d0, 1.0d0)
   char_mat(11,20)=( 0.d0,-1.0d0)
!   char_mat(11,21)=(-1.0d0, 0.d0)
!   char_mat(11,23)=( 0.d0,-1.0d0)
!   char_mat(11,24)=( 0.d0, 1.0d0)
   char_mat(11,22)=(-1.0d0, 0.d0)
   char_mat(11,23)=( 0.d0, 1.0d0)
   char_mat(11,24)=( 0.d0,-1.0d0)

   name_rap(12)="G_12-"
   char_mat(12,3) =( 0.d0,-1.0d0)
   char_mat(12,4) =( 0.d0, 1.0d0)
   char_mat(12,5) =(-1.d0, 0.d0 )
   char_mat(12,7) =( 0.d0, 1.0d0)
   char_mat(12,8) =( 0.d0,-1.0d0)
!   char_mat(12,10)=(-1.d0, 0.d0 )
!   char_mat(12,11)=( 0.d0,-1.0d0)
!   char_mat(12,12)=( 0.d0, 1.0d0)
   char_mat(12,9)=(-1.d0, 0.d0 )
   char_mat(12,11)=( 0.d0, 1.0d0)
   char_mat(12,12)=( 0.d0,-1.0d0)
   char_mat(12,13)=(-1.0d0, 0.d0)
   char_mat(12,15)=( 0.d0, 1.0d0)
   char_mat(12,16)=( 0.d0,-1.0d0)
   char_mat(12,18)=(-1.0d0, 0.d0)
   char_mat(12,19)=( 0.d0,-1.0d0)
   char_mat(12,20)=( 0.d0, 1.0d0)
!   char_mat(12,21)=(-1.0d0, 0.d0)
!   char_mat(12,23)=( 0.d0, 1.0d0)
!   char_mat(12,24)=( 0.d0,-1.0d0)
   char_mat(12,22)=(-1.0d0, 0.d0)
   char_mat(12,23)=( 0.d0,-1.0d0)
   char_mat(12,24)=( 0.d0, 1.0d0)



ELSEIF (code_group==20) THEN
!
! D_2h
!
   nclass_ref=10
   name_class(3)=" C2 "
   name_class1(3)="-C2 "
   name_class(4)=" C2'  "
   name_class1(4)="-C2'"
   name_class(5)="C2''  "
   name_class1(5)="-C2''"
   name_class(6)="i   "
   name_class(7)="-i  "
   name_class(8)=" s_v"
   name_class1(8)="-s_v"
   name_class(9)=" s_v'"
   name_class1(9)="-s_v'"
   name_class(10)="s_v''"
   name_class1(10)="-s_v''"

   nrap_ref=2

   name_rap(1)="G_5+ "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=( 0.d0, 0.d0)
   char_mat(1,4)=( 0.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 2.d0, 0.d0)
   char_mat(1,7)=(-2.d0, 0.d0)
   char_mat(1,8)=( 0.d0, 0.d0)
   char_mat(1,9)=( 0.d0, 0.d0)
   char_mat(1,10)=( 0.d0, 0.d0)

   name_rap(2)="G_5- "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=( 0.d0, 0.d0)
   char_mat(2,4)=( 0.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=(-2.d0, 0.d0)
   char_mat(2,7)=( 2.d0, 0.d0)
   char_mat(2,8)=( 0.d0, 0.d0)
   char_mat(2,9)=( 0.d0, 0.d0)
   char_mat(2,10)=( 0.d0, 0.d0)

ELSEIF (code_group==21) THEN
!
! D_3h
!
   nclass_ref=9
   name_class(3)="2C3 "
   name_class(4)="-2C3"
   name_class(5)=" 3C2'"
   name_class1(5)="-3C2'"
   name_class(6)=" s_h"
   name_class1(6)="-s_h'"
   name_class(7)="2S3 "
   name_class(8)="-2S3"
   name_class(9)=" 3s_v"
   name_class1(9)="-3s_v"

   nrap_ref=3

   name_rap(1)="G_7  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(1,8)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(1,9)=( 0.d0, 0.d0)

   name_rap(2)="G_8  "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,4)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=( 0.d0, 0.d0)
   char_mat(2,7)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(2,8)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(2,9)=( 0.d0, 0.d0)

   name_rap(3)="G_9  "
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=(-2.d0, 0.d0)
   char_mat(3,4)=( 2.d0, 0.d0)
   char_mat(3,5)=( 0.d0, 0.d0)
   char_mat(3,6)=( 0.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=( 0.d0, 0.d0)
   char_mat(3,9)=( 0.d0, 0.d0)

ELSEIF (code_group==22) THEN
!
! D_4h
!
   nclass_ref=14
   name_class(3)="2C4 "
   name_class(4)="-2C4"
   name_class(5)=" C2 "
   name_class1(5)="-C2 "
   name_class(6)=" 2C2'"
   name_class1(6)="-2C2' "
   name_class(7)="2C2''"
   name_class1(7)="-2C2''"
   name_class(8)="i   "
   name_class(9)="-i  "
   name_class(10)="2S4 "
   name_class(11)="-2S4"
   name_class(12)=" s_h"
   name_class1(12)="-s_h"
   name_class(13)=" 2s_v"
   name_class1(13)="-2s_v"
   name_class(14)=" 2s_d"
   name_class1(14)="-2s_d"

   nrap_ref=4

   name_rap(1)="G_6+  M_6+"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)
   char_mat(1,8)=( 2.d0, 0.d0)
   char_mat(1,9)=(-2.d0, 0.d0)
   char_mat(1,10)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,11)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,12)=( 0.d0, 0.d0)
   char_mat(1,13)=( 0.d0, 0.d0)
   char_mat(1,14)=( 0.d0, 0.d0)

   name_rap(2)="G_7+  M_7+"
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=( 0.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)
   char_mat(2,8)=( 2.d0, 0.d0)
   char_mat(2,9)=(-2.d0, 0.d0)
   char_mat(2,10)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,11)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,12)=( 0.d0, 0.d0)
   char_mat(2,13)=( 0.d0, 0.d0)
   char_mat(2,14)=( 0.d0, 0.d0)

   name_rap(3)="G_6-  M_6- "
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(3,4)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(3,5)=( 0.d0, 0.d0)
   char_mat(3,6)=( 0.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=(-2.d0, 0.d0)
   char_mat(3,9)=( 2.d0, 0.d0)
   char_mat(3,10)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(3,11)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(3,12)=( 0.d0, 0.d0)
   char_mat(3,13)=( 0.d0, 0.d0)
   char_mat(3,14)=( 0.d0, 0.d0)

   name_rap(4)="G_7-  M_7- "
   char_mat(4,1)=( 2.d0, 0.d0)
   char_mat(4,2)=(-2.d0, 0.d0)
   char_mat(4,3)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(4,4)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(4,5)=( 0.d0, 0.d0)
   char_mat(4,6)=( 0.d0, 0.d0)
   char_mat(4,7)=( 0.d0, 0.d0)
   char_mat(4,8)=(-2.d0, 0.d0)
   char_mat(4,9)=( 2.d0, 0.d0)
   char_mat(4,10)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(4,11)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(4,12)=( 0.d0, 0.d0)
   char_mat(4,13)=( 0.d0, 0.d0)
   char_mat(4,14)=( 0.d0, 0.d0)

ELSEIF (code_group==23) THEN
!
! D_6h
!
   nclass_ref=18
   name_class(3)="2C6 "
   name_class(4)="-2C6 "
   name_class(5)="2C3  "
   name_class(6)="-2C3 "
   name_class(7)="C2  "
   name_class1(7)="-C2 "
   name_class(8)=" 3C2'"
   name_class1(8)="-3C2'"
   name_class(9)="3C2''"
   name_class1(9)="-3C2''"
   name_class(10)=" i "
   name_class(11)="-i "
   name_class(12)="2S3"
   name_class(13)="-2S3"
   name_class(14)="2S6"
   name_class(15)="-2S6"
   name_class(16)=" s_h"
   name_class1(16)="-s_h"
   name_class(17)=" 3s_v"
   name_class1(17)="-3s_v"
   name_class(18)=" 3s_d"
   name_class1(18)="-3s_d"

   nrap_ref=6

   name_rap(1)="G_7+ "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(1,5)=( 1.d0, 0.d0)
   char_mat(1,6)=(-1.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)
   char_mat(1,8)=( 0.d0, 0.d0)
   char_mat(1,9)=( 0.d0, 0.d0)
   char_mat(1,10)=( 2.d0, 0.d0)
   char_mat(1,11)=(-2.d0, 0.d0)
   char_mat(1,12)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(1,13)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(1,14)=( 1.d0, 0.d0)
   char_mat(1,15)=(-1.d0, 0.d0)
   char_mat(1,16)=( 0.d0, 0.d0)
   char_mat(1,17)=( 0.d0, 0.d0)
   char_mat(1,18)=( 0.d0, 0.d0)

   name_rap(2)="G_8+   "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(2,5)=( 1.d0, 0.d0)
   char_mat(2,6)=(-1.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)
   char_mat(2,8)=( 0.d0, 0.d0)
   char_mat(2,9)=( 0.d0, 0.d0)
   char_mat(2,10)=( 2.d0, 0.d0)
   char_mat(2,11)=(-2.d0, 0.d0)
   char_mat(2,12)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(2,13)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(2,14)=( 1.d0, 0.d0)
   char_mat(2,15)=(-1.d0, 0.d0)
   char_mat(2,16)=( 0.d0, 0.d0)
   char_mat(2,17)=( 0.d0, 0.d0)
   char_mat(2,18)=( 0.d0, 0.d0)

   name_rap(3)="G_9+   "
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=( 0.d0, 0.d0)
   char_mat(3,4)=( 0.d0, 0.d0)
   char_mat(3,5)=(-2.d0, 0.d0)
   char_mat(3,6)=( 2.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=( 0.d0, 0.d0)
   char_mat(3,9)=( 0.d0, 0.d0)
   char_mat(3,10)=( 2.d0, 0.d0)
   char_mat(3,11)=(-2.d0, 0.d0)
   char_mat(3,12)=( 0.d0, 0.d0)
   char_mat(3,13)=( 0.d0, 0.d0)
   char_mat(3,14)=(-2.d0, 0.d0)
   char_mat(3,15)=( 2.d0, 0.d0)
   char_mat(3,16)=( 0.d0, 0.d0)
   char_mat(3,17)=( 0.d0, 0.d0)
   char_mat(3,18)=( 0.d0, 0.d0)

   name_rap(4)="G_7- "
   char_mat(4,1)=( 2.d0, 0.d0)
   char_mat(4,2)=(-2.d0, 0.d0)
   char_mat(4,3)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(4,4)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(4,5)=( 1.d0, 0.d0)
   char_mat(4,6)=(-1.d0, 0.d0)
   char_mat(4,7)=( 0.d0, 0.d0)
   char_mat(4,8)=( 0.d0, 0.d0)
   char_mat(4,9)=( 0.d0, 0.d0)
   char_mat(4,10)=(-2.d0, 0.d0)
   char_mat(4,11)=( 2.d0, 0.d0)
   char_mat(4,12)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(4,13)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(4,14)=(-1.d0, 0.d0)
   char_mat(4,15)=( 1.d0, 0.d0)
   char_mat(4,16)=( 0.d0, 0.d0)
   char_mat(4,17)=( 0.d0, 0.d0)
   char_mat(4,18)=( 0.d0, 0.d0)

   name_rap(5)="G_8-   "
   char_mat(5,1)=( 2.d0, 0.d0)
   char_mat(5,2)=(-2.d0, 0.d0)
   char_mat(5,3)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(5,4)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(5,5)=( 1.d0, 0.d0)
   char_mat(5,6)=(-1.d0, 0.d0)
   char_mat(5,7)=( 0.d0, 0.d0)
   char_mat(5,8)=( 0.d0, 0.d0)
   char_mat(5,9)=( 0.d0, 0.d0)
   char_mat(5,10)=(-2.d0, 0.d0)
   char_mat(5,11)=( 2.d0, 0.d0)
   char_mat(5,12)=CMPLX( sqrt3, 0.d0,kind=DP)
   char_mat(5,13)=CMPLX(-sqrt3, 0.d0,kind=DP)
   char_mat(5,14)=(-1.d0, 0.d0)
   char_mat(5,15)=( 1.d0, 0.d0)
   char_mat(5,16)=( 0.d0, 0.d0)
   char_mat(5,17)=( 0.d0, 0.d0)
   char_mat(5,18)=( 0.d0, 0.d0)

   name_rap(6)="G_9-   "
   char_mat(6,1)=( 2.d0, 0.d0)
   char_mat(6,2)=(-2.d0, 0.d0)
   char_mat(6,3)=( 0.d0, 0.d0)
   char_mat(6,4)=( 0.d0, 0.d0)
   char_mat(6,5)=(-2.d0, 0.d0)
   char_mat(6,6)=( 2.d0, 0.d0)
   char_mat(6,7)=( 0.d0, 0.d0)
   char_mat(6,8)=( 0.d0, 0.d0)
   char_mat(6,9)=( 0.d0, 0.d0)
   char_mat(6,10)=(-2.d0, 0.d0)
   char_mat(6,11)=( 2.d0, 0.d0)
   char_mat(6,12)=( 0.d0, 0.d0)
   char_mat(6,13)=( 0.d0, 0.d0)
   char_mat(6,14)=( 2.d0, 0.d0)
   char_mat(6,15)=(-2.d0, 0.d0)
   char_mat(6,16)=( 0.d0, 0.d0)
   char_mat(6,17)=( 0.d0, 0.d0)
   char_mat(6,18)=( 0.d0, 0.d0)

ELSEIF (code_group==24) THEN
!
! D_2d
!
   nclass_ref=7
   name_class(3)="2S4 "
   name_class(4)="-2S4"
   name_class(5)="C2  "
   name_class1(5)="-C2 "
   name_class(6)="2C2' "
   name_class1(6)="-2C2'"
   name_class(7)="2s_d "
   name_class1(7)="-2s_d"

   nrap_ref=2

   name_rap(1)="G_6  X_6  W_6"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,4)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=( 0.d0, 0.d0)

   name_rap(2)="G_7  X_7  W_7"
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,4)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=( 0.d0, 0.d0)
   char_mat(2,7)=( 0.d0, 0.d0)

ELSEIF (code_group==25) THEN
!
! D_3d
!
   nclass_ref=12
   name_class(3)="2C3 "
   name_class(4)="-2C3"
   name_class(5)="3C2'"
   name_class(6)="-3C2'"
   name_class(7)="i  "
   name_class(8)="-i  "
   name_class(9)="2S6 "
   name_class(10)="-2S6"
   name_class(11)="3s_v "
   name_class(12)="-3s_v"

   nrap_ref=6

   name_rap(1)="G_4+  L_6+"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=( 0.d0, 0.d0)
   char_mat(1,7)=( 2.d0, 0.d0)
   char_mat(1,8)=(-2.d0, 0.d0)
   char_mat(1,10)=(-1.d0, 0.d0)
   char_mat(1,11)=( 0.d0, 0.d0)
   char_mat(1,12)=( 0.d0, 0.d0)

   name_rap(2)="G_5+  L_4+"
   char_mat(2,2)=(-1.d0, 0.d0)
   char_mat(2,3)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 1.d0)
   char_mat(2,6)=( 0.d0,-1.d0)
   char_mat(2,8)=(-1.d0, 0.d0)
   char_mat(2,9)=(-1.d0, 0.d0)
   char_mat(2,11)=( 0.d0, 1.d0)
   char_mat(2,12)=( 0.d0,-1.d0)

   name_rap(3)="G_6+  L_5+"
   char_mat(3,2)=(-1.d0, 0.d0)
   char_mat(3,3)=(-1.d0, 0.d0)
   char_mat(3,5)=( 0.d0,-1.d0)
   char_mat(3,6)=( 0.d0, 1.d0)
   char_mat(3,8)=(-1.d0, 0.d0)
   char_mat(3,9)=(-1.d0, 0.d0)
   char_mat(3,11)=( 0.d0,-1.d0)
   char_mat(3,12)=( 0.d0, 1.d0)

   name_rap(4)="G_4-  L_6-"
   char_mat(4,1)=( 2.d0, 0.d0)
   char_mat(4,2)=(-2.d0, 0.d0)
   char_mat(4,4)=(-1.d0, 0.d0)
   char_mat(4,5)=( 0.d0, 0.d0)
   char_mat(4,6)=( 0.d0, 0.d0)
   char_mat(4,7)=(-2.d0, 0.d0)
   char_mat(4,8)=( 2.d0, 0.d0)
   char_mat(4,9)=(-1.d0, 0.d0)
   char_mat(4,11)=( 0.d0, 0.d0)
   char_mat(4,12)=( 0.d0, 0.d0)

   name_rap(5)="G_5-  L_4-"
   char_mat(5,2)=(-1.d0, 0.d0)
   char_mat(5,3)=(-1.d0, 0.d0)
   char_mat(5,5)=( 0.d0, 1.d0)
   char_mat(5,6)=( 0.d0,-1.d0)
   char_mat(5,7)=(-1.d0, 0.d0)
   char_mat(5,10)=(-1.d0, 0.d0)
   char_mat(5,11)=( 0.d0,-1.d0)
   char_mat(5,12)=( 0.d0, 1.d0)

   name_rap(6)="G_6-  L_5-"
   char_mat(6,2)=(-1.d0, 0.d0)
   char_mat(6,3)=(-1.d0, 0.d0)
   char_mat(6,5)=( 0.d0,-1.d0)
   char_mat(6,6)=( 0.d0, 1.d0)
   char_mat(6,7)=(-1.d0, 0.d0)
   char_mat(6,10)=(-1.d0, 0.d0)
   char_mat(6,11)=( 0.d0, 1.d0)
   char_mat(6,12)=( 0.d0,-1.d0)

ELSEIF (code_group==26) THEN
!
! S_4 ! This character table has been found to be working in at least one case. 
!       NB: The signs of the characters reported by Koster, Space groups and 
!            their representations are in the comment. They do not work at
!            least in one case. The characters in the tables in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups 
!            do not match neither those used here nor those of Koster.
!            They do not work at least in one case. 
!            Please report any problem that you might find with S_4. 
!
   nclass_ref=8
   name_class(3)=" S4^3"
   name_class(4)="-S4^3"
   name_class(5)=" C2 "
   name_class(6)="-C2 "
   name_class(7)=" S4 "
   name_class(8)="-S4 "

   nrap_ref=4

   name_rap(1)="G_5 "
   char_mat(1,3)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(1,4)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(1,5)=( 0.d0, 1.d0)
   char_mat(1,6)=( 0.d0,-1.d0)
!   char_mat(1,7)=CMPLX(-dsq2, dsq2,kind=DP)
!   char_mat(1,8)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,7)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(1,8)=CMPLX(-dsq2, dsq2,kind=DP)

   name_rap(2)="G_6 "
   char_mat(2,3)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(2,4)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(2,5)=( 0.d0,-1.d0)
   char_mat(2,6)=( 0.d0, 1.d0)
!   char_mat(2,7)=CMPLX(-dsq2,-dsq2,kind=DP)
!   char_mat(2,8)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,7)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(2,8)=CMPLX(-dsq2,-dsq2,kind=DP)

   name_rap(3)="G_7 "
   char_mat(3,3)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(3,4)=CMPLX( dsq2, dsq2,kind=DP)
   char_mat(3,5)=( 0.d0, 1.d0)
   char_mat(3,6)=( 0.d0,-1.d0)
!   char_mat(3,7)=CMPLX( dsq2,-dsq2,kind=DP)
!   char_mat(3,8)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,7)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(3,8)=CMPLX( dsq2,-dsq2,kind=DP)

   name_rap(4)="G_8 "
   char_mat(4,3)=CMPLX(-dsq2, dsq2,kind=DP)
   char_mat(4,4)=CMPLX( dsq2,-dsq2,kind=DP)
   char_mat(4,5)=( 0.d0,-1.d0)
   char_mat(4,6)=( 0.d0, 1.d0)
!   char_mat(4,7)=CMPLX( dsq2, dsq2,kind=DP)
!   char_mat(4,8)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,7)=CMPLX(-dsq2,-dsq2,kind=DP)
   char_mat(4,8)=CMPLX( dsq2, dsq2,kind=DP)

ELSEIF (code_group==27) THEN
!
! S_6    NB: The signs of the characters of the classes C3^2 -C3^2 S6 -S6 
!            are changed with respect to Koster, Space groups and 
!            their representations. They match the table in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups.
!
!
   nclass_ref=12
   name_class(3)="C3  "
   name_class(4)="-C3 "
   name_class(5)="C3^2"
   name_class(6)="-C3^2"
   name_class(7)=" i  "
   name_class(8)="-i  "
   name_class(9)="S6^5"
   name_class(10)="-S6^5"
   name_class(11)="S6"
   name_class(12)="-S6"

   nrap_ref=6

   name_rap(1)="G_4+"
   char_mat(1,3)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(1,4)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)

!   char_mat(1,5)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
!   char_mat(1,6)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(1,5)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(1,6)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)

   char_mat(1,8)=(-1.0d0, 0.d0 )
   char_mat(1,9)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(1,10)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(1,11)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
!   char_mat(1,12)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(1,11)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(1,12)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)

   name_rap(2)="G_5+"
   char_mat(2,3)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)

!   char_mat(2,5)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(2,6)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)

   char_mat(2,8)=(-1.0d0, 0.d0 )
   char_mat(2,9)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(2,10)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
!   char_mat(2,11)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(2,12)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(2,11)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(2,12)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)

   name_rap(3)="G_6+"
   char_mat(3,3)=(-1.0d0, 0.d0 )
!   char_mat(3,6)=(-1.0d0, 0.d0 )
   char_mat(3,5)=(-1.0d0, 0.d0 )
   char_mat(3,8)=(-1.0d0, 0.d0 )
   char_mat(3,9)=(-1.0d0, 0.d0 )
!   char_mat(3,12)=(-1.0d0, 0.d0 )
   char_mat(3,11)=(-1.0d0, 0.d0 )

   name_rap(4)="G_4-"
   char_mat(4,3)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(4,4)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(4,5)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
!   char_mat(4,6)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(4,5)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(4,6)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
   char_mat(4,7)=(-1.0d0, 0.d0 )
   char_mat(4,9)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
   char_mat(4,10)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
!   char_mat(4,11)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(4,12)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
   char_mat(4,11)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
   char_mat(4,12)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)

   name_rap(5)="G_5-"
   char_mat(5,3)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
   char_mat(5,4)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
!   char_mat(5,5)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(5,6)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(5,5)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
   char_mat(5,6)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
   char_mat(5,7)=(-1.0d0, 0.d0 )
   char_mat(5,9)=CMPLX(-0.5d0, sqr3d2 ,kind=DP)
   char_mat(5,10)=CMPLX( 0.5d0,-sqr3d2 ,kind=DP)
!   char_mat(5,11)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)
!   char_mat(5,12)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
   char_mat(5,11)=CMPLX(-0.5d0,-sqr3d2 ,kind=DP)
   char_mat(5,12)=CMPLX( 0.5d0, sqr3d2 ,kind=DP)

   name_rap(6)="G_6-"
   char_mat(6,3)=(-1.0d0, 0.d0 )
!   char_mat(6,6)=(-1.0d0, 0.d0 )
   char_mat(6,5)=(-1.0d0, 0.d0 )
   char_mat(6,7)=(-1.0d0, 0.d0 )
   char_mat(6,10)=(-1.0d0, 0.d0 )
!   char_mat(6,11)=(-1.0d0, 0.d0 )
   char_mat(6,12)=(-1.0d0, 0.d0 )


ELSEIF (code_group==28) THEN
!
! NB: The signs of the characters of the classes C3 -C3 C3^2 -C3^2 
!            are changed with respect to Koster, Space groups and 
!            their representations. They match the table in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups.

!
! T
!
   nclass_ref=7
   name_class(3)="3C2  "
   name_class1(3)="-3C2 "
   name_class(4)=" 4C3 "
   name_class(5)="-4C3 "
   name_class(6)=" 4C3'"
   name_class(7)="-4C3'"

   nrap_ref=3

   name_rap(1)="G_5 "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=( 0.d0, 0.d0)
   char_mat(1,5)=(-1.d0, 0.d0)
   char_mat(1,7)=(-1.d0, 0.d0)

   name_rap(2)="G_6 "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=( 0.d0, 0.d0)
!   char_mat(2,4)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(2,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)

   name_rap(3)="G_7 "
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=( 0.d0, 0.d0)
!   char_mat(3,4)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

ELSEIF (code_group==29) THEN
! NB: The signs of the characters of the some classes
!            are changed with respect to Koster, Space groups and 
!            their representations. They match the table in Koster, 
!            Dimmock, Wheeler, Statz, Properties of the 32 point groups.
!
! T_h
!
   nclass_ref=14
   name_class(3)=" 3C2 "
   name_class1(3)="-3C2 "
   name_class(4)=" 4C3 "
   name_class(5)="-4C3 "
   name_class(6)=" 4C3'"
   name_class(7)="-4C3'"
   name_class(8)="i   "
   name_class(9)="-i   "
   name_class(10)=" 3s_h"
   name_class1(10)="-3s_h"
   name_class(11)="4S6'"
   name_class(12)="-4S6'"
   name_class(13)=" 4S6 "
   name_class(14)="-4S6 "

   nrap_ref=6

   name_rap(1)="G_5+"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,3)=( 0.d0, 0.d0)
   char_mat(1,5)=(-1.d0, 0.d0)
   char_mat(1,7)=(-1.d0, 0.d0)
   char_mat(1,8)=( 2.d0, 0.d0)
   char_mat(1,9)=(-2.d0, 0.d0)
   char_mat(1,10)=( 0.d0, 0.d0)
   char_mat(1,12)=(-1.d0, 0.d0)
   char_mat(1,14)=(-1.d0, 0.d0)

   name_rap(2)="G_6+"
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,3)=( 0.d0, 0.d0)
!   char_mat(2,4)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(2,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(2,8)=( 2.d0, 0.d0)
   char_mat(2,9)=(-2.d0, 0.d0)
   char_mat(2,10)=( 0.d0, 0.d0)
!   char_mat(2,11)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(2,12)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(2,13)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(2,14)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,11)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,12)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(2,13)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,14)=CMPLX( 0.5d0, sqr3d2,kind=DP)

   name_rap(3)="G_7+"
   char_mat(3,1)=( 2.d0, 0.d0)
   char_mat(3,2)=(-2.d0, 0.d0)
   char_mat(3,3)=( 0.d0, 0.d0)
!   char_mat(3,4)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,8)=( 2.d0, 0.d0)
   char_mat(3,9)=(-2.d0, 0.d0)
   char_mat(3,10)=( 0.d0, 0.d0)
!   char_mat(3,11)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,12)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(3,13)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(3,14)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,11)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,12)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,13)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,14)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

   name_rap(4)="G_5-"
   char_mat(4,1)=( 2.d0, 0.d0)
   char_mat(4,2)=(-2.d0, 0.d0)
   char_mat(4,3)=( 0.d0, 0.d0)
   char_mat(4,5)=(-1.d0, 0.d0)
   char_mat(4,7)=(-1.d0, 0.d0)
   char_mat(4,8)=(-2.d0, 0.d0)
   char_mat(4,9)=( 2.d0, 0.d0)
   char_mat(4,10)=( 0.d0, 0.d0)
   char_mat(4,11)=(-1.d0, 0.d0)
   char_mat(4,13)=(-1.d0, 0.d0)


   name_rap(5)="G_6-"
   char_mat(5,1)=( 2.d0, 0.d0)
   char_mat(5,2)=(-2.d0, 0.d0)
   char_mat(5,3)=( 0.d0, 0.d0)
!   char_mat(5,4)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(5,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(5,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(5,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(5,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(5,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(5,8)=(-2.d0, 0.d0)
   char_mat(5,9)=( 2.d0, 0.d0)
   char_mat(5,10)=( 0.d0, 0.d0)
!   char_mat(5,11)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(5,12)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(5,13)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(5,14)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,11)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(5,12)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,13)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(5,14)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(6)="G_7-"
   char_mat(6,1)=( 2.d0, 0.d0)
   char_mat(6,2)=(-2.d0, 0.d0)
   char_mat(6,3)=( 0.d0, 0.d0)
!   char_mat(6,4)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(6,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(6,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
!   char_mat(6,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(6,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,5)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(6,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(6,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(6,8)=(-2.d0, 0.d0)
   char_mat(6,9)=( 2.d0, 0.d0)
   char_mat(6,10)=( 0.d0, 0.d0)
!   char_mat(6,11)=CMPLX(-0.5d0, sqr3d2,kind=DP)
!   char_mat(6,12)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
!   char_mat(6,13)=CMPLX( 0.5d0, sqr3d2,kind=DP)
!   char_mat(6,14)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,11)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(6,12)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,13)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(6,14)=CMPLX(-0.5d0, sqr3d2,kind=DP)

ELSEIF (code_group==30) THEN
!
! T_d
!
   nclass_ref=8
   name_class(3)="8C3  "
   name_class(4)="-8C3 "
   name_class(5)=" 3C2 "
   name_class1(5)="-3C2 "
   name_class(6)="6S4 "
   name_class(7)="-6S4 "
   name_class(8)="6s_d"
   name_class1(8)="-6s_d"

   nrap_ref=3

   name_rap(1)="G_6  P_6"
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,7)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,8)=( 0.d0, 0.d0)

   name_rap(2)="G_7  P_7"
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,4)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,7)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,8)=( 0.d0, 0.d0)

   name_rap(3)="G_8  P_8"
   char_mat(3,1)=( 4.d0, 0.d0)
   char_mat(3,2)=(-4.d0, 0.d0)
   char_mat(3,3)=(-1.d0, 0.d0)
   char_mat(3,5)=( 0.d0, 0.d0)
   char_mat(3,6)=( 0.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=( 0.d0, 0.d0)


ELSEIF (code_group==31) THEN
!
! O
!
   nclass_ref=8
   name_class(3)="8C3  "
   name_class(4)="-8C3 "
   name_class(5)=" 3C2"
   name_class1(5)="-3C2"
   name_class(6)="6C4  "
   name_class(7)="-6C4 "
   name_class(8)=" 6C2'"
   name_class1(8)="-6C2'"

   nrap_ref=3

   name_rap(1)="G_6  "
   char_mat(1,1)=(  2.d0, 0.d0)
   char_mat(1,2)=( -2.d0, 0.d0)
   char_mat(1,4)=( -1.d0, 0.d0)
   char_mat(1,5)=(  0.d0, 0.d0)
   char_mat(1,6)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,7)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,8)=(  0.d0, 0.d0)

   name_rap(2)="G_7  "
   char_mat(2,1)=(  2.d0, 0.d0)
   char_mat(2,2)=( -2.d0, 0.d0)
   char_mat(2,4)=( -1.d0, 0.d0)
   char_mat(2,5)=(  0.d0, 0.d0)
   char_mat(2,6)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,7)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,8)=(  0.d0, 0.d0)

   name_rap(3)="G_8  "
   char_mat(3,1)=(  4.d0, 0.d0)
   char_mat(3,2)=( -4.d0, 0.d0)
   char_mat(3,3)=( -1.d0, 0.d0)
   char_mat(3,5)=(  0.d0, 0.d0)
   char_mat(3,6)=(  0.d0, 0.d0)
   char_mat(3,7)=(  0.d0, 0.d0)
   char_mat(3,8)=(  0.d0, 0.d0)


ELSEIF (code_group==32) THEN
!
! O_h
!
   nclass_ref=16
   name_class(3)="8C3  "
   name_class(4)="-8C3 "
   name_class(5)=" 3C2"
   name_class1(5)="-3C2"
   name_class(6)="6C4  "
   name_class(7)="-6C4 "
   name_class(8)=" 6C2'"
   name_class1(8)="-6C2'"
   name_class(9)="i "
   name_class(10)="-i  "
   name_class(11)="8S6  "
   name_class(12)="-8S6 "
   name_class(13)=" 3s_h"
   name_class1(13)="-3s_h"
   name_class(14)="6S4  "
   name_class(15)="-6S4 "
   name_class(16)=" 6s_d"
   name_class1(16)="-6s_d"

   nrap_ref=6

   name_rap(1)="G_6+  "
   char_mat(1,1)=( 2.d0, 0.d0)
   char_mat(1,2)=(-2.d0, 0.d0)
   char_mat(1,4)=(-1.d0, 0.d0)
   char_mat(1,5)=( 0.d0, 0.d0)
   char_mat(1,6)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,7)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,8)=( 0.d0, 0.d0)
   char_mat(1,9)=( 2.d0, 0.d0)
   char_mat(1,10)=(-2.d0, 0.d0)
   char_mat(1,12)=(-1.d0, 0.d0)
   char_mat(1,13)=( 0.d0, 0.d0)
   char_mat(1,14)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(1,15)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(1,16)=( 0.d0, 0.d0)

   name_rap(2)="G_7+  "
   char_mat(2,1)=( 2.d0, 0.d0)
   char_mat(2,2)=(-2.d0, 0.d0)
   char_mat(2,4)=(-1.d0, 0.d0)
   char_mat(2,5)=( 0.d0, 0.d0)
   char_mat(2,6)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,7)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,8)=( 0.d0, 0.d0)
   char_mat(2,9)=( 2.d0, 0.d0)
   char_mat(2,10)=(-2.d0, 0.d0)
   char_mat(2,12)=(-1.d0, 0.d0)
   char_mat(2,13)=( 0.d0, 0.d0)
   char_mat(2,14)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(2,15)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(2,16)=( 0.d0, 0.d0)

   name_rap(3)="G_8+  "
   char_mat(3,1)=( 4.d0, 0.d0)
   char_mat(3,2)=(-4.d0, 0.d0)
   char_mat(3,3)=(-1.d0, 0.d0)
   char_mat(3,5)=( 0.d0, 0.d0)
   char_mat(3,6)=( 0.d0, 0.d0)
   char_mat(3,7)=( 0.d0, 0.d0)
   char_mat(3,8)=( 0.d0, 0.d0)
   char_mat(3,9)=( 4.d0, 0.d0)
   char_mat(3,10)=(-4.d0, 0.d0)
   char_mat(3,11)=(-1.d0, 0.d0)
   char_mat(3,13)=( 0.d0, 0.d0)
   char_mat(3,14)=( 0.d0, 0.d0)
   char_mat(3,15)=( 0.d0, 0.d0)
   char_mat(3,16)=( 0.d0, 0.d0)

   name_rap(4)="G_6-  "
   char_mat(4,1)=( 2.d0, 0.d0)
   char_mat(4,2)=(-2.d0, 0.d0)
   char_mat(4,4)=(-1.d0, 0.d0)
   char_mat(4,5)=( 0.d0, 0.d0)
   char_mat(4,6)=CMPLX(sqrt2, 0.d0,kind=DP)
   char_mat(4,7)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(4,8)=( 0.d0, 0.d0)
   char_mat(4,9)=(-2.d0, 0.d0)
   char_mat(4,10)=( 2.d0, 0.d0)
   char_mat(4,11)=(-1.d0, 0.d0)
   char_mat(4,13)=( 0.d0, 0.d0)
   char_mat(4,14)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(4,15)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(4,16)=( 0.d0, 0.d0)

   name_rap(5)="G_7-  "
   char_mat(5,1)=( 2.d0, 0.d0)
   char_mat(5,2)=(-2.d0, 0.d0)
   char_mat(5,4)=(-1.d0, 0.d0)
   char_mat(5,5)=( 0.d0, 0.d0)
   char_mat(5,6)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(5,7)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(5,8)=( 0.d0, 0.d0)
   char_mat(5,9)=(-2.d0, 0.d0)
   char_mat(5,10)=( 2.d0, 0.d0)
   char_mat(5,11)=(-1.d0, 0.d0)
   char_mat(5,13)=( 0.d0, 0.d0)
   char_mat(5,14)=CMPLX( sqrt2, 0.d0,kind=DP)
   char_mat(5,15)=CMPLX(-sqrt2, 0.d0,kind=DP)
   char_mat(5,16)=( 0.d0, 0.d0)


   name_rap(6)="G_8-  "
   char_mat(6,1)=( 4.d0, 0.d0)
   char_mat(6,2)=(-4.d0, 0.d0)
   char_mat(6,3)=(-1.d0, 0.d0)
   char_mat(6,5)=( 0.d0, 0.d0)
   char_mat(6,6)=( 0.d0, 0.d0)
   char_mat(6,7)=( 0.d0, 0.d0)
   char_mat(6,8)=( 0.d0, 0.d0)
   char_mat(6,9)=(-4.d0, 0.d0)
   char_mat(6,10)=( 4.d0, 0.d0)
   char_mat(6,12)=(-1.d0, 0.d0)
   char_mat(6,13)=( 0.d0, 0.d0)
   char_mat(6,14)=( 0.d0, 0.d0)
   char_mat(6,15)=( 0.d0, 0.d0)
   char_mat(6,16)=( 0.d0, 0.d0)
ELSE
   CALL errore('set_irr_rap_so','code number not allowed',1)
END IF

RETURN
END SUBROUTINE set_irr_rap_so

!--------------------------------------------------------------------------
FUNCTION is_complex_so(code)
!--------------------------------------------------------------------------
! This function receives a code of the group and provide .true. or 
! .false. if the double group HAS or HAS NOT complex irreducible 
! representations.
! The order is the following:
!
!   1  "C_1 " F    11 "D_6 " F    21 "D_3h" F    31 "O   " F
!   2  "C_i " F    12 "C_2v" F    22 "D_4h" F    32 "O_h " F 
!   3  "C_s " T    13 "C_3v" T    23 "D_6h" F 
!   4  "C_2 " T    14 "C_4v" F    24 "D_2d" F  
!   5  "C_3 " T    15 "C_6v" F    25 "D_3d" T
!   6  "C_4 " T    16 "C_2h" T    26 "S_4 " T
!   7  "C_6 " T    17 "C_3h" T    27 "S_6 " T
!   8  "D_2 " F    18 "C_4h" T    28 "T   " T
!   9  "D_3 " T    19 "C_6h" T    29 "T_h " T
!   10 "D_4 " F    20 "D_2h" F    30 "T_d " F
!
IMPLICIT NONE

INTEGER :: code
LOGICAL :: is_complex_so

LOGICAL :: complex_aux(32)

data complex_aux  / .FALSE., .FALSE., .TRUE.,  .TRUE.,  .TRUE. , &
                    .TRUE. , .TRUE. , .FALSE., .TRUE.,  .FALSE., &
                    .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE., &
                    .TRUE.,  .TRUE. , .TRUE.,  .TRUE.,  .FALSE., &
                    .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  &
                    .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., &
                    .FALSE., .FALSE.  /

IF (code < 1 .OR. code > 32 ) CALL errore('is_complex', &
                                          'code is out of range',1)

is_complex_so= complex_aux(code)

RETURN
END FUNCTION is_complex_so


!
!----------------------------------------------------------------------------
SUBROUTINE write_group_info(flag)
!----------------------------------------------------------------------------
!
! This routine writes on output the main information on the point group
! If flag is .false. writes only the character table. If flag is .true.
! writes also the elements of each class.
!
!
USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
                                 char_mat, name_rap, name_class, gname,      &
                                 elem_name
USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
                                 which_irr_so, char_mat_so, name_rap_so, &
                                 name_class_so, d_spin, name_class_so1,  &
                                 elem_name_so
USE rap_point_group_is,   ONLY : code_group_is, gname_is
USE spin_orb,             ONLY : domag
USE noncollin_module,     ONLY : noncolin
USE io_global,            ONLY : stdout

IMPLICIT NONE

INTEGER :: iclass, irot, i, idx, irap
LOGICAL :: is_complex, is_complex_so, flag

IF (noncolin) THEN
   IF (domag) THEN
      WRITE(stdout,'(/,5x,"the magnetic double point group is ", &
                     & a11," [",a11,"]")') &
                      gname, gname_is
      WRITE(stdout,'(5x,"using the double point group ",a11)') &
                      gname_is
   ELSE
      WRITE(stdout,'(/,5x,"double point group ",a11)') gname
   END IF
   WRITE(stdout,'(5x, "there are", i3," classes and",i3, &
                     &   " irreducible representations")') nclass, nrap
ELSE
   WRITE(stdout,'(/,5x,"point group ",a11)') gname
   WRITE(stdout,'(5x, "there are", i3," classes")') nclass
ENDIF
WRITE(stdout,'(5x, "the character table:")')
IF (noncolin) THEN
   WRITE(stdout,'(/,7x,12(a5,1x))') (name_class_so(irot), &
                                     irot=1,MIN(12,nclass))
   WRITE(stdout,'(7x,12(a5,1x))') (name_class_so1(irot), &
                                     irot=1,MIN(12,nclass))

   DO iclass=1,nrap
      WRITE(stdout,'(a5,12f6.2)') name_rap_so(iclass), &
               (REAL(char_mat_so(iclass,irot)), irot=1,MIN(nclass,12))
   END DO
   IF (nclass > 12 ) THEN
      WRITE(stdout,'(/,7x,12(a5,1x))') (name_class_so(irot), &
                                        irot=13,nclass)
      WRITE(stdout,'(7x,12(a5,1x))') (name_class_so1(irot), &
                                        irot=13,nclass)
      DO iclass=1,nrap
         WRITE(stdout,'(a5,12f6.2)') name_rap_so(iclass), &
                   (REAL(char_mat_so(iclass,irot)), irot=13,nclass)
      END DO
   END IF
   idx=code_group
   IF (noncolin.and.domag) idx=code_group_is
   IF (is_complex_so(idx)) THEN
      WRITE(stdout,'(/,5x,"imaginary part")')
      WRITE(stdout,'(/,7x,12(a5,1x))') (name_class_so(irot), &
                                        irot=1,MIN(12,nclass))
      WRITE(stdout,'(7x,12(a5,1x))') (name_class_so1(irot), &
                                        irot=1,MIN(12,nclass))
      DO iclass=1,nrap
         WRITE(stdout,'(a5,12f6.2)') name_rap_so(iclass), &
              (AIMAG(char_mat_so(iclass,irot)),irot=1, MIN(nclass,12))
      END DO
      IF (nclass > 12 ) THEN
         WRITE(stdout,'(/,7x,12(a5,1x))') (name_class_so(irot), &
                                        irot=13,nclass)
         WRITE(stdout,'(7x,12(a5,1x))') (name_class_so1(irot), &
                                        irot=13,nclass)
         DO iclass=1,nrap
            WRITE(stdout,'(a5,12f6.2)') name_rap_so(iclass), &
                 (AIMAG(char_mat_so(iclass,irot)),irot=13, nclass)
         END DO
      END IF
   END IF
   IF (flag) THEN
      WRITE(stdout,'(/5x, "the symmetry operations in each class and &
                          &the name of the first element:",/)')
      DO irap = 1, nclass
         DO iclass=1,nclass
            IF ( which_irr_so(iclass) /= irap ) CYCLE
            WRITE(stdout,'(5x,2a5,12i5)') &
                              name_class_so(which_irr_so(iclass)), &
                              name_class_so1(which_irr_so(iclass)), &
               (elem_so(i,iclass)*has_e(i,iclass), i=1,nelem_so(iclass))
            WRITE(stdout,'(10x,a)') elem_name_so(1,iclass)
         END DO
      ENDDO
   ENDIF
ELSE
   WRITE(stdout,'(/,7x,12(a5,1x))') (name_class(irot),irot=1,nclass)
   DO iclass=1,nclass
      WRITE(stdout,'(a5,12f6.2)') name_rap(iclass), &
         (REAL(char_mat(iclass,irot)),irot=1,nclass)
   ENDDO
   idx=code_group
   IF (noncolin.and.domag) idx=code_group_is
   IF (is_complex(idx)) THEN
      WRITE(stdout,'(5x,"imaginary part")')
      DO iclass=1,nclass
         WRITE(stdout,'(a5,12f6.2)') name_rap(iclass), &
              (AIMAG(char_mat(iclass,irot)),irot=1,nclass)
      ENDDO
   ENDIF
   IF (flag) THEN
      WRITE(stdout,'(/5x, "the symmetry operations in each class and &
                           &the name of the first element:",/)')
      DO irap = 1, nclass
         DO iclass=1,nclass
            IF (which_irr(iclass)/=irap) CYCLE
            WRITE(stdout,'(5x,a5,12i5)') name_class(which_irr(iclass)), &
               (elem(i,iclass), i=1,nelem(iclass))
!
!    The name of the first element of each class is written on output
!
            WRITE(stdout,'(10x,a)') elem_name(1,iclass)
         ENDDO
      ENDDO
   END IF
END IF
RETURN
END SUBROUTINE write_group_info

!---------------------------------------------------------------------------
SUBROUTINE find_u(s,u)
!---------------------------------------------------------------------------
!
!  This subroutine receives as input a 3x3 rotation matrix s, and gives
!  as output the matrix u which represents the same rotation in the spin
!  space. Only one of the two u matrices is given. See below for the
!  definition of the sign.
!
USE kinds,   ONLY : DP
USE constants, ONLY : pi

IMPLICIT NONE
REAL(DP) :: s(3,3)

COMPLEX(DP) :: u(2,2)

REAL(DP), PARAMETER :: eps=1.d-8
REAL(DP)  :: det, saux(3,3), ax(3), angle, cosa, sina, angle_rot
!
!  For consistency check uncomment here
!
!COMPLEX(DP) :: a, as, b, bs
!REAL(DP) :: r(3,3), r1(3,3), diff

det = s(1,1) * ( s(2,2) * s(3,3) - s(3,2) * s(2,3) )-   &
      s(1,2) * ( s(2,1) * s(3,3) - s(3,1) * s(2,3) )+   &
      s(1,3) * ( s(2,1) * s(3,2) - s(3,1) * s(2,2) )
!
!  inversion has no effect in spin space, so improper rotations are 
!  multiplied by inversion
!
IF (ABS(det+1.d0) < eps) THEN
   saux=-s
ELSE
   saux=s
ENDIF
!
! Check for identity or inversion
!
IF ((ABS(saux(1,1)-1.d0) < eps).AND. &
    (ABS(saux(2,2)-1.d0) < eps).AND. &
    (ABS(saux(3,3)-1.d0) < eps).AND. &
    (ABS(saux(1,2)) < eps).AND.(ABS(saux(2,1)) < eps) &
.AND.(ABS(saux(2,3)) < eps).AND. &
    (ABS(saux(3,2)) < eps).AND.(ABS(saux(1,3)) < eps) &
.AND.(ABS(saux(3,1)) < eps)) THEN
   u(1,1)=(1.d0,0.d0)
   u(1,2)=(0.d0,0.d0)
   u(2,1)=(0.d0,0.d0)
   u(2,2)=(1.d0,0.d0)
   RETURN
ENDIF
!
!   Find the rotation axis and the rotation angle
!
CALL versor(saux,ax)
angle=angle_rot(saux)
!write(6,'(5x,"find u",3f12.5,5x,f12.5)') ax(1), ax(2), ax(3), angle
angle=0.5d0*angle*pi/180.d0
cosa=COS(angle)
sina=SIN(angle)
!write(6,'(2f12.5)') cosa, sina
!
!  set the spin space rotation matrix elements
!
u(1,1)=CMPLX(cosa,-ax(3)*sina,kind=DP)
u(1,2)=CMPLX(-ax(2)*sina, -ax(1)*sina,kind=DP)
u(2,1)=-CONJG(u(1,2))
u(2,2)=CONJG(u(1,1))
!
!  To each 3x3 rotation one can associate two 2x2 rotation matrices in spin
!  space. This function returns the U matrix with positive cosa term
!
IF (cosa < -eps ) u=-u

!IF (ABS(cosa) < eps) THEN
!
!  Special case when cosa=0. For this rotation we must take the negative sign.
!
!   IF (ax(1)*ax(3) < -eps) u=-u
!ENDIF
!
!   Here compute the 3x3 rotation matrix starting form the axis, angle
!   or from the rotation in spin space for consistency check.
!
!angle=angle*2.d0
!cosa=COS(angle)
!sina=SIN(angle)
!r1(1,1)=1.d0+(1.d0-cosa)*(ax(1)**2-1)
!r1(1,2)=-ax(3)*sina+(1.d0-cosa)*ax(1)*ax(2)
!r1(1,3)=ax(2)*sina+(1.d0-cosa)*ax(1)*ax(3)
!r1(2,1)=ax(3)*sina+(1.d0-cosa)*ax(1)*ax(2)
!r1(2,2)=1.d0+(1.d0-cosa)*(ax(2)**2-1)
!r1(2,3)=-ax(1)*sina+(1.d0-cosa)*ax(2)*ax(3)
!r1(3,1)=-ax(2)*sina+(1.d0-cosa)*ax(1)*ax(3)
!r1(3,2)=ax(1)*sina+(1.d0-cosa)*ax(2)*ax(3)
!r1(3,3)=1.d0+(1.d0-cosa)*(ax(3)**2-1)

!a=u(1,1)
!as=u(2,2)
!b=u(1,2)
!bs=-u(2,1)

!r(1,1)=0.5d0*(a**2+as**2-b**2-bs**2)
!r(1,2)=0.5d0*(0.d0,1.d0)*(as**2+bs**2-a**2-b**2)
!r(1,3)=-(a*b+as*bs)

!r(2,1)=-0.5d0*(0.d0,1.d0)*(as**2-a**2+b**2-bs**2)
!r(2,2)=0.5d0*(a**2+b**2+as**2+bs**2)
!r(2,3)=(0.d0,1.d0)*(as*bs-a*b)

!r(3,1)=(bs*a+as*b)
!r(3,2)=(0.d0,1.d0)*(as*b-bs*a)
!r(3,3)=(a*as-b*bs)

!diff=ABS(r(1,1)-saux(1,1))+ &
!     ABS(r(1,2)-saux(1,2))+ &
!     ABS(r(1,3)-saux(1,3))+ &
!     ABS(r(2,1)-saux(2,1))+ &
!     ABS(r(2,2)-saux(2,2))+ &
!     ABS(r(2,3)-saux(2,3))+ &
!     ABS(r(3,1)-saux(3,1))+ &
!     ABS(r(3,2)-saux(3,2))+ &
!     ABS(r(3,3)-saux(3,3))


!write(6,*) diff
!write(6,'(3f15.5)') r1(1,1),r1(1,2),r1(1,3)
!write(6,'(3f15.5)') r1(2,1),r1(2,2),r1(2,3)
!write(6,'(3f15.5)') r1(3,1),r1(3,2),r1(3,3)
!write(6,*)
!write(6,'(3f15.5)') r(1,1),r(1,2),r(1,3)
!write(6,'(3f15.5)') r(2,1),r(2,2),r(2,3)
!write(6,'(3f15.5)') r(3,1),r(3,2),r(3,3)
!write(6,*)
!write(6,'(4f15.5)') u(1,1),u(1,2)
!write(6,'(4f15.5)') u(2,1),u(2,2)
!
RETURN
END SUBROUTINE find_u

!-----------------------------------------------------------------------------
FUNCTION set_e(hase,ind)
!-----------------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: set_e, hase, ind

IF (hase==-1) THEN
   set_e=ind+1
ELSE
   set_e=ind
ENDIF

RETURN
END FUNCTION set_e

!-----------------------------------------------------------------------------
SUBROUTINE check_tgroup(nsym,a,b)
!-----------------------------------------------------------------------------
!
!  This subroutine receives a set of 2x2 and 3x3 rotation matrices and 
!  checks if they are a group. 
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nsym
COMPLEX(DP), INTENT(IN) :: a(2,2,96)
REAL(DP), INTENT(IN) :: b(3,3,nsym)
REAL(DP) ::   d(3,3), b1(3,3), b2(3,3), b3(3,3)
COMPLEX(dp):: c(2,2), a1(2,2), a2(2,2), a3(2,2)
INTEGER :: done
LOGICAL :: compare_mat_so

INTEGER :: i, j, k

DO i=1,nsym
   a1(:,:)=a(:,:,i)
   b1(:,:)=b(:,:,i)
   DO j=1,nsym
      a2(:,:)=a(:,:,j)
      b2(:,:)=b(:,:,j)
      c=MATMUL(a1,a2)
      d=MATMUL(b1,b2)
      done=0
      do k=1,nsym
         a3(:,:)=a(:,:,k)
         b3(:,:)=b(:,:,k)
         IF (compare_mat_so(d,c,b3,a3)) THEN 
            done=done+1
         ENDIF
      ENDDO
      IF (done.ne.1) write(6,*) 'problem, i,j',i,j
   END DO
END DO 
RETURN
END SUBROUTINE check_tgroup

SUBROUTINE set_class_el_name_so(nsym,sname,has_e,nclass,&
                                     nelem_so,elem_so,elem_name_so)

IMPLICIT NONE
INTEGER :: nsym
CHARACTER(LEN=45) :: sname(nsym)
CHARACTER(LEN=55) :: elem_name_so(12,24)
INTEGER :: nclass, nelem_so(24), elem_so(12,24), has_e(12,24)

INTEGER :: iclass, ielem

DO iclass=1,nclass
   DO ielem=1,nelem_so(iclass)
      elem_name_so(ielem,iclass)=sname(elem_so(ielem,iclass))
      IF (has_e(ielem,iclass)==-1) elem_name_so(ielem,iclass)=&
                          TRIM(elem_name_so(ielem,iclass)) // ' E'
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_class_el_name_so
