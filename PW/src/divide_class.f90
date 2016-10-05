!
! Copyright (C) 2006-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------------
SUBROUTINE divide_class(code_group,nrot,smat,nclass,nelem,elem,which_irr)
!-----------------------------------------------------------------------------
!
! This subroutine receives as input a set of nrot 3x3 matrices smat, which
! are assumed to be the operations of the point group given by code_group.
! smat are in cartesian coordinates.
! This routine divides the group in classes and find:
!
! nclass         the number of classes of the group
! nelem(iclass)  for each class, the number of elements of the class
! elem(i,iclass) 1<i<nelem(iclass) for each class tells which matrices 
!                smat belong to that class
! which_irr(iclass) for each class gives the position of that class in the
!                character table associated with the group and provided
!                by the routine set_irr_rap. NB: changing the order of
!                the elements in the character table must be reflected in 
!                a change to which_irr. Presently the character tables 
!                are those given by P.W. Atkins, M.S. Child, and 
!                C.S.G. Phillips, "Tables for group theory". 
!                Several equivalent names for the irreducible representation
!                are given. D, G, L, S are used for Delta, Gamma, Lambda 
!                and Sigma.
!

USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: & 
          code_group,  &   ! The code of the point group
          nrot,        &   ! The number of symmetry operation
          nclass,      &   ! The number of classes
          nelem(12),   &   ! The elements of each class 
          elem(8,12),  &   ! Which elements in the smat list for each class
          which_irr(12)    ! See above 

REAL(DP) :: smat(3,3,nrot), cmat(3,3), ax(3), bx(3), cx(3), ars

INTEGER :: done(48), irot, jrot, krot, iclass, i, other, other1
INTEGER :: tipo_sym, ipol, axis, axis1, axis2, ts, iax, ibx, icx, aclass, &
           bclass, cclass, imax, imbx, imcx, amclass, bmclass, cmclass, ind2(3)
REAL(DP), PARAMETER :: eps = 1.d-7
REAL(DP) :: angle_rot, angle_rot_s, angle_vectors, ax_save(3,2:4)
LOGICAL :: compare_mat, is_axis, is_parallel, isok, isok1, done_ax(6)
!
! Divide the group in classes.
!
nclass=0
nelem=0
done=0
DO irot=1,nrot
   IF (done(irot)==0) THEN
      nclass=nclass+1
      DO jrot=1,nrot
         CALL coniug_mat(smat(1,1,jrot),smat(1,1,irot),cmat)
         DO krot=1,nrot
            IF (compare_mat(cmat,smat(1,1,krot)).AND.done(krot)==0) THEN
               nelem(nclass)=nelem(nclass)+1
               elem(nelem(nclass),nclass)=krot
               done(krot)=1 
            ENDIF 
         ENDDO
      ENDDO
   ENDIF
ENDDO
!
!  For each class we should now decide which_irr. This depends on the group
!  and on the tables of characters of the irreducible representations,
!  so we must make different things for different groups.
!
which_irr(1)=1
IF (code_group==1) THEN
   IF (nclass /= 1) CALL errore('divide_class','Wrong classes for C_1',1)
!
!  C_1 
!
ELSEIF (code_group==2.OR.code_group==3.OR.code_group==4) THEN
!
!  C_i, C_s, C_2
!
   IF (nclass /= 2) &
        CALL errore('divide_class','Wrong classes for C_i, C_s or C_2',1)
   which_irr(2)=2
ELSEIF (code_group==5) THEN
!
!  C_3
!
! The function angle_rot(smat) provides the rotation angle of the matrix smat
!
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for C_3',1)
   DO iclass=2,nclass
      IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
         which_irr(iclass)=2
      ELSE
         which_irr(iclass)=3
      ENDIF
   ENDDO

ELSEIF (code_group==6) THEN
!
!  C_4
!
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-90.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-270.d0)<eps) THEN
            which_irr(iclass)=4
         ELSE
            CALL errore('divide_class','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==7) THEN
!
!  C_6
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=5
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=6
         ELSE
            CALL errore('divide_class','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO

ELSEIF (code_group==8) THEN
!
!  D_2  
!
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for D_2',1)

   CALL versor(smat(1,1,elem(1,2)),ax)
   CALL which_c2(ax,iax)
   CALL versor(smat(1,1,elem(1,3)),bx)
   CALL which_c2(bx,ibx)
   CALL versor(smat(1,1,elem(1,4)),cx)
   CALL which_c2(cx,icx)

   CALL is_d2(iax, ibx, icx, ind2)

   which_irr(2)=ind2(1) + 1
   which_irr(3)=ind2(2) + 1
   which_irr(4)=ind2(3) + 1

ELSEIF (code_group==9) THEN
!
!  D_3
!
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for D_3',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==3) THEN
         which_irr(iclass)=2
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==10) THEN
!
!  D_4
!
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for D_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         axis1=MOD(ipol,3)+1
         axis2=MOD(ipol+1,3)+1
         IF (axis==0) call errore('divide_class','unknown D_4 axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=3
         ELSEIF (is_axis(ax,axis1).or.is_axis(ax,axis2)) THEN
            which_irr(iclass)=4
         ELSE
            which_irr(iclass)=5
         END IF
      ELSEIF (ts.ne.3) THEN
         CALL errore('divide_class','wrong sym_type',1)
      END IF
   END DO
ELSEIF (code_group==11) THEN
!
!  D_6
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps) ) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         ENDIF
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,3)) THEN
            which_irr(iclass)=4
         ELSE 
            CALL which_c2(ax, iax)
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=5
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=6
            ELSE
               CALL errore('divide_class','D_6 C_2 axis not recognized',1)
            END IF
         END IF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      END IF
   END DO
ELSEIF (code_group==12) THEN
!
!  C_2v
!
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_2v',1)
   iax=0
   ibx=0
   icx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2( ax, iax) 
         which_irr(iclass)=2
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
      which_irr(bclass)=3
      which_irr(cclass)=4
   ELSE
      CALL is_c2v(iax, icx, ibx, isok1)
      IF (.NOT.isok1) CALL errore('divide_class','problem with C_2v',1)
      which_irr(bclass)=4
      which_irr(cclass)=3
   ENDIF
ELSEIF (code_group==13) THEN
!
!  C_3v
!
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for C_3v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSEIF (ts==5) THEN
         which_irr(iclass)=3
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==14) THEN
!
!  C_4v
!
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for C_4v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax, iax)
         IF (iax < 4) THEN
!
!   x, y, or z
!
            which_irr(iclass)=4
         ELSE
!
!  all other cases
!
            which_irr(iclass)=5
         ENDIF
      ENDIF
   ENDDO

ELSEIF (code_group==15) THEN
!
!  C_6v
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_6v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         ENDIF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax, iax)
         IF (iax==2 .OR. iax==12 .OR. iax==13) THEN
            which_irr(iclass)=6
         ELSEIF (iax==1 .OR. iax==10 .OR. iax==11) THEN
            which_irr(iclass)=5
         ELSE
            CALL errore('divide_class','C_6v mirror not recognized',1)
         ENDIF
      ENDIF
   ENDDO
ELSEIF (code_group==16) THEN
!
!  C_2h
!
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_2h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=2
      ELSEIF (ts==2) THEN
         which_irr(iclass)=3
      ELSEIF (ts==5) THEN
         which_irr(iclass)=4
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==17) THEN
!
!  C_3h
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSEIF (ts==5) THEN
         which_irr(iclass)=4
      ELSEIF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=5
         ELSE
            which_irr(iclass)=6
         END IF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==18) THEN
!
!  C_4h
!
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for C_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (angle_rot(smat(1,1,elem(1,iclass)))-90.d0<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==2) THEN
         which_irr(iclass)=5
      ELSEIF (ts==5) THEN
         which_irr(iclass)=7
      ELSEIF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-90.d0)<eps) THEN
            which_irr(iclass)=8
         ELSE
            which_irr(iclass)=6
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==19) THEN
!
!  C_6h
!
   IF (nclass /= 12) CALL errore('divide_class','Wrong classes for C_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=5
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=6
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==2) THEN
         which_irr(iclass)=7
      ELSEIF (ts==5) THEN
         which_irr(iclass)=10
      ELSEIF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=11
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=12
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=8
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=9
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==20) THEN
!
!  D_2h
!
!  mirror_axis gives the normal to the mirror plane
!
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for D_2h',1)
   iax=0
   ibx=0
   icx=0
   imax=0
   imbx=0
   imcx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
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
            CALL errore('divide_class','D_2h too many C_2 axis',1)
         ENDIF 
      ELSEIF (ts==2) THEN
         which_irr(iclass)=5
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
            CALL errore('divide_class','D_2h too many mirrors',1)
         ENDIF 
      ELSE
         CALL errore('divide_class','D_2h operation not recognized',1)
      ENDIF
   ENDDO

   CALL is_d2(iax, ibx, icx, ind2)

   which_irr(aclass)=ind2(1) + 1 
   which_irr(bclass)=ind2(2) + 1 
   which_irr(cclass)=ind2(3) + 1  

   IF (imax==iax) which_irr(amclass) = which_irr(aclass) + 4
   IF (imax==ibx) which_irr(amclass) = which_irr(bclass) + 4
   IF (imax==icx) which_irr(amclass) = which_irr(cclass) + 4
   IF (imbx==iax) which_irr(bmclass) = which_irr(aclass) + 4
   IF (imbx==ibx) which_irr(bmclass) = which_irr(bclass) + 4
   IF (imbx==icx) which_irr(bmclass) = which_irr(cclass) + 4
   IF (imcx==iax) which_irr(cmclass) = which_irr(aclass) + 4
   IF (imcx==ibx) which_irr(cmclass) = which_irr(bclass) + 4
   IF (imcx==icx) which_irr(cmclass) = which_irr(cclass) + 4

ELSEIF (code_group==21) THEN
!
!  D_3h
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)>1) THEN
            which_irr(iclass)=6
         ELSE 
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==6) THEN
         which_irr(iclass)=5
      END IF
   END DO
ELSEIF (code_group==22) THEN
!
!  D_4h
!
!
!  First search the order 4 axis
!
   IF (nclass /= 10) CALL errore('divide_class','Wrong classes for D_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         IF (axis==0) call errore('divide_class','unknown D_4h axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=0
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=3
         ELSE
            DO ipol=1,3
               IF (is_axis(ax,ipol)) which_irr(iclass)=4
            ENDDO
            IF (which_irr(iclass)==0) which_irr(iclass)=5
         END IF
      ELSEIF (ts==2) THEN
         which_irr(iclass)=6
      ELSEIF (ts==5) THEN
         which_irr(iclass)=0
         CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=8
         ELSE 
            DO ipol=1,3
               IF (is_axis(ax,ipol)) which_irr(iclass)=9
            ENDDO
            IF (which_irr(iclass)==0) which_irr(iclass)=10
         END IF
      ELSEIF (ts==6) THEN
         which_irr(iclass)=7
      END IF
   END DO
ELSEIF (code_group==23) THEN
!
!  D_6h
!
   IF (nclass /= 12) CALL errore('divide_class','Wrong classes for D_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==1) THEN
            which_irr(iclass)=4
         ELSE
            CALL versor(smat(1,1,elem(1,iclass)),ax)
            CALL which_c2(ax, iax)
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=5
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=6
            ELSE
               CALL errore('divide_class','Problem with C_2 of D_6h',1)
            ENDIF
         END IF
      ELSE IF (ts==2) THEN
          which_irr(iclass)=7
      ELSE IF (ts==5) THEN
          IF (nelem(iclass)==1) THEN
             which_irr(iclass)=10
          ELSE 
             CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
             CALL which_c2(ax, iax)
             IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
                which_irr(iclass)=11
             ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
                which_irr(iclass)=12
             ELSE
                CALL errore('divide_class','Problem with mirror of D_6h',1)
             ENDIF
          END IF
      ELSE IF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
             which_irr(iclass)=9
         ELSE
             which_irr(iclass)=8
         END IF
      END IF
   END DO
ELSEIF (code_group==24) THEN
!
!  D_2d
!
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for D_2d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==6) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==1) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=5
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==25) THEN
!
!  D_3d
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_3d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==2) THEN
         which_irr(iclass)=4
      ELSE IF (ts==6) THEN
         which_irr(iclass)=5
      ELSE IF (ts==5) THEN
         which_irr(iclass)=6
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==26) THEN
!
!  S_4
!
    IF (nclass /= 4) CALL errore('divide_class','Wrong classes for S_4',1)
    DO iclass=2,nclass
       ts=tipo_sym(smat(1,1,elem(1,iclass)))
       IF (ts==4) THEN
          which_irr(iclass)=3
       ELSE IF (ts==6) THEN
          IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-90.d0)<eps) THEN
             which_irr(iclass)=4
          ELSE
             which_irr(iclass)=2
          END IF
       ELSE
          CALL errore('divide_class','wrong operation',1)
       END IF
    END DO
ELSE IF (code_group==27) THEN
!
!  S_6
!
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for S_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=4
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=6
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==28) THEN
!
!  T
!
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for T',1)
   DO iclass=2,nclass
      ars=angle_rot(smat(1,1,elem(1,iclass)))
      IF (ABS(ars-120.d0)<eps) THEN
         which_irr(iclass)=3
      ELSE IF (ABS(ars-240.d0)<eps) THEN
         which_irr(iclass)=4
      ELSE IF (ABS(ars-180.d0)<eps) THEN
         which_irr(iclass)=2
      ELSE
         CALL errore('divide_class','wrong angle',1)
      END IF
   END DO
ELSE IF (code_group==29) THEN
!
!  T_h
!
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for T_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==4) THEN
         which_irr(iclass)=2
      ELSE IF (ts==2) THEN
         which_irr(iclass)=5
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=8
         ELSE
            which_irr(iclass)=7
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=6
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==30) THEN
!
!  T_d
!
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for T_d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==6) THEN
         which_irr(iclass)=4
      ELSE IF (ts==5) THEN
         which_irr(iclass)=5
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==31) THEN
!
!  O
!
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for O',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         IF (nelem(iclass)==3) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==32) THEN
!
!  O_h
!
   IF (nclass /= 10) CALL errore('divide_class','Wrong classes for O_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         IF (nelem(iclass)==3) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=6
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)==6) THEN
            which_irr(iclass)=10
         ELSE
            which_irr(iclass)=8
         END IF
      ELSE IF (ts==6) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=7
         ELSE
            which_irr(iclass)=9
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   ENDDO
ELSE
 CALL errore('divide_class','code_group not correct',1)
ENDIF

RETURN
END SUBROUTINE divide_class

!-----------------------------------------------------------------------------
SUBROUTINE coniug_mat(a,b,c)
!-----------------------------------------------------------------------------
USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP) :: a(3,3), b(3,3), c(3,3)

c=MATMUL(a,MATMUL(b,TRANSPOSE(a)))

RETURN
END SUBROUTINE coniug_mat

!-----------------------------------------------------------------------------
FUNCTION compare_mat(a,b)
!-----------------------------------------------------------------------------
!
!  This function compares two 3x3 matrix and returns .true. if they coincide.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: a(3,3), b(3,3)
REAL(DP), PARAMETER :: eps=1.d-7
LOGICAL :: compare_mat

compare_mat=(ABS(MAXVAL(a-b))<eps).AND.(ABS(MINVAL(a-b))<eps)

END FUNCTION compare_mat

FUNCTION is_axis(ax,iflag)
!
!   This is a logical function which returns .true. if ax is the versor
!   of the axis x (iflag=1) or y (iflag=2) or z (iflag=3)
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: ax(3)
INTEGER :: iflag

LOGICAL :: is_axis

REAL(DP), PARAMETER :: eps=1.d-7

IF (iflag==1) THEN
   is_axis=ABS(ax(2))<eps.AND.ABS(ax(3))<eps
ELSE IF (iflag==2) THEN
   is_axis=ABS(ax(1))<eps.AND.ABS(ax(3))<eps
ELSE IF (iflag==3) THEN
   is_axis=ABS(ax(1))<eps.AND.ABS(ax(2))<eps
ELSE
   call errore('is_axis','iflag not allowed',1)
END IF

RETURN
END FUNCTION is_axis

!-----------------------------------------------------------------------------
SUBROUTINE versor(smat,ax)
!-----------------------------------------------------------------------------
!
!  This subroutine receives a rotation matrix and determines the rotation
!  axis. The orientation of the axis is with the tip in the hemisphere
!  z>=0. In the xy plane the axis is in the x>0 region and the positive
!  y axis is taken for z=0 and x=0.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: smat(3,3), ax(3)
REAL(DP), PARAMETER  ::  eps=1.d-7
REAL(DP) :: a1(3), norm
INTEGER :: ipol, jpol, tipo_sym, ts
!
!  Check if it is a 180 rotation
!
ts=tipo_sym(smat)
IF (ts/=3.and.ts/=4.and.ts/=6) &
     call errore('versor','called in the wrong case',1)
IF (ts==4) THEN
!
!   First the case where the axis is parallel to a coordinate axis
!
   ax=0.d0
   DO ipol=1,3
      IF (ABS(smat(ipol,ipol)-1.d0) < eps ) ax(ipol)=1.d0
   END DO
   norm=sqrt(ax(1)**2+ax(2)**2+ax(3)**2)
   IF (ABS(norm)>eps) RETURN
!
!   then the general case
!
   DO ipol=1,3
      a1(ipol)=sqrt(ABS(smat(ipol,ipol)+1.d0)/2.d0)
   END DO
   
   DO ipol=1,3
      DO jpol=ipol+1,3
         IF (ABS(a1(ipol)*a1(jpol))>eps) THEN
            a1(ipol)=0.5d0*smat(ipol,jpol)/a1(jpol)
         END IF
      END DO
   END DO

ELSE
!
!  It is not a 180 rotation: compute the rotation axis
!
   a1(1) =-smat(2,3)+smat(3,2)
   a1(2) =-smat(3,1)+smat(1,3)
   a1(3) =-smat(1,2)+smat(2,1)

END IF
!
!  The direction of the axis is arbitrarily chosen, with positive z. In the
!  xy plane with positive x, and along y with positive y.
!
IF (a1(3) < -eps ) THEN
   a1=-a1
ELSEIF (abs(a1(3))<eps .and. a1(1) < -eps ) THEN
   a1=-a1
ELSEIF (abs(a1(3))<eps .and. abs(a1(1))<eps.and.a1(2) < -eps ) THEN
   a1=-a1
ENDIF

norm=sqrt(a1(1)**2+a1(2)**2+a1(3)**2)
IF (norm<eps) call errore('versor','problem with the matrix',1)
ax=a1/norm

RETURN
END SUBROUTINE versor

!-----------------------------------------------------------------------------
SUBROUTINE mirror_axis(smat,ax)
!-----------------------------------------------------------------------------
!
!  This subroutine receives a mirror symmetry and determine the normal
!  to the mirror plane. The sign of the normal is undefined.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: smat(3,3), ax(3)
REAL(DP) :: aux_mat(3,3)

aux_mat=-smat

CALL versor(aux_mat,ax)

RETURN

END SUBROUTINE mirror_axis

!-----------------------------------------------------------------------------
FUNCTION angle_rot(smat)
!-----------------------------------------------------------------------------
!
!  This subroutine receives a rotation matrix and determine the 
!  rotation angle. 
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), PARAMETER  ::  eps=1.d-7

REAL(DP) :: smat(3,3)

REAL(DP) :: a1(3), ax(3)
REAL(DP) :: angle_rot, angle_rot1, pi, sint, cost
 
INTEGER :: tipo_sym
!
!  Check if it is a 180 rotation
!
IF (tipo_sym(smat)==4) THEN
   angle_rot=180.d0
   RETURN
END IF
pi=4.d0*atan(1.d0)
!
!  Compute the axis
!
a1(1) =-smat(2,3)+smat(3,2)
a1(2) =-smat(3,1)+smat(1,3)
a1(3) =-smat(1,2)+smat(2,1)

sint=0.5d0*sqrt(a1(1)**2+a1(2)**2+a1(3)**2)
IF (sint<eps) CALL errore('angle_rot','problem with the matrix',1)
IF (ABS(sint)> 1.0_DP+eps) CALL errore('angle_rot','problem with sint',1)
!
!  small rounding errors that make |sint|>1.0 produce NaN in the next ASIN
!  function, so we remove them
!
IF (ABS(sint) > 1.0_DP) sint=SIGN(1.0_DP, sint)
!
!  The direction of the axis is chosen in such a way that a1(3) is always
!  positive if non zero. Otherwise a1(2) is positive, or a1(1) respectively
!
ax=a1
IF (ax(3) < -eps ) THEN
   ax=-ax
ELSEIF (abs(ax(3))<eps .and. ax(2) < -eps ) THEN
   ax=-ax
ELSEIF (abs(ax(3))<eps .and. abs(ax(2))<eps.and.ax(1) < -eps ) THEN
   ax=-ax
ENDIF
IF (ABS(a1(1))>eps) THEN
   sint=SIGN(sint,a1(1)/ax(1)) 
ELSEIF (ABS(a1(2))>eps) THEN
   sint=SIGN(sint,a1(2)/ax(2)) 
ELSEIF (ABS(a1(3))>eps) THEN
   sint=SIGN(sint,a1(3)/ax(3)) 
END IF
!
!  Compute the cos of the angle
!
ax=a1/(2.d0*sint)
IF (ABS(ax(1)**2-1.d0)>eps) THEN
   cost=(smat(1,1)-ax(1)**2)/(1.d0-ax(1)**2)
ELSE IF (ABS(ax(2)**2-1.d0)>eps) THEN
   cost=(smat(2,2)-ax(2)**2)/(1.d0-ax(2)**2)
ELSE IF (ABS(ax(3)**2-1.d0)>eps) THEN
   cost=(smat(3,3)-ax(3)**2)/(1.d0-ax(3)**2)
END IF

IF (ABS(sint**2+cost**2-1.d0) > eps ) &
       CALL errore('angle_rot','problem with the matrix',1)
angle_rot1=ASIN(sint)*180.d0/pi
IF (angle_rot1 < 0.d0) THEN
   IF (cost < 0.d0) THEN
      angle_rot1=-angle_rot1+180.d0
   ELSE
      angle_rot1=360.d0+angle_rot1
   ENDIF
ELSE
   IF (cost < 0.d0) angle_rot1=-angle_rot1+180.d0
ENDIF

angle_rot=angle_rot1

RETURN
END FUNCTION angle_rot

!-----------------------------------------------------------------------------
FUNCTION angle_rot_s(smat)
!-----------------------------------------------------------------------------
!
!  This subroutine receives an improper rotation matrix and determines the 
!  rotation angle. 
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP) :: smat(3,3)

REAL(DP) :: aux_mat(3,3)
REAL(DP) :: angle_rot, angle_rot_s

aux_mat=-smat
angle_rot_s=mod(angle_rot(aux_mat)+180.0_DP,360.0_DP)

RETURN

END FUNCTION angle_rot_s


!-----------------------------------------------------------------------------
SUBROUTINE set_irr_rap(code_group,nclass_ref,char_mat,name_rap, &
                       name_class,ir_ram)
!-----------------------------------------------------------------------------
!
!  This subroutine collects the character tables of the 32 crystallographic
!  point groups.
!  Various names have been used in the litterature to identify
!  the irreducible representations. Several equivalent names are
!  collected in this routine. The first name is taken 
!  from the book of P.W. Atkins, M.S. Child, and C.S.G. Phillips, 
!  "Tables for group theory".
!   D, G, L, S are used for Delta, Gamma, Lambda and Sigma.
!   Representations which correspond to infrared or raman active modes
!   are identified with the string in ir_ram: I (infrared active), 
!   R (Raman active), I+R (Infrared and Raman active).
!   
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: nclass_ref, &    ! Output: number of irreducible representation
           code_group       ! Input: code of the group 

CHARACTER(LEN=15) :: name_rap(12)   ! Output: name of the representations
CHARACTER(LEN=5) :: name_class(12) ! Output: name of the classes
CHARACTER(LEN=3) :: ir_ram(12)

COMPLEX(DP) :: char_mat(12,12) ! Output: character matrix

REAL(DP) :: sqr3d2

sqr3d2=SQRT(3.d0)*0.5d0

char_mat=(1.d0,0.d0)

name_class(1)="E   "

ir_ram="   "

IF (code_group==1) THEN
!
! C_1
!
   nclass_ref=1

   name_rap(1)="A    "
   ir_ram(1)="I+R"

ELSEIF (code_group==2) THEN
!
! C_i
!
   nclass_ref=2
   name_class(2)="i    "

   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="A_u "
   ir_ram(2)="I"
   char_mat(2,2)=(-1.d0,0.d0)

ELSEIF (code_group==3) THEN
!
! C_s
!
   nclass_ref=2
   name_class(2)="s    "

   name_rap(1)="A'  "
   ir_ram(1)="I+R"

   name_rap(2)="A'' "
   ir_ram(2)="I+R"
   char_mat(2,2)=(-1.d0,0.d0)

ELSEIF (code_group==4) THEN
!
! C_2
!
   nclass_ref=2
   name_class(2)="C2   "

   name_rap(1)="A   "
   ir_ram(1)="I+R"

   name_rap(2)="B   "
   ir_ram(2)="I+R"
   char_mat(2,2)=(-1.d0,0.d0)

ELSEIF (code_group==5) THEN
!
! C_3
!
   nclass_ref=3
   name_class(2)="C3   "
   name_class(3)="C3^2 "

   name_rap(1)="A   "
   ir_ram(1)="I+R"

   name_rap(2)="E   "
   ir_ram(2)="I+R"
   char_mat(2,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(2,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   
   name_rap(3)="E*  "
   ir_ram(3)="I+R"
   char_mat(3,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)


ELSEIF (code_group==6) THEN
!
! C_4
!
   nclass_ref=4
   name_class(2)="C4   "
   name_class(3)="C2   "
   name_class(4)="C4^3 "

   name_rap(1)="A   "
   ir_ram(1)="I+R"

   name_rap(2)="B   "
   ir_ram(2)="R"
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)

   name_rap(3)="E   "
   ir_ram(3)="I+R"
   char_mat(3,2)=( 0.d0,1.d0)
   char_mat(3,3)=(-1.d0,0.d0)
   char_mat(3,4)=( 0.d0,-1.d0)

   name_rap(4)="E*  "
   ir_ram(4)="I+R"
   char_mat(4,2)=( 0.d0,-1.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,4)=( 0.d0,1.d0)


ELSEIF (code_group==7) THEN
!
! C_6
!
   nclass_ref=6
   name_class(2)="C6   "
   name_class(3)="C3   "
   name_class(4)="C2   "
   name_class(5)="C3^2 "
   name_class(6)="C6^5 "

   name_rap(1)="A   "
   ir_ram(1)="I+R"

   name_rap(2)="B   "
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)

   name_rap(3)="E_1 "
   ir_ram(3)="I+R"
   char_mat(3,2)=CMPLX( 0.5d0,sqr3d2,kind=DP)
   char_mat(3,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

   name_rap(4)="E_1*"
   ir_ram(4)="I+R"
   char_mat(4,2)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(4,6)=CMPLX( 0.5d0,sqr3d2,kind=DP)

   name_rap(5)="E_2 "
   ir_ram(5)="R"
   char_mat(5,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(5,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,5)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(5,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(6)="E_2*"
   ir_ram(6)="R"
   char_mat(6,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(6,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,6)=CMPLX(-0.5d0,sqr3d2,kind=DP)

ELSEIF (code_group==8) THEN
!
! D_2
!
   nclass_ref=4
   name_class(2)="C2   "
   name_class(3)="C2'  "
   name_class(4)="C2'' "

   name_rap(1)="A   "
   ir_ram(1)="R"

   name_rap(2)="B_1 "
   ir_ram(2)="I+R"
   char_mat(2,3)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)

   name_rap(3)="B_2 "
   ir_ram(3)="I+R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)

   name_rap(4)="B_3 "
   ir_ram(4)="I+R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)

ELSEIF (code_group==9) THEN
!
! D_3
!
   nclass_ref=3
   name_class(2)="2C3  "
   name_class(3)="3C2' "

   name_rap(1)="A_1 "
   ir_ram(1)="R"

   name_rap(2)="A_2 "
   ir_ram(2)="I"
   char_mat(2,3)=(-1.d0,0.d0)

   name_rap(3)="E   "
   ir_ram(3)="I+R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 0.d0,0.d0)

ELSEIF (code_group==10) THEN
!
! D_4
!
   nclass_ref=5
   name_class(2)="2C4  "
   name_class(3)="C2   "
   name_class(4)="2C2' "
   name_class(5)="2C2''"

   name_rap(1)="A_1 "
   ir_ram(1)="R"

   name_rap(2)="A_2 "
   ir_ram(2)="I"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)

   name_rap(3)="B_1 "
   ir_ram(3)="R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)

   name_rap(4)="B_2 "
   ir_ram(4)="R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)

   name_rap(5)="E   "
   ir_ram(5)="I+R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-2.d0,0.d0)
   char_mat(5,4)=( 0.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)

ELSEIF (code_group==11) THEN
!
! D_6
!
   nclass_ref=6
   name_class(2)="2C6  "
   name_class(3)="2C3  "
   name_class(4)="C2   "
   name_class(5)="3C2' "
   name_class(6)="3C2''"

   name_rap(1)="A_1 "
   ir_ram(1)="R"

   name_rap(2)="A_2 "
   ir_ram(2)="I"
   char_mat(2,5)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)

   name_rap(3)="B_1 "
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,6)=(-1.d0,0.d0)

   name_rap(4)="B_2 "
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)

   name_rap(5)="E_1 "
   ir_ram(5)="I+R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-2.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)
   char_mat(5,6)=( 0.d0,0.d0)

   name_rap(6)="E_2 "
   ir_ram(6)="R"
   char_mat(6,1)=( 2.d0,0.d0)
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,3)=(-1.d0,0.d0)
   char_mat(6,4)=( 2.d0,0.d0)
   char_mat(6,5)=( 0.d0,0.d0)
   char_mat(6,6)=( 0.d0,0.d0)

ELSEIF (code_group==12) THEN
!
! C_2v
!
   nclass_ref=4
   name_class(2)="C2   "
   name_class(3)="s_v  "
   name_class(4)="s_v' "

   name_rap(1)="A_1  D_1  S_1"
   ir_ram(1)="I+R"
 
   name_rap(2)="A_2  D_2  S_2"
   ir_ram(2)="R"
   char_mat(2,3)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)

   name_rap(3)="B_1  D_3  S_3"
   ir_ram(3)="I+R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)

   name_rap(4)="B_2  D_4  S_4"
   ir_ram(4)="I+R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)

ELSEIF (code_group==13) THEN
!
! C_3v
!
   nclass_ref=3
   name_class(2)="2C3  "
   name_class(3)="3s_v "

   name_rap(1)="A_1  L_1"
   ir_ram(1)="I+R"

   name_rap(2)="A_2  L_2"
   char_mat(2,3)=(-1.d0,0.d0)

   name_rap(3)="E    L_3"
   ir_ram(3)="I+R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 0.d0,0.d0)

ELSEIF (code_group==14) THEN
!
! C_4v
!
   nclass_ref=5
   name_class(2)="2C4  "
   name_class(3)="C2   "
   name_class(4)="2s_v "
   name_class(5)="2s_d "

   name_rap(1)="A_1  G_1 D_1"
   ir_ram(1)="I+R"

   name_rap(2)="A_2  G_2 D_1'"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)

   name_rap(3)="B_1  G_3 D_2"
   ir_ram(3)="R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)

   name_rap(4)="B_2  G_4 D_2'"
   ir_ram(4)="R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)

   name_rap(5)="E    G_5 D_5"
   ir_ram(5)="I+R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-2.d0,0.d0)
   char_mat(5,4)=( 0.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)

ELSEIF (code_group==15) THEN
!
! C_6v
!
   nclass_ref=6
   name_class(2)="2C6  "
   name_class(3)="2C3  "
   name_class(4)="C2   "
   name_class(5)="3s_v "
   name_class(6)="3s_d "

   name_rap(1)="A_1 "
   ir_ram(1)="I+R"

   name_rap(2)="A_2 "
   char_mat(2,5)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)

   name_rap(3)="B_1 "
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,6)=(-1.d0,0.d0)

   name_rap(4)="B_2 "
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)

   name_rap(5)="E_1 "
   ir_ram(5)="I+R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-2.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)
   char_mat(5,6)=( 0.d0,0.d0)

   name_rap(6)="E_2 "
   ir_ram(6)="R"
   char_mat(6,1)=( 2.d0,0.d0)
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,3)=(-1.d0,0.d0)
   char_mat(6,4)=( 2.d0,0.d0)
   char_mat(6,5)=( 0.d0,0.d0)
   char_mat(6,6)=( 0.d0,0.d0)

ELSEIF (code_group==16) THEN
!
! C_2h
!
   nclass_ref=4
   name_class(2)="C2   "
   name_class(3)="i    "
   name_class(4)="s_h  "

   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="B_g "
   ir_ram(2)="R"
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)

   name_rap(3)="A_u "
   ir_ram(3)="I"
   char_mat(3,3)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)

   name_rap(4)="B_u "
   ir_ram(4)="I"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)

ELSEIF (code_group==17) THEN
!
! C_3h
!
   nclass_ref=6
   name_class(2)="C3   "
   name_class(3)="C3^2 "
   name_class(4)="s_h  "
   name_class(5)="S3   "
   name_class(6)="S3^5 "

   name_rap(1)="A'  "
   ir_ram(1)="R"
 
   name_rap(2)="E'  "
   ir_ram(2)="I+R"
   char_mat(2,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(2,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(3)="E'* "
   ir_ram(3)="I+R"
   char_mat(3,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(3,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0,sqr3d2,kind=DP)

   name_rap(4)="A'' "
   ir_ram(4)="I"
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,6)=(-1.d0,0.d0)

   name_rap(5)="E'' "
   ir_ram(5)="R"
   char_mat(5,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(5,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,4)=(-1.d0,0.d0)
   char_mat(5,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(5,6)=CMPLX(0.5d0,sqr3d2,kind=DP)

   name_rap(6)="E''*"
   ir_ram(6)="R"
   char_mat(6,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(6,4)=(-1.d0,0.d0)
   char_mat(6,5)=CMPLX( 0.5d0,sqr3d2,kind=DP)
   char_mat(6,6)=CMPLX(0.5d0,-sqr3d2,kind=DP)


ELSEIF (code_group==18) THEN
!
! C_4h
!
   nclass_ref=8
   name_class(2)="C4   "
   name_class(3)="C2   "
   name_class(4)="C4^3 "
   name_class(5)="i    "
   name_class(6)="S4^3 "
   name_class(7)="s_h  "
   name_class(8)="S4   "

   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="B_g "
   ir_ram(2)="R"
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)
   char_mat(2,8)=(-1.d0,0.d0)

   name_rap(3)="E_g "
   ir_ram(3)="R"
   char_mat(3,2)=( 0.d0,1.d0)
   char_mat(3,3)=(-1.d0,0.d0)
   char_mat(3,4)=( 0.d0,-1.d0)
   char_mat(3,6)=( 0.d0,1.d0)
   char_mat(3,7)=(-1.d0,0.d0)
   char_mat(3,8)=( 0.d0,-1.d0)

   name_rap(4)="E_g*"
   ir_ram(4)="R"
   char_mat(4,2)=(0.d0,-1.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,4)=( 0.d0,1.d0)
   char_mat(4,6)=( 0.d0,-1.d0)
   char_mat(4,7)=(-1.d0,0.d0)
   char_mat(4,8)=( 0.d0,1.d0)

   name_rap(5)="A_u "
   ir_ram(5)="I"
   char_mat(5,5)=(-1.d0,0.d0)
   char_mat(5,6)=(-1.d0,0.d0)
   char_mat(5,7)=(-1.d0,0.d0)
   char_mat(5,8)=(-1.d0,0.d0)

   name_rap(6)="B_u "
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,4)=(-1.d0,0.d0)
   char_mat(6,5)=(-1.d0,0.d0)
   char_mat(6,7)=(-1.d0,0.d0)

   name_rap(7)="E_u "
   ir_ram(7)="I"
   char_mat(7,2)=( 0.d0,1.d0)
   char_mat(7,3)=(-1.d0,0.d0)
   char_mat(7,4)=( 0.d0,-1.d0)
   char_mat(7,5)=(-1.d0, 0.d0)
   char_mat(7,6)=( 0.d0,-1.d0)
   char_mat(7,8)=( 0.d0,1.d0)

   name_rap(8)="E_u*"
   ir_ram(8)="I"
   char_mat(8,2)=( 0.d0,-1.d0)
   char_mat(8,3)=(-1.d0,0.d0)
   char_mat(8,4)=( 0.d0,1.d0)
   char_mat(8,5)=(-1.d0, 0.d0)
   char_mat(8,6)=( 0.d0,1.d0)
   char_mat(8,8)=( 0.d0,-1.d0)

ELSEIF (code_group==19) THEN
!
! C_6h
!
   nclass_ref=12
   name_class(2)="C6   "
   name_class(3)="C3   "
   name_class(4)="C2   "
   name_class(5)="C3^2 "
   name_class(6)="C6^5 "
   name_class(7)="i    "
   name_class(8)="S3^5 "
   name_class(9)="S6^5 "
   name_class(10)="s_h "
   name_class(11)="S6  "
   name_class(12)="S3  "

   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="B_g "
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)
   char_mat(2,8)=(-1.d0,0.d0)
   char_mat(2,10)=(-1.d0,0.d0)
   char_mat(2,12)=(-1.d0,0.d0)

   name_rap(3)="E_1g"
   ir_ram(3)="R"
   char_mat(3,2)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(3,8)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(3,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,10)=(-1.d0,0.d0)
   char_mat(3,11)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,12)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

   name_rap(4)="E1g*"
   ir_ram(4)="R"
   char_mat(4,2)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(4,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(4,8)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(4,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(4,10)=(-1.d0,0.d0)
   char_mat(4,11)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(4,12)=CMPLX( 0.5d0,sqr3d2,kind=DP)

   name_rap(5)="E_2g"
   ir_ram(5)="R"
   char_mat(5,2)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,8)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,9)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,11)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(5,12)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(6)="E2g*"
   ir_ram(6)="R"
   char_mat(6,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(6,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(6,8)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,9)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(6,11)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,12)=CMPLX(-0.5d0, sqr3d2,kind=DP)

   name_rap(7)="A_u "
   ir_ram(7)="I"
   char_mat(7,7)=(-1.d0,0.d0)
   char_mat(7,8)=(-1.d0,0.d0)
   char_mat(7,9)=(-1.d0,0.d0)
   char_mat(7,10)=(-1.d0,0.d0)
   char_mat(7,11)=(-1.d0,0.d0)
   char_mat(7,12)=(-1.d0,0.d0)

   name_rap(8)="B_u "
   char_mat(8,2)=(-1.d0,0.d0)
   char_mat(8,4)=(-1.d0,0.d0)
   char_mat(8,6)=(-1.d0,0.d0)
   char_mat(8,7)=(-1.d0,0.d0)
   char_mat(8,9)=(-1.d0,0.d0)
   char_mat(8,11)=(-1.d0,0.d0)

   name_rap(9)="E_1u"
   ir_ram(9)="I"
   char_mat(9,2)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(9,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(9,4)=(-1.d0,0.d0)
   char_mat(9,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(9,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(9,7)=(-1.d0,0.d0)
   char_mat(9,8)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(9,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(9,11)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(9,12)=CMPLX(-0.5d0, sqr3d2,kind=DP)

   name_rap(10)="E1u*"
   ir_ram(10)="I"
   char_mat(10,2)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(10,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(10,4)=(-1.d0,0.d0)
   char_mat(10,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(10,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(10,7)=(-1.d0,0.d0)
   char_mat(10,8)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(10,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(10,11)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(10,12)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(11)="E_2u"
   char_mat(11,2)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(11,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(11,5)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(11,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(11,7)=(-1.d0,0.d0)
   char_mat(11,8)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(11,9)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(11,10)=(-1.d0,0.d0)
   char_mat(11,11)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(11,12)=CMPLX( 0.5d0, sqr3d2,kind=DP)

   name_rap(12)="E2u*"
   char_mat(12,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(12,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(12,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(12,6)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(12,7)=(-1.d0,0.d0)
   char_mat(12,8)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(12,9)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(12,10)=(-1.d0,0.d0)
   char_mat(12,11)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(12,12)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

ELSEIF (code_group==20) THEN
!
! D_2h
!
   nclass_ref=8
   name_class(2)="C2   "
   name_class(3)="C2'  "
   name_class(4)="C2'' "
   name_class(5)="i    "
   name_class(6)="s_v  "
   name_class(7)="s_v' "
   name_class(8)="s_v''"

   name_rap(1)="A_g "
   ir_ram(1)="R"
 
   name_rap(2)="B_1g"
   ir_ram(2)="R"
   char_mat(2,3)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,7)=(-1.d0,0.d0)
   char_mat(2,8)=(-1.d0,0.d0)

   name_rap(3)="B_2g"
   ir_ram(3)="R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,6)=(-1.d0,0.d0)
   char_mat(3,8)=(-1.d0,0.d0)

   name_rap(4)="B_3g"
   ir_ram(4)="R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,6)=(-1.d0,0.d0)
   char_mat(4,7)=(-1.d0,0.d0)

   name_rap(5)="A_u "
   char_mat(5,5)=(-1.d0,0.d0)
   char_mat(5,6)=(-1.d0,0.d0)
   char_mat(5,7)=(-1.d0,0.d0)
   char_mat(5,8)=(-1.d0,0.d0)

   name_rap(6)="B_1u"
   ir_ram(6)="I"
   char_mat(6,3)=(-1.d0,0.d0)
   char_mat(6,4)=(-1.d0,0.d0)
   char_mat(6,5)=(-1.d0,0.d0)
   char_mat(6,6)=(-1.d0,0.d0)

   name_rap(7)="B_2u"
   ir_ram(7)="I"
   char_mat(7,2)=(-1.d0,0.d0)
   char_mat(7,4)=(-1.d0,0.d0)
   char_mat(7,5)=(-1.d0,0.d0)
   char_mat(7,7)=(-1.d0,0.d0)

   name_rap(8)="B_3u"
   ir_ram(8)="I"
   char_mat(8,2)=(-1.d0,0.d0)
   char_mat(8,3)=(-1.d0,0.d0)
   char_mat(8,5)=(-1.d0,0.d0)
   char_mat(8,8)=(-1.d0,0.d0)

ELSEIF (code_group==21) THEN
!
! D_3h
!
   nclass_ref=6
   name_class(2)="2C3  "
   name_class(3)="3C2  "
   name_class(4)="s_h  "
   name_class(5)="2S3  "
   name_class(6)="3s_v "

   name_rap(1)="A'_1"
   ir_ram(1)="R"

   name_rap(2)="A'_2"
   char_mat(2,3)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)

   name_rap(3)="E'  "
   ir_ram(3)="I+R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 0.d0,0.d0)
   char_mat(3,4)=( 2.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)
   char_mat(3,6)=( 0.d0,0.d0)

   name_rap(4)="A''1"
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,6)=(-1.d0,0.d0)

   name_rap(5)="A''2"
   ir_ram(5)="I"
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-1.d0,0.d0)
   char_mat(5,5)=(-1.d0,0.d0)

   name_rap(6)="E'' "
   ir_ram(6)="R"
   char_mat(6,1)=( 2.d0,0.d0)
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,3)=( 0.d0,0.d0)
   char_mat(6,4)=(-2.d0,0.d0)
   char_mat(6,6)=( 0.d0,0.d0)

ELSEIF (code_group==22) THEN
!
! D_4h
!
   nclass_ref=10
   name_class(2)="2C4  "
   name_class(3)="C2   "
   name_class(4)="2C2' "
   name_class(5)="2C2''"
   name_class(6)="i    "
   name_class(7)="2S4  "
   name_class(8)="s_h  "
   name_class(9)="2s_v "
   name_class(10)="2s_d "


   name_rap(1)="A_1g X_1  M_1"
   ir_ram(1)="R"

   name_rap(2)="A_2g X_4  M_4"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)
   char_mat(2,9)=(-1.d0,0.d0)
   char_mat(2,10)=(-1.d0,0.d0)

   name_rap(3)="B_1g X_2  M_2"
   ir_ram(3)="R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)
   char_mat(3,7)=(-1.d0,0.d0)
   char_mat(3,10)=(-1.d0,0.d0)

   name_rap(4)="B_2g X_3  M_3"
   ir_ram(4)="R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,7)=(-1.d0,0.d0)
   char_mat(4,9)=(-1.d0,0.d0)

   name_rap(5)="E_g  X_5  M_5"
   ir_ram(5)="R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-2.d0,0.d0)
   char_mat(5,4)=( 0.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)
   char_mat(5,6)=( 2.d0,0.d0)
   char_mat(5,7)=( 0.d0,0.d0)
   char_mat(5,8)=(-2.d0,0.d0)
   char_mat(5,9)=( 0.d0,0.d0)
   char_mat(5,10)=( 0.d0,0.d0)

   name_rap(6)="A_1u X_1' M_1'"
   char_mat(6,6)=(-1.d0,0.d0)
   char_mat(6,7)=(-1.d0,0.d0)
   char_mat(6,8)=(-1.d0,0.d0)
   char_mat(6,9)=(-1.d0,0.d0)
   char_mat(6,10)=(-1.d0,0.d0)

   name_rap(7)="A_2u X_4' M_4'"
   ir_ram(7)="I"
   char_mat(7,4)=(-1.d0,0.d0)
   char_mat(7,5)=(-1.d0,0.d0)
   char_mat(7,6)=(-1.d0,0.d0)
   char_mat(7,7)=(-1.d0,0.d0)
   char_mat(7,8)=(-1.d0,0.d0)

   name_rap(8)="B_1u X_2' M_2'"
   char_mat(8,2)=(-1.d0,0.d0)
   char_mat(8,5)=(-1.d0,0.d0)
   char_mat(8,6)=(-1.d0,0.d0)
   char_mat(8,8)=(-1.d0,0.d0)
   char_mat(8,9)=(-1.d0,0.d0)

   name_rap(9)="B_2u X_3' M_3'"
   char_mat(9,2)=(-1.d0,0.d0)
   char_mat(9,4)=(-1.d0,0.d0)
   char_mat(9,6)=(-1.d0,0.d0)
   char_mat(9,8)=(-1.d0,0.d0)
   char_mat(9,10)=(-1.d0,0.d0)

   name_rap(10)="E_u  X_5' M_5'"
   ir_ram(10)="I"
   char_mat(10,1)=( 2.d0,0.d0)
   char_mat(10,2)=( 0.d0,0.d0)
   char_mat(10,3)=(-2.d0,0.d0)
   char_mat(10,4)=( 0.d0,0.d0)
   char_mat(10,5)=( 0.d0,0.d0)
   char_mat(10,6)=(-2.d0,0.d0)
   char_mat(10,7)=( 0.d0,0.d0)
   char_mat(10,8)=( 2.d0,0.d0)
   char_mat(10,9)=( 0.d0,0.d0)
   char_mat(10,10)=( 0.d0,0.d0)

ELSEIF (code_group==23) THEN
!
! D_6h
!
   nclass_ref=12
   name_class(2)="2C6  "
   name_class(3)="2C3  "
   name_class(4)="C2   "
   name_class(5)="3C2' "
   name_class(6)="3C2''"
   name_class(7)="i    "
   name_class(8)="2S3  "
   name_class(9)="2S6  "
   name_class(10)="s_h "
   name_class(11)="3s_d "
   name_class(12)="3s_v "

   name_rap(1)="A_1g"
   ir_ram(1)="R"
 
   name_rap(2)="A_2g"
   char_mat(2,5)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)
   char_mat(2,11)=(-1.d0,0.d0)
   char_mat(2,12)=(-1.d0,0.d0)

   name_rap(3)="B_1g"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,4)=(-1.d0,0.d0)
   char_mat(3,6)=(-1.d0,0.d0)
   char_mat(3,8)=(-1.d0,0.d0)
   char_mat(3,10)=(-1.d0,0.d0)
   char_mat(3,12)=(-1.d0,0.d0)

   name_rap(4)="B_2g"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,8)=(-1.d0,0.d0)
   char_mat(4,10)=(-1.d0,0.d0)
   char_mat(4,11)=(-1.d0,0.d0)

   name_rap(5)="E_1g"
   ir_ram(5)="R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-2.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)
   char_mat(5,6)=( 0.d0,0.d0)
   char_mat(5,7)=( 2.d0,0.d0)
   char_mat(5,9)=(-1.d0,0.d0)
   char_mat(5,10)=(-2.d0,0.d0)
   char_mat(5,11)=( 0.d0,0.d0)
   char_mat(5,12)=( 0.d0,0.d0)

   name_rap(6)="E_2g"
   ir_ram(6)="R"
   char_mat(6,1)=( 2.d0,0.d0)
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,3)=(-1.d0,0.d0)
   char_mat(6,4)=( 2.d0,0.d0)
   char_mat(6,5)=( 0.d0,0.d0)
   char_mat(6,6)=( 0.d0,0.d0)
   char_mat(6,7)=( 2.d0,0.d0)
   char_mat(6,8)=(-1.d0,0.d0)
   char_mat(6,9)=(-1.d0,0.d0)
   char_mat(6,10)=( 2.d0,0.d0)
   char_mat(6,11)=( 0.d0,0.d0)
   char_mat(6,12)=( 0.d0,0.d0)

   name_rap(7)="A_1u"
   char_mat(7,7)=(-1.d0,0.d0)
   char_mat(7,8)=(-1.d0,0.d0)
   char_mat(7,9)=(-1.d0,0.d0)
   char_mat(7,10)=(-1.d0,0.d0)
   char_mat(7,11)=(-1.d0,0.d0)
   char_mat(7,12)=(-1.d0,0.d0)

   name_rap(8)="A_2u"
   ir_ram(8)="I"
   char_mat(8,5)=(-1.d0,0.d0)
   char_mat(8,6)=(-1.d0,0.d0)
   char_mat(8,7)=(-1.d0,0.d0)
   char_mat(8,8)=(-1.d0,0.d0)
   char_mat(8,9)=(-1.d0,0.d0)
   char_mat(8,10)=(-1.d0,0.d0)

   name_rap(9)="B_1u"
   char_mat(9,2)=(-1.d0,0.d0)
   char_mat(9,4)=(-1.d0,0.d0)
   char_mat(9,6)=(-1.d0,0.d0)
   char_mat(9,7)=(-1.d0,0.d0)
   char_mat(9,9)=(-1.d0,0.d0)
   char_mat(9,11)=(-1.d0,0.d0)

   name_rap(10)="B_2u"
   char_mat(10,2)=(-1.d0,0.d0)
   char_mat(10,4)=(-1.d0,0.d0)
   char_mat(10,5)=(-1.d0,0.d0)
   char_mat(10,7)=(-1.d0,0.d0)
   char_mat(10,9)=(-1.d0,0.d0)
   char_mat(10,12)=(-1.d0,0.d0)

   name_rap(11)="E_1u"
   ir_ram(11)="I"
   char_mat(11,1)=( 2.d0,0.d0)
   char_mat(11,3)=(-1.d0,0.d0)
   char_mat(11,4)=(-2.d0,0.d0)
   char_mat(11,5)=( 0.d0,0.d0)
   char_mat(11,6)=( 0.d0,0.d0)
   char_mat(11,7)=(-2.d0,0.d0)
   char_mat(11,8)=(-1.d0,0.d0)
   char_mat(11,10)=( 2.d0,0.d0)
   char_mat(11,11)=( 0.d0,0.d0)
   char_mat(11,12)=( 0.d0,0.d0)

   name_rap(12)="E_2u"
   char_mat(12,1)=( 2.d0,0.d0)
   char_mat(12,2)=(-1.d0,0.d0)
   char_mat(12,3)=(-1.d0,0.d0)
   char_mat(12,4)=( 2.d0,0.d0)
   char_mat(12,5)=( 0.d0,0.d0)
   char_mat(12,6)=( 0.d0,0.d0)
   char_mat(12,7)=(-2.d0,0.d0)
   char_mat(12,10)=(-2.d0,0.d0)
   char_mat(12,11)=( 0.d0,0.d0)
   char_mat(12,12)=( 0.d0,0.d0)

ELSEIF (code_group==24) THEN
!
! D_2d
!
   nclass_ref=5
   name_class(2)="2S4  "
   name_class(3)="C2   "
   name_class(4)="2C2' "
   name_class(5)="2s_d "

   name_rap(1)="A_1  X_1  W_1"
   ir_ram(1)="R"

   name_rap(2)="A_2  X_4  W_2'"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)

   name_rap(3)="B_1  X_2  W_1'"
   ir_ram(3)="R"
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)

   name_rap(4)="B_2  X_3  W_2"
   ir_ram(4)="I+R"
   char_mat(4,2)=(-1.d0,0.d0)
   char_mat(4,4)=(-1.d0,0.d0)

   name_rap(5)="E    X_5  W_3"
   ir_ram(5)="I+R"
   char_mat(5,1)=( 2.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-2.d0,0.d0)
   char_mat(5,4)=( 0.d0,0.d0)
   char_mat(5,5)=( 0.d0,0.d0)

ELSEIF (code_group==25) THEN
!
! D_3d
!
   nclass_ref=6
   name_class(2)="2C3  "
   name_class(3)="3C2' "
   name_class(4)="i    "
   name_class(5)="2S6  "
   name_class(6)="3s_d "

   name_rap(1)="A_1g L_1"
   ir_ram(1)="R"

   name_rap(2)="A_2g L_2"
   char_mat(2,3)=(-1.d0,0.d0)
   char_mat(2,6)=(-1.d0,0.d0)

   name_rap(3)="E_g  L_3"
   ir_ram(3)="R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 0.d0,0.d0)
   char_mat(3,4)=( 2.d0,0.d0)
   char_mat(3,5)=(-1.d0,0.d0)
   char_mat(3,6)=( 0.d0,0.d0)

   name_rap(4)="A_1u L_1'"
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,6)=(-1.d0,0.d0)

   name_rap(5)="A_2u L_2'"
   ir_ram(5)="I"
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-1.d0,0.d0)
   char_mat(5,5)=(-1.d0,0.d0)

   name_rap(6)="E_u  L_3'"
   ir_ram(6)="I"
   char_mat(6,1)=( 2.d0,0.d0)
   char_mat(6,2)=(-1.d0,0.d0)
   char_mat(6,3)=( 0.d0,0.d0)
   char_mat(6,4)=(-2.d0,0.d0)
   char_mat(6,6)=( 0.d0,0.d0)

ELSEIF (code_group==26) THEN
!
! S_4
!
   nclass_ref=4
   name_class(2)="S4^3 "
   name_class(3)="C2   "
   name_class(4)="S4   "

   name_rap(1)="A    W_1"
   ir_ram(1)="R"

   name_rap(2)="B    W_3"
   ir_ram(2)="I+R"
   char_mat(2,2)=(-1.d0,0.d0)
   char_mat(2,4)=(-1.d0,0.d0)

   name_rap(3)="E    W_4"
   ir_ram(3)="I+R"
   char_mat(3,2)=( 0.d0, 1.d0)
   char_mat(3,3)=(-1.d0,0.d0)
   char_mat(3,4)=( 0.d0,-1.d0)

   name_rap(4)="E*   W_2"
   ir_ram(4)="I+R"
   char_mat(4,2)=( 0.d0,-1.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,4)=( 0.d0, 1.d0)

ELSEIF (code_group==27) THEN
!
! S_6
!
   nclass_ref=6
   name_class(2)="C3   "
   name_class(3)="C3^2 "
   name_class(4)="i    "
   name_class(5)="S6^5 "
   name_class(6)="S6   "


   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="E_g "
   ir_ram(2)="R"
   char_mat(2,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(2,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,5)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(2,6)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(3)="E_g*"
   ir_ram(3)="R"
   char_mat(3,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(3,5)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,6)=CMPLX(-0.5d0,sqr3d2,kind=DP)

   name_rap(4)="A_u "
   ir_ram(4)="I"
   char_mat(4,4)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,6)=(-1.d0,0.d0)

   name_rap(5)="E_u "
   ir_ram(5)="I"
   char_mat(5,2)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(5,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(5,4)=(-1.d0,0.d0)
   char_mat(5,5)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(5,6)=CMPLX( 0.5d0, sqr3d2,kind=DP)

   name_rap(6)="E_u*"
   ir_ram(6)="I"
   char_mat(6,2)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,3)=CMPLX(-0.5d0,sqr3d2,kind=DP)
   char_mat(6,4)=(-1.d0,0.d0)
   char_mat(6,5)=CMPLX( 0.5d0,sqr3d2,kind=DP)
   char_mat(6,6)=CMPLX( 0.5d0,-sqr3d2,kind=DP)


ELSEIF (code_group==28) THEN
!
! T
!
   nclass_ref=4
   name_class(2)="3C2  "
   name_class(3)="4C3  "
   name_class(4)="4C3' "

   name_rap(1)="A   "
   ir_ram(1)="R"

   name_rap(2)="E   "
   ir_ram(2)="R"
   char_mat(2,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(3)="E*  "
   ir_ram(3)="R"
   char_mat(3,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)

   name_rap(4)="T   "
   ir_ram(4)="I+R"
   char_mat(4,1)=( 3.0d0,0.d0)
   char_mat(4,3)=( 0.0d0,0.d0)
   char_mat(4,4)=( 0.0d0,0.d0)
   char_mat(4,2)=(-1.0d0,0.d0)

ELSEIF (code_group==29) THEN
!
! T_h
!
   nclass_ref=8
   name_class(2)="3C2  "
   name_class(3)="4C3  "
   name_class(4)="4C3' "
   name_class(5)="i    "
   name_class(6)="3s_h "
   name_class(7)="4S6^5"
   name_class(8)="4S6  "

   name_rap(1)="A_g "
   ir_ram(1)="R"

   name_rap(2)="E_g "
   ir_ram(2)="R"
   char_mat(2,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(2,7)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(2,8)=CMPLX(-0.5d0,-sqr3d2,kind=DP)

   name_rap(3)="E_g*"
   ir_ram(3)="R"
   char_mat(3,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(3,7)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(3,8)=CMPLX(-0.5d0, sqr3d2,kind=DP)

   name_rap(4)="T_g "
   ir_ram(4)="R"
   char_mat(4,1)=( 3.0d0,0.d0)
   char_mat(4,2)=(-1.0d0,0.d0)
   char_mat(4,3)=( 0.0d0,0.d0)
   char_mat(4,4)=( 0.0d0,0.d0)
   char_mat(4,5)=( 3.0d0,0.d0)
   char_mat(4,6)=(-1.0d0,0.d0)
   char_mat(4,7)=( 0.0d0,0.d0)
   char_mat(4,8)=( 0.0d0,0.d0)

   name_rap(5)="A_u "
   char_mat(5,5)=(-1.0d0,0.d0)
   char_mat(5,6)=(-1.0d0,0.d0)
   char_mat(5,7)=(-1.0d0,0.d0)
   char_mat(5,8)=(-1.0d0,0.d0)

   name_rap(6)="E_u "
   char_mat(6,3)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(6,4)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(6,6)=(-1.0d0,0.d0)
   char_mat(6,5)=(-1.0d0,0.d0)
   char_mat(6,7)=CMPLX( 0.5d0,-sqr3d2,kind=DP)
   char_mat(6,8)=CMPLX( 0.5d0, sqr3d2,kind=DP)

   name_rap(7)="E_u*"
   char_mat(7,3)=CMPLX(-0.5d0,-sqr3d2,kind=DP)
   char_mat(7,4)=CMPLX(-0.5d0, sqr3d2,kind=DP)
   char_mat(7,5)=(-1.0d0,0.d0)
   char_mat(7,6)=(-1.0d0,0.d0)
   char_mat(7,7)=CMPLX( 0.5d0, sqr3d2,kind=DP)
   char_mat(7,8)=CMPLX( 0.5d0,-sqr3d2,kind=DP)

   name_rap(8)="T_u "
   ir_ram(8)="I"
   char_mat(8,1)=( 3.0d0,0.d0)
   char_mat(8,2)=(-1.0d0,0.d0)
   char_mat(8,3)=( 0.0d0,0.d0)
   char_mat(8,4)=( 0.0d0,0.d0)
   char_mat(8,5)=(-3.0d0,0.d0)
   char_mat(8,7)=( 0.0d0,0.d0)
   char_mat(8,8)=( 0.0d0,0.d0)


ELSEIF (code_group==30) THEN
!
! T_d
!
   nclass_ref=5
   name_class(2)="8C3  "
   name_class(3)="3C2  "
   name_class(4)="6S4  "
   name_class(5)="6s_d "

   name_rap(1)="A_1  G_1  P_1"
   ir_ram(1)="R"

   name_rap(2)="A_2  G_2  P_2"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)

   name_rap(3)="E    G_12 P_3"
   ir_ram(3)="R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 2.d0,0.d0)
   char_mat(3,4)=( 0.d0,0.d0)
   char_mat(3,5)=( 0.d0,0.d0)

   name_rap(4)="T_1  G_25 P_5"
   char_mat(4,1)=( 3.d0,0.d0)
   char_mat(4,2)=( 0.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)

   name_rap(5)="T_2  G_15 P_4"
   ir_ram(5)="I+R"
   char_mat(5,1)=( 3.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-1.d0,0.d0)

ELSEIF (code_group==31) THEN
!
! O
!
   nclass_ref=5
   name_class(2)="8C3  "
   name_class(3)="3C2  "
   name_class(4)="6C4  "
   name_class(5)="6C2' "

   name_rap(1)="A_1 "
   ir_ram(1)="R"

   name_rap(2)="A_2 "
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)

   name_rap(3)="E   "
   ir_ram(3)="R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 2.d0,0.d0)
   char_mat(3,4)=( 0.d0,0.d0)
   char_mat(3,5)=( 0.d0,0.d0)

   name_rap(4)="T_1 "
   ir_ram(4)="I"
   char_mat(4,1)=( 3.d0,0.d0)
   char_mat(4,2)=( 0.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)

   name_rap(5)="T_2 "
   ir_ram(5)="R"
   char_mat(5,1)=( 3.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-1.d0,0.d0)

ELSEIF (code_group==32) THEN
!
! O_h
!
   nclass_ref=10
   name_class(2)="8C3  "
   name_class(3)="3C2  "
   name_class(4)="6C4  "
   name_class(5)="6C2' "
   name_class(6)="i    "
   name_class(7)="8S6  "
   name_class(8)="3s_h "
   name_class(9)="6S4  "
   name_class(10)="6s_d "

   name_rap(1)="A_1g G_1   G_1+"
   ir_ram(1)="R"

   name_rap(2)="A_2g G_2   G_2+"
   char_mat(2,4)=(-1.d0,0.d0)
   char_mat(2,5)=(-1.d0,0.d0)
   char_mat(2,9)=(-1.d0,0.d0)
   char_mat(2,10)=(-1.d0,0.d0)

   name_rap(3)="E_g  G_12  G_3+"
   ir_ram(3)="R"
   char_mat(3,1)=( 2.d0,0.d0)
   char_mat(3,2)=(-1.d0,0.d0)
   char_mat(3,3)=( 2.d0,0.d0)
   char_mat(3,4)=( 0.d0,0.d0)
   char_mat(3,5)=( 0.d0,0.d0)
   char_mat(3,6)=( 2.d0,0.d0)
   char_mat(3,7)=(-1.d0,0.d0)
   char_mat(3,8)=( 2.d0,0.d0)
   char_mat(3,9)=( 0.d0,0.d0)
   char_mat(3,10)=( 0.d0,0.d0)

   name_rap(4)="T_1g G_15' G_4+"
   char_mat(4,1)=( 3.d0,0.d0)
   char_mat(4,2)=( 0.d0,0.d0)
   char_mat(4,5)=(-1.d0,0.d0)
   char_mat(4,3)=(-1.d0,0.d0)
   char_mat(4,6)=( 3.d0,0.d0)
   char_mat(4,7)=( 0.d0,0.d0)
   char_mat(4,8)=(-1.d0,0.d0)
   char_mat(4,10)=(-1.d0,0.d0)


   name_rap(5)="T_2g G_25' G_5+"
   ir_ram(5)="R"
   char_mat(5,1)=( 3.d0,0.d0)
   char_mat(5,2)=( 0.d0,0.d0)
   char_mat(5,3)=(-1.d0,0.d0)
   char_mat(5,4)=(-1.d0,0.d0)
   char_mat(5,6)=( 3.d0,0.d0)
   char_mat(5,7)=( 0.d0,0.d0)
   char_mat(5,8)=(-1.d0,0.d0)
   char_mat(5,9)=(-1.d0,0.d0)

   name_rap(6)="A_1u G_1'  G_1-"
   char_mat(6,6)=(-1.d0,0.d0)
   char_mat(6,7)=(-1.d0,0.d0)
   char_mat(6,8)=(-1.d0,0.d0)
   char_mat(6,9)=(-1.d0,0.d0)
   char_mat(6,10)=(-1.d0,0.d0)

   name_rap(7)="A_2u G_2'  G_2-"
   char_mat(7,4)=(-1.d0,0.d0)
   char_mat(7,5)=(-1.d0,0.d0)
   char_mat(7,6)=(-1.d0,0.d0)
   char_mat(7,7)=(-1.d0,0.d0)
   char_mat(7,8)=(-1.d0,0.d0)

   name_rap(8)="E_u  G_12' G_3-"
   char_mat(8,1)=( 2.d0,0.d0)
   char_mat(8,2)=(-1.d0,0.d0)
   char_mat(8,3)=( 2.d0,0.d0)
   char_mat(8,4)=( 0.d0,0.d0)
   char_mat(8,5)=( 0.d0,0.d0)
   char_mat(8,6)=(-2.d0,0.d0)
   char_mat(8,8)=(-2.d0,0.d0)
   char_mat(8,9)=( 0.d0,0.d0)
   char_mat(8,10)=( 0.d0,0.d0)

   name_rap(9)="T_1u G_15  G_4-"
   ir_ram(9)="I"
   char_mat(9,1)=( 3.d0,0.d0)
   char_mat(9,2)=( 0.d0,0.d0)
   char_mat(9,3)=(-1.d0,0.d0)
   char_mat(9,5)=(-1.d0,0.d0)
   char_mat(9,6)=(-3.d0,0.d0)
   char_mat(9,7)=( 0.d0,0.d0)
   char_mat(9,9)=(-1.d0,0.d0)

   name_rap(10)="T_2u G_25  G_5-"
   char_mat(10,1)=( 3.d0,0.d0)
   char_mat(10,2)=( 0.d0,0.d0)
   char_mat(10,3)=(-1.d0,0.d0)
   char_mat(10,4)=(-1.d0,0.d0)
   char_mat(10,6)=(-3.d0,0.d0)
   char_mat(10,7)=( 0.d0,0.d0)
   char_mat(10,10)=(-1.d0,0.d0)
ELSE
   CALL errore('set_irr_rap','code number not allowed',1)
END IF

RETURN
END SUBROUTINE set_irr_rap

!--------------------------------------------------------------------------
FUNCTION is_complex(code)
!--------------------------------------------------------------------------
! This function receives a code of the group and provide .true. or 
! .false. if the group HAS or HAS NOT complex irreducible 
! representations.
! The order is the following:
!
!   1  "C_1 " F    11 "D_6 " F    21 "D_3h" F    31 "O   " F
!   2  "C_i " F    12 "C_2v" F    22 "D_4h" F    32 "O_h " F 
!   3  "C_s " F    13 "C_3v" F    23 "D_6h" F 
!   4  "C_2 " F    14 "C_4v" F    24 "D_2d" F
!   5  "C_3 " T    15 "C_6v" F    25 "D_3d" F
!   6  "C_4 " T    16 "C_2h" F    26 "S_4 " T
!   7  "C_6 " T    17 "C_3h" T    27 "S_6 " T
!   8  "D_2 " F    18 "C_4h" T    28 "T   " T
!   9  "D_3 " F    19 "C_6h" T    29 "T_h " T
!   10 "D_4 " F    20 "D_2h" F    30 "T_d " F
!
IMPLICIT NONE

INTEGER :: code
LOGICAL :: is_complex

LOGICAL :: complex_aux(32)

data complex_aux  / .FALSE., .FALSE., .FALSE., .FALSE., .TRUE. , &
                    .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., &
                    .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                    .FALSE., .TRUE. , .TRUE. , .TRUE. , .FALSE., &
                    .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                    .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., &
                    .FALSE., .FALSE.  /

IF (code < 1 .OR. code > 32 ) CALL errore('is_complex', &
                                          'code is out of range',1)

is_complex= complex_aux(code)

RETURN
END FUNCTION is_complex

FUNCTION is_parallel(a,b)
!
!  This function returns true if a(3) and b(3) are parallel vectors
!
USE kinds, ONLY : DP
IMPLICIT none
LOGICAL :: is_parallel
REAL(DP) :: a(3), b(3)
REAL(DP) :: cross

cross=(a(2)*b(3)-a(3)*b(2))**2+(a(3)*b(1)-a(1)*b(3))**2+(a(1)*b(2)-a(2)*b(1))**2

is_parallel=(ABS(cross)< 1.d-6)

RETURN
END FUNCTION is_parallel

FUNCTION angle_vectors(ax,bx)
!
!  This function returns the angle, in degrees between two vectors
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
IMPLICIT none
REAL(DP) :: angle_vectors
REAL(DP) :: ax(3), bx(3)
REAL(DP) :: cosangle, moda, modb

moda=sqrt(ax(1)**2+ax(2)**2+ax(3)**2)
modb=sqrt(bx(1)**2+bx(2)**2+bx(3)**2)

IF (moda<1.d-12.OR.modb<1.d-12) &
   CALL errore('angle vectors','zero module vector',1)

cosangle = (ax(1)*bx(1)+ax(2)*bx(2)+ax(3)*bx(3))/moda/modb
angle_vectors = acos(cosangle) * 180.d0 / pi

RETURN
END FUNCTION angle_vectors

SUBROUTINE set_class_el_name(nsym,sname,nclass,nelem,elem,elem_name)

IMPLICIT NONE
INTEGER :: nsym
CHARACTER(LEN=45) :: sname(nsym)
CHARACTER(LEN=55) :: elem_name(8,12)
INTEGER :: nclass, nelem(12), elem(8,12)

INTEGER :: iclass, ielem

DO iclass=1,nclass
   DO ielem=1,nelem(iclass)
      elem_name(ielem,iclass)=sname(elem(ielem,iclass))
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_class_el_name

SUBROUTINE which_c2( a, ia)
!
!   This routine gives a code to identify the direction of a C_2 axis
!   x    1   y=z,  x=0  4   x=-z,  y=0  7   y= m x,  z=0  10   y=-x/m, z=0 13
!   y    2   y=-z, x=0  5   y=x,   z=0  8   y= -m x, z=0  11
!   z    3   x=z,  y=0  6   y=-x,  z=0  9   y= x/m,  z=0  12
!
!   m=sqrt(3.)
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: a(3)
INTEGER, INTENT(OUT) :: ia
REAL(DP) :: epsil=1.D-7, sqr3=SQRT(3.0_DP)
LOGICAL :: is_axis

ia=0
IF (is_axis(a, 1)) THEN
    ia=1
ELSEIF (is_axis(a, 2)) THEN
    ia=2
ELSEIF (is_axis(a, 3)) THEN
    ia=3
ELSEIF (ABS(a(1)) < epsil) THEN
    IF (ABS(a(2)-a(3))< epsil) THEN
       ia=4
    ELSEIF (ABS(a(2)+a(3))< epsil) THEN
       ia=5 
    ELSE
       CALL errore('which_c2','c2 axis not recognized',1)
    ENDIF
ELSEIF (ABS(a(2)) < epsil) THEN
    IF (ABS(a(1)-a(3))< epsil) THEN
       ia=6
    ELSEIF (ABS(a(1)+a(3))< epsil) THEN
       ia=7 
    ELSE
       CALL errore('which_c2','c2 axis not recognized',2)
    ENDIF
ELSEIF (ABS(a(3)) < epsil) THEN
    IF (ABS(a(1)-a(2))< epsil) THEN
       ia=8
    ELSEIF (ABS(a(1)+a(2))< epsil) THEN
       ia=9 
    ELSEIF (ABS(a(1)-a(2)/sqr3) < epsil) THEN
       ia=10
    ELSEIF (ABS(a(1)+a(2)/sqr3) < epsil) THEN
       ia=11
    ELSEIF (ABS(a(1)-a(2)*sqr3) < epsil) THEN
       ia=12
    ELSEIF (ABS(a(1)+a(2)*sqr3) < epsil) THEN
       ia=13
    ELSE
       CALL errore('which_c2','c2 axis not recognized',3)
    ENDIF
ELSE
    CALL errore('which_c2','c2 axis not recognized',4)
ENDIF

RETURN
END SUBROUTINE which_c2

SUBROUTINE is_c2v(iax,ibx,icx,isok)
!
!  This subroutine receives the indeces of a C_2 axis iax, with the convention
!  of the routine which_c2, and of the perpendicular of two mirrors, ibx and
!  icx with the same convention. It returns .true. in isok if the combination
!  corresponds to a known possibility for C_2v. Usually isok .FALSE. means
!  that the two mirrors are not in the correct order.
!
!   Note: the order of the two mirrors in C_2v is defined by the   
!         condition to have the same double group multiplication table of the 
!         isomorphous group D_2. Only the order of one D_2 is arbitrary,
!         all the C_2v and D_2 are ordered for isomorphism of the double
!         groups
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: iax, ibx, icx
LOGICAL, INTENT(OUT) :: isok

isok=.FALSE.
isok = isok .OR. ( iax==1 .AND. ibx==2 .AND. icx==3  )
isok = isok .OR. ( iax==1 .AND. ibx==4 .AND. icx==5  )
isok = isok .OR. ( iax==2 .AND. ibx==3 .AND. icx==1  )
isok = isok .OR. ( iax==2 .AND. ibx==7 .AND. icx==6  ) 
isok = isok .OR. ( iax==3 .AND. ibx==1 .AND. icx==2  ) 
isok = isok .OR. ( iax==3 .AND. ibx==8 .AND. icx==9  )
isok = isok .OR. ( iax==4 .AND. ibx==5 .AND. icx==1  )
isok = isok .OR. ( iax==5 .AND. ibx==1 .AND. icx==4  )
isok = isok .OR. ( iax==6 .AND. ibx==2 .AND. icx==7  )
isok = isok .OR. ( iax==7 .AND. ibx==6 .AND. icx==2  )
isok = isok .OR. ( iax==8 .AND. ibx==3 .AND. icx==9  )
isok = isok .OR. ( iax==9 .AND. ibx==8 .AND. icx==3  )
isok = isok .OR. ( iax==3 .AND. ibx==11.AND. icx==12 )
isok = isok .OR. ( iax==3 .AND. ibx==13.AND. icx==10 )
isok = isok .OR. ( iax==12.AND. ibx==3 .AND. icx==11 )
isok = isok .OR. ( iax==13.AND. ibx==10.AND. icx==3  )
isok = isok .OR. ( iax==10.AND. ibx==3 .AND. icx==13 )
isok = isok .OR. ( iax==11.AND. ibx==12 .AND. icx==3 )

RETURN
END SUBROUTINE is_c2v

SUBROUTINE is_d2(iax, ibx, icx, ind2)
!
!   This routine receives as input the indices of three c2 axes, with the
!   same codes as in which_c2, and gives as output the array ind2(3) of intergers,
!   that gives the positions of iax, ibx, icx in the list C_2, C_2', C_2''
!   ind2(1) says which is the position of iax, ecc.
!   For instance ind2 = 3, 2, 1  says that iax is C_2'', ibx is C_2' and icx is 
!   C_2. If on output ind2 = 0, 0, 0 means that iax, ibx, and icx does not belong
!   to a possible D_2
!
!   Note: this order is arbitrary for one D_2, all the others should be
!         isomorphous with the same double group multiplication table.
!

IMPLICIT NONE
INTEGER, INTENT(IN) :: iax, ibx, icx
INTEGER, INTENT(OUT) :: ind2(3)

   ind2=0
   IF (iax==1) THEN
      IF (ibx==2) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSEIF (ibx==3) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==4) THEN
         ind2(1)=1
         ind2(2)=2
         ind2(3)=3
      ELSEIF (ibx==5) THEN
         ind2(1)=1
         ind2(2)=3
         ind2(3)=2
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',1)
      ENDIF
   ELSEIF (iax==2) THEN
      IF (ibx==1) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSEIF (ibx==3) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==6) THEN
         ind2(1)=1
         ind2(2)=3
         ind2(3)=2
      ELSEIF (ibx==7) THEN
         ind2(1)=1
         ind2(2)=2
         ind2(3)=3
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',2)
      ENDIF
   ELSEIF (iax==3) THEN
      IF (ibx==1 .OR. ibx==9 .OR. ibx==11 .OR. ibx==13) THEN
         ind2(1)=1
         ind2(2)=2
         ind2(3)=3
      ELSEIF (ibx==2 .OR. ibx==8 .OR. ibx==10 .OR. ibx==12) THEN
         ind2(1)=1
         ind2(2)=3
         ind2(3)=2
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',3)
      ENDIF
   ELSEIF (iax==4) THEN
      IF (ibx==1) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==5) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',4)
      ENDIF
   ELSEIF (iax==5) THEN
      IF (ibx==1) THEN
         ind2(1)=3
         ind2(2)=1
         ind2(3)=2
      ELSEIF (ibx==4) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',4)
      ENDIF
   ELSEIF (iax==6) THEN
      IF (ibx==2) THEN
         ind2(1)=3
         ind2(2)=1
         ind2(3)=2
      ELSEIF (ibx==7) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',5)
      ENDIF
   ELSEIF (iax==7) THEN
      IF (ibx==2) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==6) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',6)
      ENDIF
   ELSEIF (iax==8) THEN
      IF (ibx==3) THEN
         ind2(1)=3
         ind2(2)=1
         ind2(3)=2
      ELSEIF (ibx==9) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',7)
      ENDIF
   ELSEIF (iax==9) THEN
      IF (ibx==3) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==8) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',8)
      ENDIF
   ELSEIF (iax==10) THEN
      IF (ibx==3) THEN
         ind2(1)=3
         ind2(2)=1
         ind2(3)=2
      ELSEIF (ibx==13) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',9)
      ENDIF
   ELSEIF (iax==11) THEN
      IF (ibx==3) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==12) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',10)
      ENDIF
   ELSEIF (iax==12) THEN
      IF (ibx==3) THEN
         ind2(1)=3
         ind2(2)=1
         ind2(3)=2
      ELSEIF (ibx==11) THEN
         ind2(1)=3
         ind2(2)=2
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',11)
      ENDIF
   ELSEIF (iax==13) THEN
      IF (ibx==3) THEN
         ind2(1)=2
         ind2(2)=1
         ind2(3)=3
      ELSEIF (ibx==10) THEN
         ind2(1)=2
         ind2(2)=3
         ind2(3)=1
      ELSE
         CALL errore('is_d2','D_2 problem with C_2 axis',12)
      ENDIF
   ELSE
      CALL errore('is_d2','D_2 problem with C_2 axis',9)
   END IF

RETURN
END SUBROUTINE is_d2

