!
! Copyright (C) 2002-2010 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ggenb (b1b, b2b, b3b, nr1b ,nr2b, nr3b, nr1bx ,nr2bx, nr3bx, gcutb )
!-----------------------------------------------------------------------
   !
   ! As ggen, for the box grid. A "b" is appended to box variables.
   ! The documentation for ggen applies
   !
   USE kinds, ONLY: DP
   use gvecb, only: ngb, ngbt, ngbl, ngbx, gb, gxb, glb, npb, nmb, mill_b
   use io_global, only: stdout, ionode
   use control_flags, only: iprsta
!
   implicit none
!
   integer nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
   REAL(DP) b1b(3), b2b(3), b3b(3), gcutb
!
   integer, allocatable:: idx(:), iglb(:)
   integer n1pb, n2pb, n3pb, n1mb, n2mb, n3mb
   integer it, icurr, nr1m1, nr2m1, nr3m1, ir, ig, i,j,k, itv(3), ip
   REAL(DP) t(3), g2
!
      nr1m1=nr1b-1
      nr2m1=nr2b-1
      nr3m1=nr3b-1
      ngb=0
!
!     first step : count the number of vectors with g2 < gcutb
!
!     exclude space with x<0
!
      do i= 0,nr1m1
         do j=-nr2m1,nr2m1
!
!     exclude plane with x=0, y<0
!
            if(i.eq.0.and.j.lt.0) go to 10
!
            do k=-nr3m1,nr3m1
!
!     exclude line with x=0, y=0, z<0
!
               if(i.eq.0.and.j.eq.0.and.k.lt.0) go to 20
               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1b(ir) + DBLE(j)*b2b(ir) + DBLE(k)*b3b(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2.gt.gcutb) go to 20
               ngb=ngb+1
 20            continue
            end do
 10         continue
         end do
      end do
!
!     second step: allocate space
!
      allocate(gxb(3,ngb))
      allocate(gb(ngb))
      allocate(npb(ngb))
      allocate(nmb(ngb))
      allocate(iglb(ngb))
      allocate(mill_b(3,ngb))
      allocate(idx(ngb))
!
!     third step : find the vectors with g2 < gcutb
!
      ngb=0
!
!     exclude space with x<0
!
      do i= 0,nr1m1
         do j=-nr2m1,nr2m1
!
!     exclude plane with x=0, y<0
!
            if(i.eq.0.and.j.lt.0) go to 15
!
            do k=-nr3m1,nr3m1
!
!     exclude line with x=0, y=0, z<0
!
               if(i.eq.0.and.j.eq.0.and.k.lt.0) go to 25
               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1b(ir) + DBLE(j)*b2b(ir) + DBLE(k)*b3b(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2.gt.gcutb) go to 25
               ngb=ngb+1
               gb(ngb)=g2
               mill_b(1,ngb)=i
               mill_b(2,ngb)=j
               mill_b(3,ngb)=k
 25            continue
            end do
 15         continue
         end do
      end do

      IF( iprsta > 3 ) THEN
        WRITE( stdout,*)
        WRITE( stdout,170) ngb
 170    format(' ggenb: # of gb vectors < gcutb ngb = ',i6)
      END IF

      idx(1)=0
      call hpsort (ngb,gb,idx)

      do ig=1,ngb-1
         icurr=ig
 30      if(idx(icurr).ne.ig) then
            itv=mill_b(:,icurr)
            mill_b(:,icurr)=mill_b(:,idx(icurr))
            mill_b(:,idx(icurr))=itv

            it=icurr
            icurr=idx(icurr)
            idx(it)=it
            if(idx(icurr).eq.ig) then
               idx(icurr)=icurr
               goto 35
            endif
            goto 30
         endif
 35      continue
      end do
!
      deallocate(idx)
!
! costruct fft indexes (n1b,n2b,n3b) for the box grid
!
      do ig=1,ngb
         i=mill_b(1,ig)
         j=mill_b(2,ig)
         k=mill_b(3,ig)
         n1pb=i+1
         n2pb=j+1
         n3pb=k+1
!
! n1pb,n2pb,n3pb: indexes of G
! negative indexes are refolded (note that by construction i.ge.0)
!
         if(i.lt.0) n1pb=n1pb+nr1b
         if(j.lt.0) n2pb=n2pb+nr2b
         if(k.lt.0) n3pb=n3pb+nr3b
!
! n1mb,n2mb,n3mb: indexes of -G
!
         if(i.eq.0) then
            n1mb=1
         else
            n1mb=nr1b-n1pb+2
         end if
         if(j.eq.0) then
            n2mb=1
         else
            n2mb=nr2b-n2pb+2
         end if
         if(k.eq.0) then
            n3mb=1
         else
            n3mb=nr3b-n3pb+2
         end if
!
! conversion from (i,j,k) index to combined 1-d ijk index:
! ijk = 1 + (i-1)+(j-1)*ix+(k-1)*ix*jx
! where the (i,j,k) array is assumed to be dimensioned (ix,jx,kx)
!
         npb(ig) = n1pb+(n2pb-1)*nr1bx+(n3pb-1)*nr1bx*nr2bx
         nmb(ig) = n1mb+(n2mb-1)*nr1bx+(n3mb-1)*nr1bx*nr2bx
      end do
!
! shells of G - first calculate their number and position
!

      CALL gshcount( ngbl, iglb, ngb, gb, -1.0d0, -1.0d0 )

      IF( iprsta > 3 ) THEN
        WRITE( stdout,180) ngbl
 180    format(' ggenb: # of gb shells  < gcutb ngbl= ',i6)
      END IF
!
! then allocate the array glb
!
      allocate(glb(ngbl))
!
! and finally fill glb with the values of the shells
!
      glb(iglb(1))=gb(1)
      do ig=2,ngb
         if(iglb(ig).ne.iglb(ig-1)) glb(iglb(ig))=gb(ig)
      end do
!
! calculation of G-vectors
!
      do ig=1,ngb
         i=mill_b(1,ig)
         j=mill_b(2,ig)
         k=mill_b(3,ig)
         gxb(1,ig)=i*b1b(1)+j*b2b(1)+k*b3b(1)
         gxb(2,ig)=i*b1b(2)+j*b2b(2)+k*b3b(2)
         gxb(3,ig)=i*b1b(3)+j*b2b(3)+k*b3b(3)
      end do
!
      DEALLOCATE (iglb)
      return
end subroutine ggenb



!-----------------------------------------------------------------------
      subroutine gcalb( alatb, b1b_ , b2b_ , b3b_  )
!-----------------------------------------------------------------------
!
      USE kinds, ONLY: DP
      use gvecb
!
      implicit none
      REAL(DP), intent(in) :: alatb, b1b_ (3), b2b_ (3), b3b_ (3)
      REAL(DP) :: b1b(3), b2b(3), b3b(3)
!
      integer i, i1,i2,i3,ig

      b1b = b1b_ * alatb
      b2b = b2b_ * alatb
      b3b = b3b_ * alatb
!
!     calculation of gxb(3,ngbx)
!
      do ig=1,ngb
         i1=mill_b(1,ig)
         i2=mill_b(2,ig)
         i3=mill_b(3,ig)
         gxb(1,ig)=i1*b1b(1)+i2*b2b(1)+i3*b3b(1)
         gxb(2,ig)=i1*b1b(2)+i2*b2b(2)+i3*b3b(2)
         gxb(3,ig)=i1*b1b(3)+i2*b2b(3)+i3*b3b(3)
         gb(ig)=gxb(1,ig)**2 + gxb(2,ig)**2 + gxb(3,ig)**2
      enddo
!
      return
      end subroutine gcalb
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE gshcount( ngl, igl, ng, gg, gcuts, gcutw )
!-------------------------------------------------------------------------

  USE kinds,     ONLY: DP

  IMPLICIT NONE

  INTEGER :: ngl
  INTEGER :: igl(*)
  INTEGER :: ng
  REAL(DP) :: gg(*), gcuts, gcutw

  INTEGER :: ig

      ngl=1
      igl(1)=ngl
      do ig=2,ng
         if(abs(gg(ig)-gg(ig-1)).gt.1.e-6)then
            ngl=ngl+1
            !!! if (gg(ig).lt.gcuts) ngsl=ngl
            !!! if (gg(ig).lt.gcutw) ngwl=ngl
         endif
         igl(ig)=ngl
      end do

  RETURN
END SUBROUTINE gshcount

!-------------------------------------------------------------------------
      subroutine gcal( bg )
!-----------------------------------------------------------------------
!   calculates the values of g-vectors to be assigned to the lattice
!   points generated in subroutine ggen. these values are derived
!   from the actual values of lattice parameters, with fixed number
!   of plane waves and a cut-off function to keep energy cut-off fixed.
!
!      g=i*b1+j*b2+k*b3,
!
!   where b1,b2,b3 are the vectors defining the reciprocal lattice,
!   i go from 1 to +(nr-1) and j,k go from -(nr-1) to +(nr-1).
!
!   the g's are in units of 2pi/a.
!
      USE kinds,     ONLY: DP
      use gvect, only: ngm, gg, g, mill
      implicit none
!
      REAL(DP), INTENT (IN) :: bg(3,3)
!
      integer i1,i2,i3,ig
!
!     calculation of g(3,ng)
!
      do ig=1,ngm
         i1=mill(1,ig)
         i2=mill(2,ig)
         i3=mill(3,ig)
         g(:,ig)=i1*bg(:,1)+i2*bg(:,2)+i3*bg(:,3)
         gg(ig)=g(1,ig)**2 + g(2,ig)**2 + g(3,ig)**2
      enddo
 
      return
      end subroutine gcal

!=----------------------------------------------------------------------------=!

        SUBROUTINE newgb( omega, alat, at )
!
!     re-generation of little box g-vectors
!
          USE kinds, ONLY: DP
          USE grid_dimensions, only: nr1, nr2, nr3
          USE smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
          USE small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
          USE constants, ONLY: pi

          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: at(3,3), omega, alat

          INTEGER :: i
          REAL(DP) :: alatb, b1b(3),b2b(3),b3b(3)

          IF ( nr1b == 0 .OR. nr2b == 0 .OR. nr3b == 0 ) return
          alatb  = alat / nr1*nr1b
          tpibab = 2.d0*pi / alatb
          do i=1,3
            a1b(i)=at(i,1)*alat/nr1*nr1b
            a2b(i)=at(i,2)*alat/nr2*nr2b
            a3b(i)=at(i,3)*alat/nr3*nr3b
          enddo

          omegab=omega/nr1*nr1b/nr2*nr2b/nr3*nr3b
!
          call recips( a1b, a2b, a3b, b1b, b2b, b3b )
          !
          call gcalb( alatb, b1b, b2b, b3b )
!
          do i=1,3
            ainvb(1,i)=b1b(i)
            ainvb(2,i)=b2b(i)
            ainvb(3,i)=b3b(i)
          end do

          RETURN
        END SUBROUTINE newgb

!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!

        SUBROUTINE ecutoffs_setup( ecutwfc_, ecutrho_, ecfixed_, qcutz_, &
                                   q2sigma_, refg_ )
 
          USE kinds,           ONLY: DP
          USE constants,       ONLY: eps8
          USE gvecw,           ONLY: ecutwfc
          USE gvecw,           ONLY: ecfixed, qcutz, q2sigma
          USE gvect,           ONLY: ecutrho
          USE gvecs,           ONLY: ecuts, dual, doublegrid
          use betax,           only: mmx, refg
          USE pseudopotential, only: tpstab
          USE control_flags,   only: thdyn
          USE io_global,       only: stdout, ionode
          USE uspp,            only: okvan

          IMPLICIT NONE
          REAL(DP), INTENT(IN) ::  ecutwfc_, ecutrho_, ecfixed_, qcutz_, &
                                   q2sigma_, refg_

          ecutwfc = ecutwfc_

          IF ( ecutrho_ <= 0.D0 ) THEN
             !
             dual = 4.D0
             !
          ELSE
             !
             dual = ecutrho_ / ecutwfc
             !
             IF ( dual <= 1.D0 ) &
                CALL errore( ' ecutoffs_setup ', ' invalid dual? ', 1 )
             !
          END IF

          doublegrid = ( dual > 4.D0 )
          IF ( doublegrid .AND. .NOT. okvan ) &
             CALL errore( 'setup', 'No USPP: set ecutrho=4*ecutwfc', 1 )
          ecutrho = dual * ecutwfc
          !
          IF ( doublegrid ) THEN
             !
             ecuts = 4.D0 * ecutwfc
             !
          ELSE
             !
             ecuts = ecutrho
             !
          END IF

          !
          ecfixed = ecfixed_
          qcutz   = qcutz_
          q2sigma = q2sigma_

          IF( refg_ < 0.0001d0 ) THEN
             tpstab = .FALSE.
             refg = 0.05d0
          ELSE
             refg = refg_
          END IF

          IF( thdyn ) THEN
             !  ... a larger table is used when cell is moving to allow 
             !  ... large volume fluctuation
             mmx  = NINT( 2.0d0 * ecutrho / refg )
          ELSE
             mmx  = NINT( 1.2d0 * ecutrho / refg )
          END IF

             mmx  = NINT( 2.0d0 * ecutrho / refg ) ! debug

          RETURN
        END SUBROUTINE ecutoffs_setup


        SUBROUTINE gcutoffs_setup( alat, tk_inp, nk_inp, kpoints_inp )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

          USE kinds, ONLY: DP
          USE gvecw, ONLY: ecutwfc,  gcutw
          USE gvect, ONLY: ecutrho,  gcutm
          USE gvecs, ONLY: ecuts, gcutms
          USE gvecb, ONLY: ecutb, gcutb
          USE gvecw, ONLY: ekcut, gkcut
          USE constants, ONLY: eps8, pi

          IMPLICIT NONE

! ...     declare subroutine arguments
          REAL(DP), INTENT(IN) :: alat
          LOGICAL, INTENT(IN) :: tk_inp
          INTEGER, INTENT(IN) :: nk_inp
          REAL(DP), INTENT(IN) :: kpoints_inp(3,*)

! ...     declare other variables
          INTEGER   :: i
          REAL(DP) :: kcut, ksq
          REAL(DP) :: tpiba

!  end of declarations
!  ----------------------------------------------

! ...   Set Values for the cutoff


          IF( alat < eps8 ) THEN
            CALL errore(' cut-off setup ', ' alat too small ', 0)
          END IF

          tpiba = 2.0d0 * pi / alat 

          ! ...  Constant cutoff simulation parameters

          gcutw = ecutwfc / tpiba**2  ! wave function cut-off
          gcutm = ecutrho / tpiba**2  ! potential cut-off
          gcutms= ecuts   / tpiba**2  ! smooth mesh cut-off

          kcut = 0.0_DP
          IF ( tk_inp ) THEN
! ...       augment plane wave cutoff to include all k+G's
            DO i = 1, nk_inp
! ...         calculate modulus
              ksq = kpoints_inp( 1, i ) ** 2 + kpoints_inp( 2, i ) ** 2 + kpoints_inp( 3, i ) ** 2
              IF ( ksq > kcut ) kcut = ksq
            END DO
          END IF

          gkcut = ( sqrt( kcut ) + sqrt( gcutw ) ) ** 2

          ekcut = gkcut * tpiba ** 2

          RETURN
        END SUBROUTINE gcutoffs_setup

!  ----------------------------------------------

      SUBROUTINE cutoffs_print_info()

        !  Print out informations about different cut-offs

        USE gvecw, ONLY: ecutwfc,  gcutw
        USE gvect, ONLY: ecutrho,  gcutm
        USE gvecw, ONLY: ecfixed, qcutz, q2sigma
        USE gvecw, ONLY: ekcut, gkcut
        USE gvecs, ONLY: ecuts, gcutms
        USE gvecb, ONLY: ecutb, gcutb
        use betax, only: mmx, refg
        USE io_global, ONLY: stdout

        WRITE( stdout, 100 ) ecutwfc, ecutrho, ecuts, sqrt(gcutw), &
                             sqrt(gcutm), sqrt(gcutms)
        IF( qcutz > 0.0d0 ) THEN
          WRITE( stdout, 150 ) qcutz, q2sigma, ecfixed
        END IF

        WRITE( stdout,200) refg, mmx

100     FORMAT(/,3X,'Energy Cut-offs',/ &
                ,3X,'---------------',/ &
                ,3X,'Ecutwfc = ',F6.1,' Ry,   ', 3X,'Ecutrho = ',F6.1,' Ry,   ', 3X,'Ecuts = ',F6.1,' Ry',/ &
                ,3X,'Gcutwfc = ',F6.1,'     , ', 3X,'Gcutrho = ',F6.1,'       ', 3X,'Gcuts = ',F6.1)
150     FORMAT(  3X,'modified kinetic energy functional, with parameters:',/,   &
                 3X,'ecutz = ',f8.4,'  ecsig = ', f7.4,'  ecfix = ',f6.2)
200     FORMAT(  3X,'NOTA BENE: refg, mmx = ', f10.6,I6 )

        RETURN
      END SUBROUTINE cutoffs_print_info

!  ----------------------------------------------

      SUBROUTINE orthogonalize_info( )
        USE control_flags, ONLY: ortho_eps, ortho_max
        USE io_global, ONLY: stdout
        IMPLICIT NONE
           WRITE(stdout, 585)
           WRITE(stdout, 511) ortho_eps, ortho_max
  511   FORMAT(   3X,'Orthog. with lagrange multipliers : eps = ',E10.2, ',  max = ',I3)
  585   FORMAT(   3X,'Eigenvalues calculated without the kinetic term contribution')
        RETURN
      END SUBROUTINE orthogonalize_info


!  ----------------------------------------------


      SUBROUTINE electrons_print_info( )

          USE kinds, ONLY: DP
          USE electrons_base, ONLY: nbnd, nspin, nel, nelt, nupdwn, iupdwn, &
                                    f, qbac
          USE io_global, ONLY: stdout
          USE ions_base, ONLY: zv, nsp, na

          IMPLICIT NONE
          INTEGER :: i,is

          IF( nspin == 1) THEN
            WRITE(stdout,6) nelt, nbnd
            WRITE(stdout,7) ( f( i ), i = 1, nbnd )
          ELSE
            WRITE(stdout,8) nelt
            WRITE(stdout,9) nel(1)
            WRITE(stdout,7) ( f( i ), i = 1, nupdwn(1))
            WRITE(stdout,10) nel(2)
            WRITE(stdout,7) ( f( i ), i = iupdwn(2), ( iupdwn(2) + nupdwn(2) - 1 ) )
          END IF

         qbac=0.
         do is=1,nsp
           qbac=qbac+na(is)*zv(is)
         end do
         qbac=qbac-nelt
         if(qbac.ne.0) write(stdout,11) qbac


6         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Number of Electron = ',I5,', of States = ',I5,/ &
                  ,3X,'Occupation numbers :')
7         FORMAT(2X,10F5.2)
8         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Local Spin Density calculation',/ &
                  ,3X,'Number of Electron = ',I5)
9         FORMAT(  3X,'Spins up   = ', I5, ', occupations: ')
10        FORMAT(  3X,'Spins down = ', I5, ', occupations: ')
11        FORMAT(/,3X,'WARNING: system charge = ',F12.6)
          RETURN
      END SUBROUTINE electrons_print_info


!  ----------------------------------------------


      SUBROUTINE exch_corr_print_info()

        USE funct, ONLY: get_iexch, get_icorr, get_igcx, get_igcc, write_dft_name
        USE io_global, ONLY: stdout

        IMPLICIT NONE

        CHARACTER(LEN = 60) :: exch_info
        CHARACTER(LEN = 60) :: corr_info
        CHARACTER(LEN = 60) :: exgc_info
        CHARACTER(LEN = 60) :: cogc_info

        WRITE(stdout,800)

          ! ...     iexch => Exchange functional form
          ! ...     icorr => Correlation functional form
          ! ...     igcx  => Gradient Correction to the Exchange potential
          ! ...     igcc  => Gradient Correction to the Correlation potential

          SELECT CASE ( get_iexch() )
            CASE (0)
              exch_info = 'NONE'
            CASE (1)
              exch_info = 'SLATER'
            CASE (2)
              exch_info = 'SLATER (alpha=1)'
            CASE DEFAULT
              exch_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_icorr() )
            CASE (0)
              corr_info = 'NONE'
            CASE (1)
              corr_info = 'PERDEW AND ZUNGER'
            CASE (2)
              corr_info = 'VOSKO, WILK AND NUSAIR'
            CASE (3)
              corr_info = 'LEE, YANG, AND PARR'
            CASE (4)
              corr_info = 'PERDEW AND WANG'
            CASE (9)
              corr_info = 'PADE APPROXIMATION'
            CASE DEFAULT
              corr_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_igcx() )
            CASE (0)
              exgc_info = 'NONE'
            CASE (1)
              exgc_info = 'BECKE'
            CASE (2)
              exgc_info = 'PERDEW'
            CASE (3)
              exgc_info = 'PERDEW BURKE ERNZERHOF'
            CASE (7)
              exgc_info = 'META-TPSS'
            CASE DEFAULT
              exgc_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_igcc() )
            CASE (0)
              cogc_info = 'NONE'
            CASE (1)
              cogc_info = 'PERDEW'
            CASE (2)
              cogc_info = 'LEE, YANG AND PARR'
            CASE (3)
              cogc_info = 'PERDEW AND WANG'
            CASE (4)
              cogc_info = 'PERDEW BURKE ERNZERHOF'
            CASE (6)
              cogc_info = 'META-TPSS'
            CASE DEFAULT
              cogc_info = 'UNKNOWN'
          END SELECT

          WRITE(stdout,910)
          WRITE(stdout,fmt='(5X,"Exchange functional: ",A)') exch_info
          WRITE(stdout,fmt='(5X,"Correlation functional: ",A)') corr_info
          IF( ( get_igcx() > 0 ) .OR. ( get_igcc() > 0 ) ) THEN
            WRITE(stdout,810)
            WRITE(stdout,fmt='(5X,"Exchange functional: ",A)') exgc_info
            WRITE(stdout,fmt='(5X,"Correlation functional: ",A)') cogc_info
          END IF

        call write_dft_name

800 FORMAT(//,3X,'Exchange and correlations functionals',/ &
             ,3X,'-------------------------------------')
810 FORMAT(   3X,'Using Generalized Gradient Corrections with')
910 FORMAT(   3X,'Using Local Density Approximation with')

        RETURN
      END SUBROUTINE exch_corr_print_info



!  ----------------------------------------------



       SUBROUTINE ions_print_info( )
            
         !  Print info about input parameter for ion dynamic

         USE io_global,     ONLY: ionode, stdout
         USE control_flags, ONLY: tranp, amprp, tnosep, tolp, tfor, tsdp, tzerop, &
                                  tv0rd, taurdr, nv0rd, nbeg, tcp, tcap
         USE ions_base,     ONLY: tau_srt, if_pos, ind_srt, nsp, na, &
                                  pmass, nat, fricp, greasp, rcmax
         USE ions_nose,     ONLY: tempw, ndega
         USE constants,     ONLY: amu_au

         IMPLICIT NONE
              
         integer is, ia, k, ic, isa
         LOGICAL :: ismb( 3 ) 
                
         WRITE( stdout, 50 ) 

         IF( .NOT. tfor ) THEN
           WRITE( stdout, 518 )
         ELSE
           WRITE( stdout, 520 )
           IF( tsdp ) THEN
             WRITE( stdout, 521 )
           ELSE
             WRITE( stdout, 522 )
           END IF
           WRITE( stdout, 523 ) ndega
           WRITE( stdout, 524 ) fricp, greasp
           IF( tzerop ) then
             IF( tv0rd ) THEN
               WRITE( stdout, 850 ) nv0rd
             ELSE
               WRITE( stdout, 635 )
             ENDIF 
           ENDIF
         END IF 
              
         DO is = 1, nsp
           IF( tranp(is) ) THEN
             WRITE( stdout,510)
             WRITE( stdout,512) is, amprp(is)
           END IF
         END DO

         WRITE(stdout,660) 
         isa = 0
         DO IS = 1, nsp
           WRITE(stdout,1000) is, na(is), pmass(is), pmass(is) / amu_au, rcmax(is)
           DO IA = 1, na(is)
             isa = isa + 1
             WRITE(stdout,1010) ( tau_srt(k,isa), K = 1,3 )
           END DO
         END DO    

         IF ( ( nbeg > -1 ) .AND. ( .NOT. taurdr ) ) THEN
            WRITE(stdout,661)
         ELSE
            WRITE(stdout,662)
         ENDIF

         IF( tfor ) THEN

            IF( ANY( ( if_pos( 1:3, 1:nat ) == 0 )  ) ) THEN

              WRITE(stdout,1020)
              WRITE(stdout,1022)

              DO isa = 1, nat
                ia = ind_srt( isa )
                ismb( 1 ) = ( if_pos(1,ia) /= 0 )
                ismb( 2 ) = ( if_pos(2,ia) /= 0 )
                ismb( 3 ) = ( if_pos(3,ia) /= 0 )
                IF( .NOT. ALL( ismb ) ) THEN
                  WRITE( stdout, 1023 ) isa, ( ismb(k), K = 1, 3 )
                END IF
              END DO

            ELSE

              WRITE(stdout,1021)

            END IF
         END IF

         IF( tfor ) THEN
           if( ( tcp .or. tcap .or. tnosep ) .and. tsdp ) then
             call errore(' ions_print_info',' t contr. for ions when tsdp=.t.',1)
           endif
           IF(.not. tcp .and. .not. tcap .and. .not. tnosep ) THEN
              WRITE( stdout,550)
           ELSE IF( tcp .and. tcap ) then
             call errore(' ions_print_info',' tcp and tcap both true',1)
           ELSE IF( tcp .and. tnosep ) then
             call errore(' ions_print_info',' tcp and tnosep both true',1)
           ELSE IF(tcap .and. tnosep ) then
             call errore(' ions_print_info',' tcap and tnosep both true',1)
           ELSE IF(tcp) THEN
             WRITE( stdout,555) tempw,tolp
           ELSE IF(tcap) THEN
             WRITE( stdout,560) tempw,tolp
           ELSE IF(tnosep) THEN
             WRITE( stdout,595)
           ELSE
             WRITE( stdout,550)
           END IF
         END IF

   50 FORMAT(//,3X,'Ions Simulation Parameters',/ &
               ,3X,'--------------------------')

  510 FORMAT(   3X,'Initial random displacement of ionic coordinates',/, & 
                3X,' specie  amplitude')
  512 FORMAT(   3X,I7,2X,F9.6)

  518 FORMAT(   3X,'Ions are not allowed to move')
  520 FORMAT(   3X,'Ions are allowed to move')
  521 FORMAT(   3X,'Ions dynamics with steepest descent')
  522 FORMAT(   3X,'Ions dynamics with newton equations')
  523 format(   3X,'the temperature is computed for ',i5,' degrees of freedom')
  524 format(   3X,'ion dynamics with fricp = ',f7.4,' and greasp = ',f7.4)
  550 FORMAT(   3X,'Ionic temperature is not controlled')
  555 FORMAT(   3X,'Ionic temperature control via ', &
                   'rescaling of velocities :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  560 FORMAT(   3X,'Ionic temperature control via ', &
                   'canonical velocities rescaling :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  595 FORMAT(   3X,'Ionic temperature control via nose thermostat')
  635 FORMAT(   3X,'Zero initial momentum for ions')

  660 FORMAT(   3X,'Ionic position (from input)', /, &
                3X,'sorted by specie, and converted to real a.u. coordinates')
  661 FORMAT(   3X,'Ionic position will be re-read from restart file')
  662 FORMAT(   3X,'Ionic position read from input file')

  850 FORMAT(   3X,'Initial ion velocities read from unit : ',I4)

 1000 FORMAT(3X,'Species ',I3,' atoms = ',I4,' mass = ',F12.2, ' (a.u.), ', &
               & F12.2, ' (amu)', ' rcmax = ', F6.2, ' (a.u.)' )
 1010 FORMAT(3X,3(1X,F12.6))
 1020 FORMAT(/,3X,'NOT all atoms are allowed to move ')
 1021 FORMAT(/,3X,'All atoms are allowed to move')
 1022 FORMAT(  3X,' indx  ..x.. ..y.. ..z..')
 1023 FORMAT(  3X,I4,3(1X,L5))



         RETURN
       END SUBROUTINE ions_print_info


!  ----------------------------------------------

        subroutine cell_print_info( )

          USE constants, ONLY: au_gpa
          USE control_flags, ONLY: thdyn, tsdc, tzeroc, tbeg, nbeg, tpre
          USE control_flags, ONLY: tnoseh
          USE io_global, ONLY: stdout
          USE cell_base, ONLY: press, frich, greash, wmass

          IMPLICIT NONE

          WRITE(stdout,545 )
          IF ( tpre ) WRITE( stdout, 600 )
          IF ( tbeg ) THEN
            WRITE(stdout,546)
          ELSE
            WRITE(stdout,547)
            IF( nbeg > -1 ) WRITE( stdout, 548 )
          END IF

          IF( .NOT. thdyn ) THEN
            WRITE( stdout,525)
            WRITE( stdout,606)
          ELSE
            IF( tsdc ) THEN
              WRITE( stdout,526)
            ELSE
              IF( frich /= 0.0d0 ) THEN
                WRITE( stdout,602) frich, greash
              ELSE
                WRITE( stdout,527)
              END IF
              IF( tnoseh ) then
                WRITE( stdout,604) 
              ELSE
                WRITE( stdout,565)
              END IF
              IF( tzeroc ) THEN
                WRITE( stdout,563)
              ENDIF
            END IF
            WRITE( stdout,530) press * au_gpa, wmass
          END IF


 545     FORMAT(//,3X,'Cell Dynamics Parameters (from STDIN)',/ &
                  ,3X,'-------------------------------------')
 546     FORMAT(   3X,'Simulation cell read from STDIN')
 547     FORMAT(   3X,'Starting cell generated from CELLDM')
 548     FORMAT(   3X,'Cell parameters will be re-read from restart file')
 525     FORMAT(   3X,'Constant VOLUME Molecular dynamics')
 606     format(   3X,'cell parameters are not allowed to move')
 526     FORMAT(   3X,'Volume dynamics with steepest descent')
 527     FORMAT(   3X,'Volume dynamics with newton equations')
 530     FORMAT(   3X,'Constant PRESSURE Molecular dynamics:',/ &
                  ,3X,'External pressure (GPa) = ',F11.2,/ &
                  ,3X,'Volume mass             = ',F11.2)
 563     FORMAT(   3X,'Zero initial momentum for cell variables')
 565     FORMAT(   3X,'Volume dynamics: the temperature is not controlled')
 604     format(   3X,'cell parameters dynamics with nose` temp. control' )

 600  format( 3X, 'internal stress tensor calculated')
 602  format( 3X, 'cell parameters dynamics with frich = ',f7.4,            &
     &        3X, 'and greash = ',f7.4 )

        return
      end subroutine cell_print_info


!----------------------------------------------
SUBROUTINE gmeshinfo( )
!----------------------------------------------
   !
   !   Print out the number of g vectors for the different mesh
   !
   USE kinds,     ONLY: DP
   USE mp_global, ONLY: nproc_bgrp, intra_bgrp_comm
   USE io_global, ONLY: ionode, ionode_id, stdout
   USE mp,        ONLY: mp_max, mp_gather
   use gvecb,     only: ngb
   USE gvecw,     only: ngw_g, ngw, ngwx
   USE gvecs,     only: ngms_g, ngms, ngsx
   USE gvect,     only: ngm, ngm_g, ngmx

   IMPLICIT NONE

   INTEGER :: ip, ng_snd(3), ng_rcv( 3, nproc_bgrp )
   INTEGER :: ierr, min_val, max_val, i
   REAL(DP) :: avg_val

   IF(ionode) THEN
      WRITE( stdout,*)
      WRITE( stdout,*) '  Reciprocal Space Mesh'
      WRITE( stdout,*) '  ---------------------'
   END IF

   ng_snd(1) = ngm_g
   ng_snd(2) = ngm
   ng_snd(3) = ngmx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   IF(ionode) THEN
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1000)
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
   END IF
   !
   ng_snd(1) = ngms_g
   ng_snd(2) = ngms
   ng_snd(3) = ngsx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   ierr = 0
   !
   IF(ionode) THEN
      WRITE( stdout,1001)
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
      IF( min_val < 1 ) ierr = ip
   END IF
   !
   CALL mp_max( ierr, intra_bgrp_comm )
   !
   IF( ierr > 0 ) &
      CALL errore( " gmeshinfo ", " Wow! some processors have no G-vectors ", ierr )
   !
   ng_snd(1) = ngw_g
   ng_snd(2) = ngw
   ng_snd(3) = ngwx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   IF(ionode) THEN
      WRITE( stdout,1002)
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
      IF( min_val < 1 ) ierr = ip
   END IF
   !
   CALL mp_max( ierr, intra_bgrp_comm )
   !
   IF( ierr > 0 ) &
      CALL errore( " gmeshinfo ", " Wow! some processors have no G-vectors ", ierr )
   !
   IF(ionode .AND. ngb > 0 ) THEN
      WRITE( stdout,1050)
      WRITE( stdout,1060) ngb
   END IF

   1000    FORMAT(3X,'Large Mesh',/, &
           '     Global(ngm_g)    MinLocal       MaxLocal      Average') 
   1001    FORMAT(3X,'Smooth Mesh',/, &
           '     Global(ngms_g)   MinLocal       MaxLocal      Average') 
   1002    FORMAT(3X,'Wave function Mesh',/, &
           '     Global(ngw_g)    MinLocal       MaxLocal      Average') 
   1011    FORMAT(  3I15, F15.2 )
   1050    FORMAT(/,3X,'Small box Mesh')
   1060    FORMAT( 3X, 'ngb = ', I12, ' not distributed to processors' )

   RETURN

END SUBROUTINE gmeshinfo

!----------------------------------------------
SUBROUTINE constraint_info()
!----------------------------------------------
   USE kinds,              ONLY: DP
   USE constraints_module, ONLY: nconstr, constr_tol, &
                                 constr_type, constr, constr_target
   USE io_global,          ONLY: ionode, stdout
   USE control_flags,      ONLY: lconstrain
   !
   IMPLICIT NONE
   !
   INTEGER :: ic
   !
   IF( lconstrain .AND. ionode ) THEN
      !
      WRITE( stdout, 10 ) 
      WRITE( stdout, 20 ) nconstr, constr_tol
      !
      DO ic = 1, nconstr
         !
         IF( constr_type( ic ) == 3 ) THEN
            !
            ! distance
            !
            WRITE( stdout, 30 ) ic
            WRITE( stdout, 40 ) NINT( constr(1,ic) ), &
                                NINT( constr(2,ic) ), constr_target(ic)
            !
         END IF
         !
      END DO
      !
   END IF
   !
10 FORMAT( 3X, "Using constrained dynamics")
20 FORMAT( 3X, "number of constrain and tolerance: ", I5, D10.2)
30 FORMAT( 3X, "constrain ", I5, " type distance ")
40 FORMAT( 3X, "  atoms ", I5, I5, " target dist ", F10.5)
   !
END SUBROUTINE constraint_info


SUBROUTINE new_atomind_constraints()
   !
   USE kinds,              ONLY: DP
   USE constraints_module, ONLY: constr
   USE ions_base,          ONLY: ind_bck
   !
   IMPLICIT NONE
   !
   INTEGER  :: ic, ia
   INTEGER  :: iaa
   REAL(DP) :: aa
   !
   !  Substitute the atom index given in the input file
   !  with the new atom index, after the sort in the
   !  atomic coordinates.
   !
   DO ic = 1, SIZE( constr, 2 )
      DO ia = 1, SIZE( constr, 1 )
         IF( constr( ia, ic ) > 0.0d0 ) THEN
            iaa = NINT( constr( ia, ic ) )
            aa  = DBLE( ind_bck( iaa ) )
            constr( ia, ic ) = aa
         END IF
      END DO
   END DO
   !
   RETURN
   !
END SUBROUTINE new_atomind_constraints


SUBROUTINE compute_stress_x( stress, detot, h, omega )
   USE kinds, ONLY : DP
   IMPLICIT NONE
   REAL(DP), INTENT(OUT) :: stress(3,3)
   REAL(DP), INTENT(IN)  :: detot(3,3), h(3,3), omega
   integer :: i, j
   do i=1,3 
      do j=1,3
         stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+              &
     &                      detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
      enddo
   enddo
   return
END SUBROUTINE compute_stress_x
!-----------------------------------------------------------------------
subroutine formf( tfirst, eself )
  !-----------------------------------------------------------------------

  !computes (a) the self-energy eself of the ionic pseudocharges;
  !         (b) the form factors of: (i) pseudopotential (vps),
  !             (ii) ionic pseudocharge (rhops)
  !         also calculated the derivative of vps with respect to
  !         g^2 (dvps)
  ! 
  USE kinds,           ONLY : DP
  use mp,              ONLY : mp_sum
  use control_flags,   ONLY : iprint, tpre, iprsta
  use io_global,       ONLY : stdout
  use mp_global,       ONLY : intra_bgrp_comm
  use gvecs,           ONLY : ngms
  use cell_base,       ONLY : omega, tpiba2, tpiba
  use ions_base,       ONLY : rcmax, zv, nsp, na
  use local_pseudo,    ONLY : vps, vps0, rhops, dvps, drhops
  use atom,            ONLY : rgrid
  use uspp_param,      ONLY : upf, oldvan
  use pseudo_base,     ONLY : compute_rhops, formfn, formfa, compute_eself
  use pseudopotential, ONLY : tpstab, vps_sp, dvps_sp
  use cp_interfaces,   ONLY : build_pstab
  use splines,         ONLY : spline
  use gvect, ONLY : gstart, gg
  use constants,       ONLY : autoev
  !
  implicit none
  logical      :: tfirst
  real(DP)    :: eself, DeltaV0
  !
  real(DP)    :: vpsum, rhopsum
  integer      :: is, ig
  REAL(DP)    :: cost1, xg

  call start_clock( 'formf' )
  !
  IF( .NOT. ALLOCATED( rgrid ) ) &
     CALL errore( ' formf ', ' rgrid not allocated ', 1 )
  IF( .NOT. ALLOCATED( upf ) ) &
     CALL errore( ' formf ', ' upf not allocated ', 1 )
  !
  ! calculation of gaussian selfinteraction
  !
  eself = compute_eself( na, zv, rcmax, nsp )

  if( tfirst .or. ( iprsta >= 4 ) )then
     WRITE( stdout, 1200 ) eself
  endif
  !
  1200 format(/,3x,'formf: eself=',f12.5)
  !
  IF( tpstab ) THEN
     !
     CALL build_pstab( )
     !
  END IF
  !
  do is = 1, nsp

     IF( tpstab ) THEN
        !
        !  Use interpolation table, with cubic spline
        !
        cost1 = 1.0d0/omega
        !
        IF( gstart == 2 ) THEN
           vps (1,is) =  vps_sp(is)%y(1) * cost1
           dvps(1,is) = dvps_sp(is)%y(1) * cost1
        END IF
        !
        DO ig = gstart, ngms
           xg = SQRT( gg(ig) ) * tpiba
           vps (ig,is) = spline(  vps_sp(is), xg ) * cost1
           dvps(ig,is) = spline( dvps_sp(is), xg ) * cost1
        END DO
        !
     ELSE

        call formfn( rgrid(is)%r, rgrid(is)%rab, &
                     upf(is)%vloc(1:rgrid(is)%mesh), zv(is), rcmax(is), gg, &
                     omega, tpiba2, rgrid(is)%mesh, ngms, oldvan(is), tpre, &
                     vps(:,is), vps0(is), dvps(:,is) )

! obsolete BHS form
! call formfa( vps(:,is), dvps(:,is), rc1(is), rc2(is), wrc1(is), wrc2(is), &
!              rcl(:,is,lloc(is)), al(:,is,lloc(is)), bl(:,is,lloc(is)),    &
!              zv(is), rcmax(is), g, omega, tpiba2, ngms, gstart, tpre )

     END IF
     !
     !     fourier transform of local pp and gaussian nuclear charge
     !
     call compute_rhops( rhops(:,is), drhops(:,is), zv(is), rcmax(is), gg, &
                         omega, tpiba2, ngms, tpre )

     if( tfirst .or. ( iprsta >= 4 ) )then
        vpsum = SUM( vps( 1:ngms, is ) )
        rhopsum = SUM( rhops( 1:ngms, is ) )
        call mp_sum( vpsum, intra_bgrp_comm )
        call mp_sum( rhopsum, intra_bgrp_comm )
        WRITE( stdout,1250) vps(1,is),rhops(1,is)
        WRITE( stdout,1300) vpsum,rhopsum
     endif
     !
  end do
  ! 
  ! ... DeltaV0 is the shift to be applied to eigenvalues
  ! ... in order to align them to other plane wave codes
  !
  DeltaV0 = 0.0_dp
  DO is = 1, nsp
     !
     ! ...  na(is)/omega is the structure factor at G=0
     !
     DeltaV0 = DeltaV0 + na(is) / omega * vps0(is)
  END DO
  !
   write(6,'("   Delta V(G=0): ",f10.6,"Ry, ",f11.6,"eV")') &
         deltaV0, deltaV0*autoev
  !
  call stop_clock( 'formf' )
  !
  1250 format(3x,'formf:     vps(g=0)=',f12.7,'     rhops(g=0)=',f12.7)
  1300 format(3x,'formf: sum_g vps(g)=',f12.7,' sum_g rhops(g)=',f12.7)
  !
  return
end subroutine formf
!
!-----------------------------------------------------------------------
SUBROUTINE newnlinit()
  !-----------------------------------------------------------------------
  !
  ! ... this routine calculates arrays beta, qq, qgb, rhocb
  ! ... and derivatives w.r.t. cell parameters dbeta
  ! ... See also comments in nlinit
  !
  use control_flags,    ONLY : tpre
  use pseudopotential,  ONLY : tpstab
  use cp_interfaces,    ONLY : interpolate_beta, interpolate_qradb, &
                               exact_beta, check_tables, exact_qradb
  !
  IMPLICIT NONE
  !
  LOGICAL :: recompute_table
  ! 
  ! ... initialization for vanderbilt species
  !
  IF( tpstab ) THEN

     recompute_table = tpre .AND. check_tables()
     !
     IF ( recompute_table ) &
        CALL errore( ' newnlinit', &
                  'interpolation tables recalculation, not implemented yet', 1 )
     !
     !     initialization that is common to all species
     !
     CALL interpolate_beta( tpre )
     !
     CALL interpolate_qradb( tpre )
     !
  ELSE
     !
     ! ... this is mainly for testing
     !
     CALL exact_beta( tpre )
     !
     CALL exact_qradb( tpre )
     !
  END IF
  !
  ! ... non-linear core-correction   ( rhocb(ig,is) )
  !
  CALL core_charge_ftr( tpre )
  !
  RETURN
  !
END SUBROUTINE newnlinit
!
!-----------------------------------------------------------------------
subroutine nlfh_x( stress, bec_bgrp, dbec, lambda )
  !-----------------------------------------------------------------------
  !
  !     contribution to the internal stress tensor due to the constraints
  !
  USE kinds,             ONLY : DP
  use uspp,              ONLY : nkb, qq
  use uspp_param,        ONLY : nh, nhm, nvb, ish
  use ions_base,         ONLY : na
  use electrons_base,    ONLY : nbspx, nbsp, nudx, nspin, nupdwn, iupdwn, ibgrp_g2l
  use cell_base,         ONLY : omega, h
  use constants,         ONLY : pi, fpi, au_gpa
  use io_global,         ONLY : stdout
  use control_flags,     ONLY : iprsta
  USE cp_main_variables, ONLY : descla, la_proc, nlam, nlax
  USE descriptors,       ONLY : nlar_ , nlac_ , ilar_ , ilac_ , nlax_
  USE mp,                ONLY : mp_sum
  USE mp_global,         ONLY : intra_bgrp_comm, inter_bgrp_comm

!
  implicit none

  REAL(DP), INTENT(INOUT) :: stress(3,3) 
  REAL(DP), INTENT(IN)    :: bec_bgrp( :, : ), dbec( :, :, :, : )
  REAL(DP), INTENT(IN)    :: lambda( :, :, : )
!
  INTEGER  :: i, j, ii, jj, inl, iv, jv, ia, is, iss, nss, istart
  INTEGER  :: jnl, ir, ic, nr, nc, nx, ibgrp_i
  REAL(DP) :: fpre(3,3), TT, T1, T2
  !
  REAL(DP), ALLOCATABLE :: tmpbec(:,:), tmpdh(:,:), temp(:,:), bec(:,:,:)
  !
  ALLOCATE( bec( nkb, nlax, nspin ) )
  !
  IF( la_proc ) THEN
     DO iss = 1, nspin
        nss = nupdwn( iss )
        istart = iupdwn( iss )
        ic = descla( ilac_ , iss )
        nc = descla( nlac_ , iss )
        DO i=1,nc
           ibgrp_i = ibgrp_g2l( i+istart-1+ic-1 )
           IF( ibgrp_i > 0 ) THEN
              bec( :, i, iss ) = bec_bgrp( :, ibgrp_i )
           ELSE
              bec( :, i, iss ) = 0.0d0
           END IF
        END DO
     END DO
  ELSE
     bec   = 0.0d0
  END IF

  CALL mp_sum( bec, inter_bgrp_comm )
  !
  IF( la_proc ) THEN
     nx=descla( nlax_ , 1 ) 
     IF( nspin == 2 ) nx = MAX( nx , descla( nlax_ , 2 ) )
     ALLOCATE ( tmpbec(nhm,nx), tmpdh(nx,nhm), temp(nx,nx) )
  END IF
  !
  fpre = 0.d0
  !
  do ii=1,3

     do jj=1,3

        do is=1,nvb

           do ia=1,na(is)

              do iss = 1, nspin
                 !
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 !
                 IF( la_proc ) THEN

                    nr = descla( nlar_ , iss )
                    nc = descla( nlac_ , iss )
                    ir = descla( ilar_ , iss )
                    ic = descla( ilac_ , iss )

                    tmpbec = 0.d0
                    tmpdh  = 0.d0
!
                    do iv=1,nh(is)
                       do jv=1,nh(is)
                          inl=ish(is)+(jv-1)*na(is)+ia
                          if(abs(qq(iv,jv,is)).gt.1.e-5) then
                             do i = 1, nc
                                !tmpbec(iv,i) = tmpbec(iv,i) +  qq(iv,jv,is) * bec(inl, i + istart - 1 + ic - 1 )
                                tmpbec(iv,i) = tmpbec(iv,i) +  qq(iv,jv,is) * bec( inl, i, iss  )
                             end do
                          endif
                       end do
                    end do

                    do iv=1,nh(is)
                       inl=ish(is)+(iv-1)*na(is)+ia
                       do i = 1, nr
                          tmpdh(i,iv) = dbec( inl, i + (iss-1)*nlax, ii, jj )
                       end do
                    end do

                    if(nh(is).gt.0)then

                       CALL dgemm &
                       ( 'N', 'N', nr, nc, nh(is), 1.0d0, tmpdh, nx, tmpbec, nhm, 0.0d0, temp, nx )

                       do j = 1, nc
                          do i = 1, nr
                             fpre(ii,jj) = fpre(ii,jj) + 2D0 * temp( i, j ) * lambda(i,j,iss)
                          end do
                       end do
                    endif

                 END IF
                 !
              end do
              !
           end do
           !
        end do
        !
     end do
     !
  end do

  CALL mp_sum( fpre, intra_bgrp_comm )

  do i=1,3
     do j=1,3
        stress(i,j)=stress(i,j)+ &
                    (fpre(i,1)*h(j,1)+fpre(i,2)*h(j,2)+fpre(i,3)*h(j,3))/omega
     enddo
  enddo

  IF( la_proc ) THEN
     DEALLOCATE ( tmpbec, tmpdh, temp )
  END IF

  DEALLOCATE( bec )


  IF( iprsta > 2 ) THEN
     WRITE( stdout,*) 
     WRITE( stdout,*) "constraints contribution to stress"
     WRITE( stdout,5555) ((-fpre(i,j),j=1,3),i=1,3)
     fpre = MATMUL( fpre, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
     WRITE( stdout,5555) ((fpre(i,j),j=1,3),i=1,3)
     WRITE( stdout,*) 
  END IF
!

5555  FORMAT(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

  return
end subroutine nlfh_x


!-----------------------------------------------------------------------
subroutine nlinit
  !-----------------------------------------------------------------------
  !
  !     this routine allocates and initalizes arrays beta, qq, qgb,
  !     rhocb, and derivatives w.r.t. cell parameters dbeta
  !
  !       beta(ig,l,is) = 4pi/sqrt(omega) y^r(l,q^)
  !                               int_0^inf dr r^2 j_l(qr) betar(l,is,r)
  !
  !       Note that beta(g)_lm,is = (-i)^l*beta(ig,l,is) (?)
  !
  !       qq_ij=int_0^r q_ij(r)=omega*qg(g=0)
  !
  !     beta and qradb are first calculated on a fixed linear grid in |G|
  !     (betax, qradx) then calculated on the box grid by interpolation
  !     (this is done in routine newnlinit)
  !
      use kinds,           ONLY : dp
      use control_flags,   ONLY : iprint, tpre
      use io_global,       ONLY : stdout, ionode
      use gvecw,           ONLY : ngw
      use core,            ONLY : rhocb, nlcc_any, allocate_core
      use constants,       ONLY : pi, fpi
      use ions_base,       ONLY : na, nsp
      use uspp,            ONLY : aainit, beta, qq, dvan, nhtol, nhtolm, indv, dbeta
      use uspp_param,      ONLY : upf, lmaxq, nbetam, lmaxkb, nhm, nh, ish, nvb
      use atom,            ONLY : rgrid
      use qgb_mod,         ONLY : qgb, dqgb
      use gvecb,           ONLY : ngb
      use gvect,           ONLY : ngm
      use betax,           ONLY : qradx, dqradx, refg, betagx, mmx, dbetagx
      use cp_interfaces,   ONLY : pseudopotential_indexes, compute_dvan, &
                                  compute_betagx, compute_qradx
      USE grid_dimensions, ONLY : nrxx

!
      implicit none
!
      integer  is, il, l, ir, iv, jv, lm, ind, ltmp, i0
      real(dp), allocatable:: fint(:), jl(:),  jltmp(:), djl(:),    &
     &              dfint(:)
      real(dp) xg, xrg, fac


      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( //, &
                3X,'Pseudopotentials initialization',/, &
                3X,'-------------------------------' )
      END IF

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' nlinit ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' nlinit ', ' upf not allocated ', 1 )
      !
      !   initialize indexes
      !
      CALL pseudopotential_indexes( )
      !
      !   initialize array ap
      !
      call aainit( lmaxkb + 1 )
      !
      CALL allocate_core( nrxx, ngm, ngb, nsp )
      !
      !
      allocate( beta( ngw, nhm, nsp ) )
      allocate( qgb( ngb, nhm*(nhm+1)/2, nsp ) )
      allocate( qq( nhm, nhm, nsp ) )
      qq  (:,:,:) =0.d0
      IF (tpre) THEN
         allocate( dqgb( ngb, nhm*(nhm+1)/2, nsp, 3, 3 ) )
         allocate( dbeta( ngw, nhm, nsp, 3, 3 ) )
      END IF
      !
      !     initialization for vanderbilt species
      !
      CALL compute_qradx( tpre )
      !    
      !     initialization that is common to all species
      !   
      WRITE( stdout, fmt="(//,3X,'Common initialization' )" )

      do is = 1, nsp
         WRITE( stdout, fmt="(/,3X,'Specie: ',I5)" ) is
         !     fac converts ry to hartree
         fac=0.5d0
         do iv = 1, nh(is)
            WRITE( stdout,901) iv, indv(iv,is), nhtol(iv,is)
         end do
 901     format(2x,i2,'  indv= ',i2,'   ang. mom= ',i2)
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    dion '
         do iv = 1, upf(is)%nbeta
            WRITE( stdout,'(8f9.4)') ( fac*upf(is)%dion(iv,jv), jv = 1, upf(is)%nbeta )
         end do
         !
      end do
      !
      !   calculation of array  betagx(ig,iv,is)
      !
      call compute_betagx( tpre )
      !
      !   calculate array  dvan(iv,jv,is)
      !
      call compute_dvan()
      !
      ! newnlinit stores qgb and qq, calculates arrays  beta  rhocb
      ! and derivatives wrt cell dbeta
      !
      call newnlinit()

      return
end subroutine nlinit

!-------------------------------------------------------------------------
subroutine qvan2b(ngy,iv,jv,is,ylm,qg,qradb)
  !--------------------------------------------------------------------------
  !
  !     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use gvecb,         ONLY : ngb
  use uspp_param,    ONLY : lmaxq, nbetam
  use ions_base,     ONLY : nsp
! 
  implicit none
  !
  integer,     intent(in)  :: ngy, iv, jv, is
  real(DP),    intent(in)  :: ylm( ngb, lmaxq*lmaxq )
  real(DP),    intent(in)  :: qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp )
  complex(DP), intent(out) :: qg( ngb )
!
  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  ! ijvs is the packed index for (ivs,jvs)
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2b ', ' wrong dimensions', MAX(ivl,jvl))
  !
  qg(:) = (0.d0, 0.d0)
  !
  !     lpx = max number of allowed y_lm
  !     lp  = composite lm to indentify them
  !
  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' qvan2b ',' lp out of bounds ',lp)
     !
     !     extraction of angular momentum l from lp:  
     !     l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig=(0.d0,-1.d0)**(l-1)
     sig=sig*ap(lp,ivl,jvl)
     do ig=1,ngy
        qg(ig)=qg(ig)+sig*ylm(ig,lp)*qradb(ig,ijvs,l,is)
     end do
  end do

  return
end subroutine qvan2b

!-------------------------------------------------------------------------
subroutine dqvan2b(ngy,iv,jv,is,ylm,dylm,dqg,dqrad,qradb)
  !--------------------------------------------------------------------------
  !
  !     dq(i,j) derivatives wrt to h(i,j) of q(g,l,k) calculated in qvan2b
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use gvecb,         ONLY : ngb
  use uspp_param,    ONLY : lmaxq, nbetam
  use ions_base,     ONLY : nsp

  implicit none

  integer,     intent(in)  :: ngy, iv, jv, is
  REAL(DP),    INTENT(IN)  :: ylm( ngb, lmaxq*lmaxq ), dylm( ngb, lmaxq*lmaxq, 3, 3 )
  complex(DP), intent(out) :: dqg( ngb, 3, 3 )
  REAL(DP),    INTENT(IN)  :: dqrad( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp, 3, 3 )
  real(DP),    intent(in)  :: qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp )

  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig, z1, z2, zfac
  !
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  !
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  !
  ! ijvs is the packed index for (ivs,jvs)
  !
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  !
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  dqg(:,:,:) = (0.d0, 0.d0)

  !  lpx = max number of allowed y_lm
  !  lp  = composite lm to indentify them

  z1 = 0.0d0
  z2 = 0.0d0
  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' dqvan2b ',' lp out of bounds ',lp)

     !  extraction of angular momentum l from lp:  
     !  l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig = (0.0d0,-1.0d0)**(l-1)
     !
     sig = sig * ap( lp, ivl, jvl ) 
     !
     do ij=1,3
        do ii=1,3
           do ig=1,ngy
              zfac = ylm(ig,lp) * dqrad(ig,ijvs,l,is,ii,ij)
              zfac = zfac - dylm(ig,lp,ii,ij) * qradb(ig,ijvs,l,is)
              dqg(ig,ii,ij) = dqg(ig,ii,ij) +  sig * zfac
           end do
        end do
     end do
  end do
  !
  ! WRITE(6,*) 'DEBUG dqvan2b: ', z1, z2
  !
  return
end subroutine dqvan2b

!-----------------------------------------------------------------------
subroutine dylmr2_( nylm, ngy, g, gg, ainv, dylm )
  !-----------------------------------------------------------------------
  !
  ! temporary CP interface for PW routine dylmr2
  ! dylmr2  calculates d Y_{lm} /d G_ipol
  ! dylmr2_ calculates G_ipol \sum_k h^(-1)(jpol,k) (dY_{lm} /dG_k)
  !
  USE kinds, ONLY: DP

  implicit none
  !
  integer,   intent(IN)  :: nylm, ngy
  real(DP), intent(IN)  :: g (3, ngy), gg (ngy), ainv(3,3)
  real(DP), intent(OUT) :: dylm (ngy, nylm, 3, 3)
  !
  integer :: ipol, jpol, lm, ig
  real(DP), allocatable :: dylmaux (:,:,:)
  !
  allocate ( dylmaux(ngy,nylm,3) )
  !
  dylmaux(:,:,:) = 0.d0
  !
  do ipol =1,3
     call dylmr2 (nylm, ngy, g, gg, dylmaux(1,1,ipol), ipol)
  enddo
  !
  do ipol =1,3
     do jpol =1,3
        do lm=1,nylm
           do ig = 1, ngy
              dylm (ig,lm,ipol,jpol) = (dylmaux(ig,lm,1) * ainv(jpol,1) + & 
                                        dylmaux(ig,lm,2) * ainv(jpol,2) + &
                                        dylmaux(ig,lm,3) * ainv(jpol,3) ) &
                                       * g(ipol,ig)
           end do
        end do
     end do
  end do
  !
  deallocate ( dylmaux )
  !
  return
  !
end subroutine dylmr2_


SUBROUTINE print_lambda_x( lambda, n, nshow, ccc, iunit )
    USE kinds, ONLY : DP
    USE io_global,         ONLY: stdout, ionode
    USE cp_main_variables, ONLY: collect_lambda, descla
    USE electrons_base,    ONLY: nudx
    IMPLICIT NONE
    real(DP), intent(in) :: lambda(:,:,:), ccc
    integer, intent(in) :: n, nshow
    integer, intent(in), optional :: iunit
    !
    integer :: nnn, j, un, i, is
    real(DP), allocatable :: lambda_repl(:,:)
    if( present( iunit ) ) then
      un = iunit
    else
      un = stdout
    end if
    nnn = min( nudx, nshow )
    ALLOCATE( lambda_repl( nudx, nudx ) )
    IF( ionode ) WRITE( un,*)
    DO is = 1, SIZE( lambda, 3 )
       CALL collect_lambda( lambda_repl, lambda(:,:,is), descla(:,is) )
       IF( ionode ) THEN
          WRITE( un,3370) '    lambda   nudx, spin = ', nudx, is
          IF( nnn < n ) WRITE( un,3370) '    print only first ', nnn
          DO i=1,nnn
             WRITE( un,3380) (lambda_repl(i,j)*ccc,j=1,nnn)
          END DO
       END IF
    END DO
    DEALLOCATE( lambda_repl )
3370   FORMAT(26x,a,2i4)
3380   FORMAT(9f8.4)
    RETURN
END SUBROUTINE print_lambda_x
!-----------------------------------------------------------------------
FUNCTION n_atom_wfc_x( )
!----------------------------------------------------------------------------
  !
  ! ... Find max number of bands needed
  !
  USE ions_base,        ONLY : na, nsp
  USE kinds,            ONLY : DP
  USE uspp_param,       ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER  :: n_atom_wfc_x
  INTEGER  :: is, n
  !
  n_atom_wfc_x = 0
  !
  DO is = 1, nsp
     !
     DO n = 1, upf(is)%nwfc
        !
        IF ( upf(is)%oc(n) >= 0.D0 ) THEN
           !
           n_atom_wfc_x = n_atom_wfc_x + na(is) * ( 2*upf(is)%lchi(n) + 1 )
           !
        END IF
        !
     END DO
     !
  END DO
  !
  RETURN
END FUNCTION

!
!
!
!-----------------------------------------------------------------------
      SUBROUTINE denlcc( nnr, nspin, vxcr, sfac, drhocg, dcc )
!-----------------------------------------------------------------------
!
! derivative of non linear core correction exchange energy wrt cell 
! parameters h 
! Output in dcc
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE gvect, ONLY: gstart, g, gg
      USE gvecs,              ONLY: ngms
      USE gvect,              ONLY: ngm, nl
      USE cell_base,          ONLY: omega, ainv, tpiba2
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE uspp_param,         ONLY: upf
      USE grid_dimensions,    ONLY: nr1, nr2, nr3
      USE fft_interfaces,     ONLY: fwfft
      USE fft_base,           ONLY: dfftp

      IMPLICIT NONE

      ! input

      INTEGER, INTENT(IN)   :: nnr, nspin
      REAL(DP)              :: vxcr( nnr, nspin )
      COMPLEX(DP)           :: sfac( ngms, nsp )
      REAL(DP)              :: drhocg( ngm, nsp )

      ! output

      REAL(DP), INTENT(OUT) ::  dcc(3,3)

      ! local

      INTEGER     :: i, j, ig, is
      COMPLEX(DP) :: srhoc
      REAL(DP)    :: vxcc
      !
      COMPLEX(DP), ALLOCATABLE :: vxc( : )
!
      dcc = 0.0d0
      !
      ALLOCATE( vxc( nnr ) )
      !
      vxc(:) = vxcr(:,1)
      !
      IF( nspin > 1 ) vxc(:) = vxc(:) + vxcr(:,2)
      !
      CALL fwfft( 'Dense', vxc, dfftp )
      !
      DO i=1,3
         DO j=1,3
            DO ig = gstart, ngms
               srhoc = 0.0d0
               DO is = 1, nsp
                 IF( upf(is)%nlcc ) srhoc = srhoc + sfac( ig, is ) * drhocg( ig, is )
               ENDDO
               vxcc = DBLE( CONJG( vxc( nl( ig ) ) ) * srhoc ) / SQRT( gg(ig) * tpiba2 )
               dcc(i,j) = dcc(i,j) + vxcc * &
     &                      2.d0 * tpiba2 * g(i,ig) *                  &
     &                    (g(1,ig)*ainv(j,1) +                         &
     &                     g(2,ig)*ainv(j,2) +                         &
     &                     g(3,ig)*ainv(j,3) )
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( vxc )

      dcc = dcc * omega

      CALL mp_sum( dcc( 1:3, 1:3 ), intra_bgrp_comm )

      RETURN
      END SUBROUTINE denlcc



!-----------------------------------------------------------------------
      SUBROUTINE dotcsc( eigr, cp, ngw, n )
!-----------------------------------------------------------------------
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: na, nsp, nat
      USE io_global,          ONLY: stdout
      USE gvect, ONLY: gstart
      USE uspp,               ONLY: nkb, qq
      USE uspp_param,         ONLY: nh, ish, nvb
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm, nbgrp
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ngw, n
      COMPLEX(DP) ::  eigr(ngw,nat), cp(ngw,n)
! local variables
      REAL(DP) rsum, csc(n) ! automatic array
      COMPLEX(DP) temp(ngw) ! automatic array
 
      REAL(DP), ALLOCATABLE::  becp(:,:)
      INTEGER i,kmax,nnn,k,ig,is,ia,iv,jv,inl,jnl
!
      IF( nbgrp > 1 ) &
         CALL errore( ' dotcsc ', ' parallelization over bands not yet implemented ', 1 )
!
      ALLOCATE(becp(nkb,n))
!
!     < beta | phi > is real. only the i lowest:
!
      nnn = MIN( 12, n )

      DO i = nnn, 1, -1
         kmax = i
         CALL nlsm1(i,1,nvb,eigr,cp,becp)
!
         DO k=1,kmax
            DO ig=1,ngw
               temp(ig)=CONJG(cp(ig,k))*cp(ig,i)
            END DO
            csc(k)=2.d0*DBLE(SUM(temp))
            IF (gstart == 2) csc(k)=csc(k)-DBLE(temp(1))
         END DO

         CALL mp_sum( csc( 1:kmax ), intra_bgrp_comm )

         DO k=1,kmax
            rsum=0.d0
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     DO ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        rsum = rsum +                                    &
     &                   qq(iv,jv,is)*becp(inl,i)*becp(jnl,k)
                     END DO
                  END DO
               END DO
            END DO
            csc(k)=csc(k)+rsum
         END DO
!
         WRITE( stdout,'("dotcsc =",12f18.15)') (csc(k),k=1,i)
!
      END DO
      WRITE( stdout,*)
!
      DEALLOCATE(becp)
!
      RETURN
      END SUBROUTINE dotcsc


!
!-----------------------------------------------------------------------
   FUNCTION enkin_x( c, f, n )
!-----------------------------------------------------------------------
      !
      ! calculation of kinetic energy term
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: pi, fpi
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE gvecw,              ONLY: ggp
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE cell_base,          ONLY: tpiba2

      IMPLICIT NONE

      REAL(DP)                :: enkin_x

      ! input

      INTEGER,     INTENT(IN) :: n
      COMPLEX(DP), INTENT(IN) :: c( :, : )
      REAL(DP),    INTENT(IN) :: f( : )
      !
      ! local

      INTEGER  :: ig, i
      REAL(DP) :: sk(n)  ! automatic array
      !
      DO i=1,n
         sk(i)=0.0d0
         DO ig=gstart,ngw
            sk(i)=sk(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*ggp(ig)
         END DO
      END DO

      CALL mp_sum( sk(1:n), intra_bgrp_comm )

      enkin_x=0.0d0
      DO i=1,n
         enkin_x=enkin_x+f(i)*sk(i)
      END DO

      ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
      ! ... multiplicative factor (2 pi/alat)**2 is required

      enkin_x = enkin_x * tpiba2
!
      RETURN
   END FUNCTION enkin_x
!            
!
!-----------------------------------------------------------------------
      SUBROUTINE initbox ( tau0, taub, irb, ainv, at, alat )
!-----------------------------------------------------------------------
!
!     sets the indexes irb and positions taub for the small boxes 
!     around atoms
!
      USE kinds,                    ONLY: DP
      USE ions_base,                ONLY: nsp, na, nat
      USE grid_dimensions,          ONLY: nr1, nr2, nr3
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
      USE control_flags,            ONLY: iprsta
      USE io_global,                ONLY: stdout
      USE mp_global,                ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm
      USE fft_base,                 ONLY: dfftb, dfftp, fft_dlay_descriptor
      USE fft_types,                ONLY: fft_box_set
      USE uspp_param,               ONLY: nvb

      IMPLICIT NONE
! input
      REAL(DP), INTENT(in)  :: tau0(3,nat)
! output
      INTEGER,  INTENT(out) :: irb(3,nat)
      REAL(DP), INTENT(out) :: taub(3,nat)
! input
      REAL(DP), INTENT(in)  :: ainv(3,3)
      REAL(DP), INTENT(in)  :: at(3,3), alat
! local
      REAL(DP) :: x(3), xmod
      INTEGER  :: nr(3), nrb(3), xint, is, ia, i, isa
!
      IF ( nr1b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 1)
      IF ( nr2b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 2)
      IF ( nr3b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 3)

      nr (1)=nr1
      nr (2)=nr2
      nr (3)=nr3
      nrb(1)=nr1b
      nrb(2)=nr2b
      nrb(3)=nr3b
!
      isa = 0
      DO is=1,nsp
         DO ia=1,na(is)
           isa = isa + 1
!
            DO i=1,3
!
! bring atomic positions to crystal axis
!
               x(i) = ainv(i,1)*tau0(1,isa) +                         &
     &                ainv(i,2)*tau0(2,isa) +                         &
     &                ainv(i,3)*tau0(3,isa)
!
! bring x in the range between 0 and 1
!
               x(i) = MOD(x(i),1.d0)
               IF (x(i).LT.0.d0) x(i)=x(i)+1.d0
!
! case of nrb(i) even
!
               IF (MOD(nrb(i),2).EQ.0) THEN
!
! find irb = index of the grid point at the corner of the small box
!           (the indices of the small box run from irb to irb+nrb-1)
!
                  xint=INT(x(i)*nr(i))
                  irb (i,isa)=xint+1-nrb(i)/2+1
                  IF(irb(i,isa).LT.1) irb(i,isa)=irb(i,isa)+nr(i)
!
! x(i) are the atomic positions in crystal coordinates, where the
! "crystal lattice" is the small box lattice and the origin is at
! the corner of the small box. Used to calculate phases exp(iG*taub)
!
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+nrb(i)/2-1)/nr(i)
               ELSE
!
! case of nrb(i) odd - see above for comments
!
                  xint=NINT(x(i)*nr(i))
                  irb (i,isa)=xint+1-(nrb(i)-1)/2
                  IF(irb(i,isa).LT.1) irb(i,isa)=irb(i,isa)+nr(i)
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+(nrb(i)-1)/2)/nr(i)
               END IF
            END DO
!
! bring back taub in cartesian coordinates
!
            DO i=1,3
               taub(i,isa)=(x(1)*at(i,1) + x(2)*at(i,2) + x(3)*at(i,3))*alat
            END DO
         END DO
      END DO

      ! initialize FFT descriptor

      CALL fft_box_set( dfftb, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, &
                        nat, irb, dfftp%npp, dfftp%ipp )

      IF( iprsta > 2 ) THEN
           isa = 1
           DO is=1,nsp
              WRITE( stdout, '( /, 2x, "species= ", i2 )' ) is
              DO ia=1,na(is)
                 WRITE( stdout,2000) ia, (irb(i,isa),i=1,3)
2000             FORMAT(2x, 'atom= ', i3, ' irb1= ', i3, ' irb2= ', i3, ' irb3= ', i3)
                 isa = isa + 1
               END DO
            END DO
      ENDIF

#ifdef __PARA
      ! 
      ! for processor that do not call fft on the box
      ! artificially start the clock
      ! 
      CALL start_clock( 'fftb' )
      CALL stop_clock( 'fftb' )
      !
#endif
!
      RETURN
   END SUBROUTINE initbox

!-------------------------------------------------------------------------
      SUBROUTINE nlfl_bgrp( bec_bgrp, becdr_bgrp, lambda, fion )
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
! 
!
      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE ions_base,         ONLY: na, nsp, nat
      USE uspp,              ONLY: nhsa=>nkb, qq
      USE uspp_param,        ONLY: nhm, nh, ish, nvb
      USE electrons_base,    ONLY: nspin, iupdwn, nupdwn, nbspx_bgrp, ibgrp_g2l, i2gupdwn_bgrp, nbspx, &
                                   iupdwn_bgrp, nupdwn_bgrp
      USE constants,         ONLY: pi, fpi
      USE cp_main_variables, ONLY: nlam, nlax, descla, la_proc
      USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , la_myr_ , la_myc_ 
      USE mp,                ONLY: mp_sum
      USE mp_global,         ONLY: intra_bgrp_comm, inter_bgrp_comm
!
      IMPLICIT NONE
      REAL(DP) bec_bgrp(nhsa,nbspx_bgrp), becdr_bgrp(nhsa,nbspx_bgrp,3), lambda(nlam,nlam,nspin)
      REAL(DP) fion(3,nat)
!
      INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc, ibgrp_i
      INTEGER :: n1, n2, m1, m2
      REAL(DP), ALLOCATABLE :: temp(:,:), tmpbec(:,:),tmpdr(:,:) 
      REAL(DP), ALLOCATABLE :: fion_tmp(:,:)
      REAL(DP), ALLOCATABLE :: bec(:,:,:)
      REAL(DP), ALLOCATABLE :: becdr(:,:,:,:)
      REAL(DP), ALLOCATABLE :: bec_g(:,:)
      REAL(DP), ALLOCATABLE :: becdr_g(:,:,:)
      !
      CALL start_clock( 'nlfl' )
      !
      ALLOCATE( fion_tmp( 3, nat ) )
      !
      fion_tmp = 0.0d0
      !
      ALLOCATE( temp( nlax, nlax ), tmpbec( nhm, nlax ), tmpdr( nlax, nhm ) )
      ALLOCATE( bec( nhsa, nlax, nspin ), becdr( nhsa, nlax, nspin, 3 ) )

      ! redistribute bec, becdr according to the ortho subgroup
      ! this is required because they are combined with "lambda" matrixes
      
      IF( la_proc ) THEN
         DO iss = 1, nspin
            nss = nupdwn( iss )
            istart = iupdwn( iss )
            ic = descla( ilac_ , iss )
            nc = descla( nlac_ , iss )
            DO i=1,nc
               ibgrp_i = ibgrp_g2l( i+istart-1+ic-1 )
               IF( ibgrp_i > 0 ) THEN
                  bec( :, i, iss ) = bec_bgrp( :, ibgrp_i )
               ELSE
                  bec( :, i, iss ) = 0.0d0
               END IF
            END DO
            ir = descla( ilar_ , iss )
            nr = descla( nlar_ , iss )
            DO i=1,nr
               ibgrp_i = ibgrp_g2l( i+istart-1+ir-1 )
               IF( ibgrp_i > 0 ) THEN
                  becdr(:,i,iss,1) = becdr_bgrp( :, ibgrp_i, 1 )
                  becdr(:,i,iss,2) = becdr_bgrp( :, ibgrp_i, 2 )
                  becdr(:,i,iss,3) = becdr_bgrp( :, ibgrp_i, 3 )
               ELSE
                  becdr(:,i,iss,1) = 0.0d0
                  becdr(:,i,iss,2) = 0.0d0
                  becdr(:,i,iss,3) = 0.0d0
               END IF
            END DO
         END DO
      ELSE
         bec   = 0.0d0
         becdr = 0.0d0
      END IF

      CALL mp_sum( bec, inter_bgrp_comm )
      CALL mp_sum( becdr, inter_bgrp_comm )
      !
      DO k=1,3
         isa = 0
         DO is=1,nvb
            DO ia=1,na(is)
               isa = isa + 1
               !
               DO iss = 1, nspin
                  !
                  nss = nupdwn( iss )
                  istart = iupdwn( iss )
                  !
                  tmpbec = 0.d0
                  tmpdr  = 0.d0
                  !
                  IF( la_proc ) THEN
                     ! tmpbec distributed by columns
                     ic = descla( ilac_ , iss )
                     nc = descla( nlac_ , iss )
                     DO iv=1,nh(is)
                        DO jv=1,nh(is)
                           inl=ish(is)+(jv-1)*na(is)+ia
                           IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                              DO i=1,nc
                                 tmpbec(iv,i)=tmpbec(iv,i) + qq(iv,jv,is)*bec(inl,i,iss)
                              END DO
                           ENDIF
                        END DO
                     END DO
                     ! tmpdr distributed by rows
                     ir = descla( ilar_ , iss )
                     nr = descla( nlar_ , iss )
                     DO iv=1,nh(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        DO i=1,nr
                           tmpdr(i,iv) = becdr( inl, i, iss, k )
                        END DO
                     END DO
                  END IF
                  !
                  IF(nh(is).GT.0)THEN
                     !
                     IF( la_proc ) THEN
                        ir = descla( ilar_ , iss )
                        ic = descla( ilac_ , iss )
                        nr = descla( nlar_ , iss )
                        nc = descla( nlac_ , iss )
                        CALL dgemm( 'N', 'N', nr, nc, nh(is), 1.0d0, tmpdr, nlax, tmpbec, nhm, 0.0d0, temp, nlax )
                        DO j = 1, nc
                           DO i = 1, nr
                              fion_tmp(k,isa) = fion_tmp(k,isa) + 2D0 * temp( i, j ) * lambda( i, j, iss )
                           END DO
                        END DO
                     END IF
!
                  ENDIF

               END DO
!
            END DO
         END DO
      END DO
      !
      DEALLOCATE( bec, becdr )
      DEALLOCATE( temp, tmpbec, tmpdr )
      !
      CALL mp_sum( fion_tmp, intra_bgrp_comm )
      !
      fion = fion + fion_tmp
      !
      DEALLOCATE( fion_tmp )
      !
      CALL stop_clock( 'nlfl' )
      !
      RETURN

      END SUBROUTINE nlfl_bgrp


!
!-----------------------------------------------------------------------
      SUBROUTINE pbc(rin,a1,a2,a3,ainv,rout)
!-----------------------------------------------------------------------
!
!     brings atoms inside the unit cell
!
      USE kinds,  ONLY: DP

      IMPLICIT NONE
! input
      REAL(DP) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
! output
      REAL(DP) rout(3)
! local
      REAL(DP) x,y,z
!
! bring atomic positions to crystal axis
!
      x = ainv(1,1)*rin(1)+ainv(1,2)*rin(2)+ainv(1,3)*rin(3)
      y = ainv(2,1)*rin(1)+ainv(2,2)*rin(2)+ainv(2,3)*rin(3)
      z = ainv(3,1)*rin(1)+ainv(3,2)*rin(2)+ainv(3,3)*rin(3)
!
! bring x,y,z in the range between -0.5 and 0.5
!
      x = x - NINT(x)
      y = y - NINT(y)
      z = z - NINT(z)
!
! bring atomic positions back in cartesian axis
!
      rout(1) = x*a1(1)+y*a2(1)+z*a3(1)
      rout(2) = x*a1(2)+y*a2(2)+z*a3(2)
      rout(3) = x*a1(3)+y*a2(3)+z*a3(3)
!
      RETURN
      END SUBROUTINE pbc

!
!-------------------------------------------------------------------------
      SUBROUTINE prefor(eigr,betae)
!-----------------------------------------------------------------------
!
!     input :        eigr =  e^-ig.r_i
!     output:        betae_i,i(g) = (-i)**l beta_i,i(g) e^-ig.r_i 
!
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nsp, na, nat
      USE gvecw,      ONLY : ngw
      USE uspp,       ONLY : nkb, beta, nhtol
      USE uspp_param, ONLY : nh, ish
!
      IMPLICIT NONE
      COMPLEX(DP) :: eigr( ngw, nat )
      COMPLEX(DP) :: betae( ngw, nkb )
!
      INTEGER     :: is, iv, ia, inl, ig, isa
      COMPLEX(DP) :: ci
!
      CALL start_clock( 'prefor' )
      isa = 0
      DO is=1,nsp
         DO iv=1,nh(is)
            ci=(0.0d0,-1.0d0)**nhtol(iv,is)
            DO ia=1,na(is)
               inl=ish(is)+(iv-1)*na(is)+ia
               DO ig=1,ngw
                  betae(ig,inl)=ci*beta(ig,iv,is)*eigr(ig,ia+isa)
               END DO
            END DO
         END DO
         isa = isa + na(is)
      END DO
      CALL stop_clock( 'prefor' )
!
      RETURN
      END SUBROUTINE prefor
!
!
!-------------------------------------------------------------------------
      SUBROUTINE s_wfc(n_atomic_wfc1,becwfc,betae,wfc,swfc) !#@@@ Changed n_atomic_wfc to n_atomic_wfc1
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
      USE kinds, ONLY: DP
      USE ions_base, ONLY: na
      USE uspp, ONLY: nhsa => nkb, nhsavb=>nkbus, qq
      USE uspp_param, ONLY: nh, nvb, ish
      USE gvecw, ONLY: ngw
      USE constants, ONLY: pi, fpi
      IMPLICIT NONE
! input
      INTEGER, INTENT(in)         :: n_atomic_wfc1
      COMPLEX(DP), INTENT(in) :: betae(ngw,nhsa),                   &
     &                               wfc(ngw,n_atomic_wfc1)
      REAL(DP), INTENT(in)    :: becwfc(nhsa,n_atomic_wfc1)
! output
      COMPLEX(DP), INTENT(out):: swfc(ngw,n_atomic_wfc1)
! local
      INTEGER is, iv, jv, ia, inl, jnl, i
      REAL(DP) qtemp(nhsavb,n_atomic_wfc1)
!
      swfc = wfc
!
      IF (nvb.GT.0) THEN
         qtemp=0.d0
         DO is=1,nvb
            DO iv=1,nh(is)
               DO jv=1,nh(is)
                  IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                     DO ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        DO i=1,n_atomic_wfc1
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*becwfc(jnl,i)
                        END DO
                     END DO
                  ENDIF
               END DO
            END DO
         END DO
!
         CALL dgemm &
              ('N','N',2*ngw,n_atomic_wfc1,nhsavb,1.0d0,betae,2*ngw,&
               qtemp,nhsavb,1.0d0,swfc,2*ngw)
!
      END IF
!
!      swfc=swfc+wfc
!
      RETURN
      END SUBROUTINE s_wfc


!-----------------------------------------------------------------------
      subroutine ldaU_init
!-----------------------------------------------------------------------
!
      USE constants,        ONLY: autoev
      use ldaU_cp,          ONLY: n_atomic_wfc, atomwfc,lda_plus_u, Hubbard_U
      use ldaU_cp,          ONLY: Hubbard_lmax, Hubbard_l, ns, vupsi
      use input_parameters, ONLY: lda_plus_u_ => lda_plus_u
      use input_parameters, ONLY: Hubbard_U_ => Hubbard_U
      use ions_base,        only: na, nsp, nat, atm
      use gvecw,            only: ngw
      use electrons_base,   only: nspin, nx => nbspx
      USE uspp_param,       ONLY: upf
      USE step_penalty,     ONLY: vpsi_pen, A_pen, sigma_pen, alpha_pen, step_pen
      use input_parameters, ONLY: step_pen_ => step_pen
      use input_parameters, ONLY: A_pen_ => A_pen
      use input_parameters, ONLY: sigma_pen_ => sigma_pen
      use input_parameters, ONLY: alpha_pen_ => alpha_pen
      !
      implicit none
      integer is, nb, l
      integer, external :: set_Hubbard_l

! allocate vupsi

      lda_plus_u = lda_plus_u_
      IF ( .NOT.lda_plus_u ) RETURN
      allocate(vupsi(ngw,nx))

      vupsi=(0.0d0,0.0d0)
      allocate(vpsi_pen(ngw,nx)) ! step_constraint 
      vpsi_pen=(0.0d0,0.0d0)
      n_atomic_wfc=0

      do is=1,nsp
         !
         Hubbard_U( is ) = Hubbard_U_( is )/autoev
         !
         do nb = 1,upf(is)%nwfc
            l = upf(is)%lchi(nb)
            n_atomic_wfc = n_atomic_wfc + (2*l+1)*na(is)
         end do
         !
      end do
!
      allocate(atomwfc(ngw,n_atomic_wfc))

      if (lda_plus_u) then
         Hubbard_lmax = -1
         do is=1,nsp
            if (Hubbard_U(is).ne.0.d0) then 
!                Hubbard_l(is)=2
               Hubbard_l(is) = set_Hubbard_l( atm(is) )
                Hubbard_lmax = max(Hubbard_lmax,Hubbard_l(is))
               write (6,*) ' HUBBARD L FOR TYPE ',atm(is),' IS ',&
     &                       Hubbard_l(is)
            end if
         end do
         write (6,*) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
         if (Hubbard_lmax.eq.-1) call errore                            &
     &        ('setup','lda_plus_u calculation but Hubbard_l not set',1)
      end if
      l = 2 * Hubbard_lmax + 1
      allocate(ns(nat,nspin,l,l))
      step_pen=step_pen_
      A_pen=A_pen_
      sigma_pen=sigma_pen_
      alpha_pen=alpha_pen_
      return
      end subroutine ldaU_init
!
!-----------------------------------------------------------------------
integer function set_Hubbard_l(psd) result (hubbard_l)
!-----------------------------------------------------------------------
!
implicit none
character*3 :: psd
!
! TRANSITION METALS
!
if (psd.eq.'V'  .or. psd.eq.'Cr' .or. psd .eq.'Mn' .or. psd.eq.'Fe' .or. &
    psd.eq.'Co' .or. psd.eq.'Ni' .or. psd .eq.'Cu'.or. psd .eq.'Fe1'.or. &
    psd .eq.'Fe2' ) then
    hubbard_l = 2
!
! RARE EARTHS
!
elseif (psd .eq.'Ce') then
   hubbard_l =  3
!
! OTHER ELEMENTS
!
elseif (psd .eq.'H') then
   hubbard_l =  0
elseif (psd .eq.'O') then
   hubbard_l = 1
else
   hubbard_l = -1
   call errore ('set_Hubbard_l','pseudopotential not yet inserted', 1)
endif
return
end function set_Hubbard_l
!
!-----------------------------------------------------------------------
      subroutine new_ns(c,eigr,betae,hpsi,hpsi_pen,forceh)
!-----------------------------------------------------------------------
!
! This routine computes the on site occupation numbers of the Hubbard ions.
! It also calculates the contribution of the Hubbard Hamiltonian to the
! electronic potential and to the forces acting on ions.
!
      use control_flags,      ONLY: tfor, tprnfor
      use kinds,              ONLY: DP        
      use ions_base,          only: na, nat, nsp
      use gvecw,              only: ngw
      use gvect, only: gstart
      USE uspp,               ONLY: nhsa=>nkb
      USE uspp_param,         ONLY: upf
      use electrons_base,     only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE ldaU_cp,               ONLY: lda_plus_u, Hubbard_U, Hubbard_l
      USE ldaU_cp,               ONLY: n_atomic_wfc, ns, e_hubbard
      USE step_penalty,       ONLY: E_pen, A_pen, sigma_pen, alpha_pen
      USE step_penalty,       ONLY: step_pen
      USE dspev_module,       only: dspev_drv
      USE mp_global,       only: nbgrp
!
      implicit none
#ifdef __PARA
      include 'mpif.h'
#endif
      integer, parameter :: ldmx = 7
      complex(DP), intent(in) :: c(ngw,nx), eigr(ngw,nat),      &
     &                               betae(ngw,nhsa)
      complex(DP), intent(out) :: hpsi(ngw,nx), hpsi_pen(ngw,nx)
      real(DP) forceh(3,nat), force_pen(3,nat)
!
      complex(DP), allocatable:: wfc(:,:), swfc(:,:),dphi(:,:,:),   &
     &                               spsi(:,:)
      real(DP), allocatable   :: becwfc(:,:), bp(:,:),              &
     &                               dbp(:,:,:), wdb(:,:,:)
      real(DP), allocatable   :: dns(:,:,:,:)
      real(DP), allocatable   :: e(:), z(:,:),                      &
     &                               proj(:,:), temp(:)
      real(DP), allocatable   :: ftemp1(:), ftemp2(:)
      real(DP)                :: lambda(ldmx), somma, ntot, nsum,   &
     &                           nsuma, x_value, g_value, step_value
      real(DP) :: f1 (ldmx*ldmx), vet (ldmx, ldmx)
      integer is, ia, iat, nb, isp, l, m, m1, m2, k, i, counter, err, ig
      integer iv, jv, inl, jnl,alpha,alpha_a,alpha_s,ipol
      integer, allocatable ::  offset (:,:)
      complex(DP) :: tempsi
!
      if( nbgrp > 1 ) &
         call errore(' new_ns ', ' parallelization over bands not yet implemented ', 1 )
!
      allocate(wfc(ngw,n_atomic_wfc))
      allocate(ftemp1(ldmx))
      allocate(ftemp2(ldmx))
!
! calculate wfc = atomic states
!
!!!      call ewfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
      allocate(becwfc(nhsa,n_atomic_wfc))
!!!      call nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)
!
      allocate(swfc(ngw,n_atomic_wfc))
!!!      call s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate proj = <c|S|wfc>
!
      allocate(proj(n,n_atomic_wfc))
      CALL projwfc_hub( c, nx, eigr, betae, n, n_atomic_wfc,            &
     & wfc, becwfc, swfc, proj ) !#@
!
      allocate(offset(nsp,nat))
      counter = 0
      do is = 1, nsp
         do ia = 1, na(is)
            do i = 1, upf(is)%nwfc
               l = upf(is)%lchi(i)
               if (l.eq.Hubbard_l(is)) offset (is,ia) = counter
               counter = counter + 2 * l + 1
            end do
         end do
      end do
      if (counter.ne.n_atomic_wfc)                                      &
     &                 call errore ('new_ns','nstart<>counter',1)
      ns(:,:,:,:) = 0.d0
      iat = 0
      do is = 1,nsp
         do ia = 1,na(is)
            iat = iat + 1
            if (Hubbard_U(is).ne.0.d0) then 
               k = offset(is,ia)
               do m1 = 1, 2*Hubbard_l(is) + 1
                  do m2 = m1, 2*Hubbard_l(is) + 1
                     do i = 1,n
!                      write(6,*) i,ispin(i),f(i)
                      ns(iat,ispin(i),m1,m2) = ns(iat,ispin(i),m1,m2) + &
     &                               f(i) * proj(i,k+m2) * proj(i,k+m1)
                     end do
!                     ns(iat,:,m2,m1) = ns(iat,:,m1,m2)
                     ns(iat,1,m2,m1) = ns(iat,1,m1,m2)
                     ns(iat,2,m2,m1) = ns(iat,2,m1,m2)
                  end do
               end do
            end if
         end do
      end do
      if (nspin.eq.1) ns = 0.5d0 * ns
! Contributions to total energy
      e_hubbard = 0.d0
      iat = 0
      do is = 1,nsp
         do ia = 1,na(is)
            iat=iat + 1
            if (Hubbard_U(is).ne.0.d0) then
                k = offset(is,ia)
                do isp = 1,nspin
                   do m1 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard + 0.5d0 * Hubbard_U(is) *    &
     &                           ns(iat,isp,m1,m1)
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        e_hubbard = e_hubbard - 0.5d0 * Hubbard_U(is) * &
     &                              ns(iat,isp,m1,m2) * ns(iat,isp,m2,m1)
                     end do
                   end do
                end do
             end if
         end do
       end do
       if (nspin.eq.1) e_hubbard = 2.d0*e_hubbard
!       if (nspin.eq.1) e_lambda = 2.d0*e_lambda
!
!      Calculate the potential and forces on wavefunctions due to U
!
      hpsi(:,:)=(0.d0,0.d0)
      iat=0
      do is = 1, nsp
         do ia=1, na(is)
            iat = iat + 1
            if (Hubbard_U(is).ne.0.d0) then
               do i=1, n
                  do m1 = 1, 2 * Hubbard_l(is) + 1
                     tempsi = proj (i,offset(is,ia)+m1)
                     do m2 = 1, 2 * Hubbard_l(is) + 1
                        tempsi = tempsi - 2.d0 * ns(iat,ispin(i),m1,m2)*&
     &                                proj (i,offset(is,ia)+m2)
                     enddo
                     tempsi = tempsi * Hubbard_U(is)/2.d0*f(i)
                     call zaxpy (ngw,tempsi,swfc(1,offset(is,ia)+m1),1, &
     &                           hpsi(1,i),1)
                  enddo
               enddo
            endif
         enddo
      enddo
!
!      Calculate the potential and energy due to constraint
!
      hpsi_pen(:,:)=0.d0

      if (step_pen) then
         iat=0
         E_pen=0
         do is = 1,nsp
            do ia = 1, na(is)
               nsuma = 0.d0
               iat = iat + 1
              if (Hubbard_U(is).ne.0.d0) then
               do isp = 1, nspin
                  if (A_pen(iat,isp).ne.0.0) then
                     k = 0
                     f1=0.0
                     do m1 = 1, 2 * Hubbard_l(is) + 1
                        do m2 = m1, 2 * Hubbard_l(is) + 1
                           k = k + 1
                           f1 (k) = ns (iat, isp, m2, m1)
                        enddo
                     enddo
                     CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1, f1, lambda, vet, ldmx  )
                     x_value=alpha_pen(iat)-lambda(2*Hubbard_l(is)+1)
                     call stepfn(A_pen(iat,isp),sigma_pen(iat),x_value, &
     &                           g_value,step_value)
                     do i=1, n
                        do m1 = 1, 2 * Hubbard_l(is) + 1
                           do m2 = 1, 2 * Hubbard_l(is) + 1
                              tempsi=-1.d0*f(i)*proj (i,offset(is,ia)+m1)* &
     &                       vet(m1,2*Hubbard_l(is)+1)*vet(m2,2*Hubbard_l(is)+1)*g_value
                              call ZAXPY (ngw,tempsi,swfc(1,offset(is,ia)+ &
     &                           m2),1,hpsi_pen(1,i),1)
                           enddo
                        enddo
                     end do
                     E_pen=E_pen+step_value
                  end if
               enddo
              endif
            enddo
         enddo
      endif

!
! Calculate the contribution to forces on ions due to U and constraint
!
      forceh=0.d0

      force_pen=0.d0

      if ((tfor).or.(tprnfor)) then
        allocate (bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3))
        allocate(dns(nat,nspin,ldmx,ldmx))
        allocate (spsi(ngw,n))
!
        call nlsm1 (n,1,nsp,eigr,c,bp)
        call s_wfc(n,bp,betae,c,spsi)
        call nlsm2_bgrp( ngw, nhsa, eigr, c, dbp, nx, n )
        call nlsm2_bgrp( ngw, nhsa, eigr, wfc, wdb, n_atomic_wfc, n_atomic_wfc )
!
        alpha=0
        do alpha_s = 1, nsp
         do alpha_a = 1, na(alpha_s)
            alpha=alpha+1
            do ipol = 1,3
               call dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,      &
     &                    offset,c,wfc,eigr,betae,proj,ipol,dns)
               iat=0
               do is = 1, nsp
                  do ia=1, na(is)
                     iat = iat + 1
                     if (Hubbard_U(is).ne.0.d0) then
                        do isp = 1,nspin
                           do m2 = 1,2*Hubbard_l(is) + 1
                              forceh(ipol,alpha) = forceh(ipol,alpha) -            &
     &                        Hubbard_U(is) * 0.5d0 * dns(iat,isp,m2,m2)
                              do m1 = 1,2*Hubbard_l(is) + 1
                                 forceh(ipol,alpha) = forceh(ipol,alpha) +         &
     &                           Hubbard_U(is)*ns(iat,isp,m2,m1)*       &
     &                           dns(iat,isp,m1,m2)
                              end do
                           end do
                        end do
                     end if
! Occupation constraint add here
                     if (step_pen) then
                        do isp = 1, nspin
                           if ((A_pen(iat,isp).ne.0.0).and.           &
      &                       (Hubbard_U(is).ne.0.d0)) then
                              k = 0
                              f1=0.0
                              do m1 = 1, 2 * Hubbard_l(is) + 1
                                 do m2 = m1, 2 * Hubbard_l(is) + 1
                                    k = k + 1
                                    f1 (k) = ns (iat, isp, m2, m1)
                                 enddo
                              enddo
                              CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1, f1, lambda, vet, ldmx  )
                              x_value=alpha_pen(iat)-lambda(2*Hubbard_l(is)+1)
                              call stepfn(A_pen(iat,isp),sigma_pen(iat),x_value,g_value,&
     &                             step_value)
                              do m1 = 1,2*Hubbard_l(is) + 1
                                 do m2 = 1,2*Hubbard_l(is) + 1
                                    force_pen(ipol,alpha) =                &
     &                              force_pen(ipol,alpha) +                &
     &                              g_value * dns(iat,isp,m1,m2)           &
     &                              *vet(m1,2*Hubbard_l(is)+1)*vet(m2,2*Hubbard_l(is)+1)
                                 end do
                              end do
                           endif
                        end do
                     end if
                  end do
               end do
            end do
         end do
        end do
        if (nspin.eq.1) then
           forceh = 2.d0 * forceh
           force_pen=2.d0 * force_pen
        end if

        forceh = forceh + force_pen
!
        deallocate ( wfc, becwfc, spsi, proj, offset, swfc, dns, bp, dbp, wdb)
      end if
      return
      end subroutine new_ns
!
!
!
!-----------------------------------------------------------------------
      subroutine write_ns
!-----------------------------------------------------------------------
!
! This routine computes the occupation numbers on atomic orbitals.
! It also write the occupation number in the output file.
!
      USE kinds,            only: DP
      USE constants,        ONLY: autoev
      use electrons_base,   only: nspin
      use electrons_base,   only: n => nbsp 
      use ions_base,        only: na, nat, nsp
      use gvecw,            only: ngw
      USE ldaU_cp,             ONLY: lda_plus_u, Hubbard_U, Hubbard_l
      USE ldaU_cp,             ONLY: n_atomic_wfc, ns, e_hubbard
      USE ldaU_cp,             ONLY: Hubbard_lmax
      use dspev_module,     only : dspev_drv
      USE step_penalty,     ONLY: step_pen, A_pen, sigma_pen, alpha_pen

      implicit none

  integer :: is, isp, ia, m1, m2, ldim, iat, err, k
! cpunter on atoms type
! counter on spin component
! counter on atoms
! counter on wavefn
! counters on d components
  integer, parameter :: ldmx = 7
  real(DP), allocatable   :: ftemp1(:), ftemp2(:)
  real(DP) :: f1 (ldmx * ldmx), vet (ldmx, ldmx)
  real(DP) :: lambda (ldmx), nsum, nsuma
  write (*,*) 'enter write_ns'

  if ( 2 * Hubbard_lmax + 1 .gt. ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

  if (step_pen) then
     do isp=1,nspin
        write (6,'(6(a,i2,a,i2,a,f8.4,6x))') &
        ('A_pen(',is,',',isp,') =', A_pen(is,isp),is=1,nsp)
     enddo
     write (6,'(6(a,i2,a,f8.4,6x))') &
           ('sigma_pen(',is,') =', sigma_pen(is), is=1,nsp)
     write (6,'(6(a,i2,a,f8.4,6x))') &
        ('alpha_pen(',is,') =', alpha_pen(is), is=1,nsp)
  endif

  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('U(',is,') =', Hubbard_U(is) * autoev, is=1,nsp)
!  write (6,'(6(a,i2,a,f8.4,6x))') &
!        ('alpha(',is,') =', Hubbard_alpha(is) * autoev, is=1,nsp)
      nsum = 0.d0
      allocate(ftemp1(ldmx))
      allocate(ftemp2(ldmx))
      iat = 0
      write(6,*) 'nsp',nsp
      do is = 1,nsp
         do ia = 1, na(is)
            nsuma = 0.d0
            iat = iat + 1
!        if (iat.eq.1) then
            if (Hubbard_U(is).ne.0.d0) then
               do isp = 1, nspin
                   do m1 = 1, 2 * Hubbard_l(is) + 1
                      nsuma = nsuma + ns (iat, isp, m1, m1)
                   end do
               end do
               if (nspin.eq.1) nsuma = 2.d0 * nsuma
               write(6,'(a,1x,i2,2x,a,f11.7)') 'atom', iat,              &
     &                                      ' Tr[ns(na)]= ',nsuma
               nsum = nsum + nsuma
!
               do isp = 1, nspin

                  k = 0
                  do m1 = 1, 2 * Hubbard_l(is) + 1
                     do m2 = m1, 2 * Hubbard_l(is) + 1
                        k = k + 1
                        f1 ( k ) = ns (iat, isp, m2, m1)
                     enddo
                  enddo

                  CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1, f1, lambda, vet, ldmx  )

                  write(6,'(a,1x,i2,2x,a,1x,i2)') 'atom', iat, 'spin', isp
                  write(6,'(a,7f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,&
     &                                2 * Hubbard_l(is) + 1)
                  write(6,*) 'eigenvectors'
                  do m2 = 1, 2*Hubbard_l(is)+1
                     write(6,'(i2,2x,7(f10.7,1x))') m2,(real(vet(m1,m2)),&
     &                            m1=1,2 * Hubbard_l(is) + 1)
                  end do
                  write(6,*) 'occupations'
                  do m1 = 1, 2*Hubbard_l(is)+1
                     write (6,'(7(f6.3,1x))') (ns(iat,isp,m1,m2),m2=1,    &
     &                     2*Hubbard_l(is)+1)
                  end do
               end do
            end if
!        end if
         end do
      end do
      deallocate ( ftemp1, ftemp2)
      return
      end subroutine write_ns
!-----------------------------------------------------------------------
      subroutine genatwfc(n_atomic_wfc,atwfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space, in the same order as used in new_ns
!
      use ions_base,          only: na, nsp
      use gvecw,              only: ngw
      use gvect, only: gg, g, gstart
      use cell_base,          only: omega, tpiba
      use constants,          only: fpi
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
      USE kinds,              ONLY: DP
!
      implicit none
      integer, intent(in) :: n_atomic_wfc
      complex(DP), intent(out):: atwfc(ngw,n_atomic_wfc)
!
      integer natwfc, is, ia, ir, nb, l, m, lm, i, lmax_wfc, ig
      real(DP), allocatable::  ylm(:,:), q(:), jl(:), vchi(:),        &
     &     chiq(:), gxn(:,:)
!
      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' genatwfc ', ' rgrid not allocated ', 1 )
!
      allocate(q(ngw))
      allocate(gxn(3,ngw))
      allocate(chiq(ngw))
!
      do ig=1,ngw
         q(ig) = sqrt(gg(ig))*tpiba
      end do
      if (gstart == 2) gxn(1,:)=0.0d0
      do ig=gstart,ngw
         gxn(:,ig) = g(:,ig)/sqrt(gg(ig)) !ik<=>ig
      end do
!
      natwfc=0
!#@@@@
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX (lmax_wfc, MAXVAL ( upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, g, gg, ylm)
!#@@@@

      do is = 1, nsp
         ALLOCATE  ( jl(rgrid(is)%mesh), vchi(rgrid(is)%mesh) )
         do ia=1,na(is)
!
!   radial fourier transform of the chi functions
!   NOTA BENE: chi is r times the radial part of the atomic wavefunction
!              bess requires l+1, not l, on input
!
            do nb = 1,upf(is)%nwfc
               l = upf(is)%lchi(nb)
               do i=1,ngw
                  call sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
                  do ir=1,rgrid(is)%mesh
                     vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
                  enddo
                  call simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
               enddo
!
!   multiply by angular part and structure factor
!   NOTA BENE: the factor i^l MUST be present!!!
!
               do m = 1,2*l+1
                  lm = l**2 + m
!                  call ylmr2b(lm,ngw,ngw,gxn,ylm)
                  natwfc = natwfc + 1
                  atwfc(:,natwfc) = (0.d0,1.d0)**l * ylm(:,lm)*chiq(:)
               enddo
            enddo
         end do
         DEALLOCATE  ( vchi, jl )
      end do
!
      do i = 1,natwfc
        call dscal(2*ngw,fpi/sqrt(omega),atwfc(1,i),1)
      end do
!
      if (natwfc.ne.n_atomic_wfc)                                       &
     &     call errore('atomic_wfc','unexpected error',natwfc)
!
      deallocate(ylm)
      deallocate(chiq)
      deallocate(gxn)
      deallocate(q)
!
      return
      end subroutine genatwfc
!
!-------------------------------------------------------------------------
      subroutine dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,         &
     &                  offset,c,wfc,                                   &
     &                  eigr,betae,                                     &
     &                  proj,ipol,dns)
!-----------------------------------------------------------------------
!
! This routine computes the derivative of the ns with respect to the ionic
! displacement tau(alpha,ipol) used to obtain the Hubbard contribution to the
! atomic forces.
!
      use ions_base, only: na, nat, nsp
      use gvecw, only: ngw
      use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE uspp,           ONLY: nhsa=>nkb
      USE ldaU_cp,           ONLY: Hubbard_U, Hubbard_l
      USE ldaU_cp,           ONLY: n_atomic_wfc, ns
      USE kinds,          ONLY: DP
!
      implicit none
      integer, parameter :: ldmx = 7
      integer ibnd,is,i,ia,counter, m1,m2, l, iat, alpha, ldim
! input
      integer,      intent(in) :: offset(nsp,nat)
      integer,      intent(in) :: alpha_a,alpha_s,ipol
      real(DP),     intent(in) :: wfc(ngw,n_atomic_wfc),  c(2,ngw,nx),  &
     &                            eigr(2,ngw,nat),betae(2,ngw,nhsa),    &
     &                            becwfc(nhsa,n_atomic_wfc),            &
     &                            bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3)
      real(DP),     intent(in) :: proj(n,n_atomic_wfc)
      complex (DP), intent(in) :: spsi(ngw,n)
! output
      real (DP),   intent(out) :: dns(nat,nspin,ldmx,ldmx)
!
!     dns !derivative of ns(:,:,:,:) w.r.t. tau
!
      real (DP),   allocatable :: dproj(:,:)
!
!     dproj(n,n_atomic_wfc) ! derivative of proj(:,:) w.r.t. tau 
!
      allocate (dproj(n,n_atomic_wfc) )
!
      dns(:,:,:,:) = 0.d0
!
          call dprojdtau(c,wfc,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,     &
     &                   alpha_s,ipol,offset(alpha_s,alpha_a),dproj)
!
! compute the derivative of occupation numbers (the quantities dn(m1,m2))
! of the atomic orbitals. They are real quantities as well as n(m1,m2)
!
      iat=0
      do is=1,nsp
         do ia = 1,na(is)
            iat=iat+1
            if (Hubbard_U(is).ne.0.d0) then
               ldim = 2*Hubbard_l(is) + 1
               do m1 = 1, ldim
                  do m2 = m1, ldim
                     do ibnd = 1,n
                        dns(iat,ispin(ibnd),m1,m2) =                    &
     &                  dns(iat,ispin(ibnd),m1,m2) +                    &
     &                   f(ibnd)*REAL(  proj(ibnd,offset(is,ia)+m1) *   &
     &                   (dproj(ibnd,offset(is,ia)+m2))  +              &
     &                         dproj(ibnd,offset(is,ia)+m1)  *          &
     &                         (proj(ibnd,offset(is,ia)+m2)) )
                     end do
                     dns(iat,:,m2,m1) = dns(iat,:,m1,m2)
                  end do
               end do
            end if
         end do
      end do
!
      deallocate (dproj)
      return
      end subroutine dndtau
!
!
!-----------------------------------------------------------------------
      subroutine dprojdtau(c,wfc,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,    &
     &                     alpha_s,ipol,offset,dproj)
!-----------------------------------------------------------------------
!
! This routine computes the first derivative of the projection
! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
!
      use ions_base, only: na, nat
      use gvecw, only: ngw
      use gvect, only: g, gstart
      use electrons_base, only: n => nbsp, nx => nbspx
!      use gvec
!      use constants
      USE uspp,           ONLY: nhsa=>nkb, qq
      USE ldaU_cp,           ONLY: Hubbard_U, Hubbard_l
      USE ldaU_cp,           ONLY: n_atomic_wfc
      use cell_base,      ONLY: tpiba
      USE uspp_param,     only: nh, ish
      use mp_global,      only: intra_bgrp_comm
      use mp,             only: mp_sum
      USE kinds,          ONLY: DP
!
       implicit none
       integer, parameter :: ldmx = 7
       integer alpha_a, alpha_s,ipol, offset
! input: the displaced atom
! input: the component of displacement
! input: the offset of the wfcs of the atom "alpha_a,alpha_s"
       complex (DP), intent(in) :: spsi(ngw,n),                     &
     &                  c(ngw,nx), eigr(ngw,nat)
! input: the atomic wfc
! input: S|evc>
       real(DP), intent(in) ::becwfc(nhsa,n_atomic_wfc),            &
     &                            wfc(2,ngw,n_atomic_wfc),              &
     &            bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3)
       real(DP), intent(out) :: dproj(n,n_atomic_wfc)
! output: the derivative of the projection
!
      integer i,ig,m1,ibnd,iwf,ia,is,iv,jv,ldim,alpha,l,m,k,inl
!
      real(kind=8), allocatable :: gk(:)
!
      complex (DP), allocatable :: dwfc(:,:)
      real (DP), allocatable :: betapsi(:,:),                       &
     &                              dbetapsi(:,:),                      &
     &                              wfcbeta(:,:),wfcdbeta(:,:),temp(:)
!      dwfc(ngw,ldmx),             ! the derivative of the atomic d wfc
!      betapsi(nh,n),              ! <beta|evc>
!      dbetapsi(nh,n),             ! <dbeta|evc>
!      wfcbeta(n_atomic_wfc,nh),   ! <wfc|beta>
!      wfcdbeta(n_atomic_wfc,nh),  ! <wfc|dbeta>
      ldim = 2 * Hubbard_l(alpha_s) + 1
      allocate ( dwfc(ngw,ldmx),betapsi(nh(alpha_s),n))
      allocate ( dbetapsi(nh(alpha_s),n),                               &
     &           wfcbeta(n_atomic_wfc,nh(alpha_s)))
      allocate (wfcdbeta(n_atomic_wfc,nh(alpha_s)) )
      dproj(:,:)=0.d0
!
! At first the derivative of the atomic wfc is computed
!
!
      allocate(gk(ngw))
      allocate(temp(ngw))
!
      if (Hubbard_U(alpha_s).ne.0.d0) then
!
         do ig=1,ngw
            gk(ig)=g(ipol,ig)*tpiba 
!
            do m1=1,ldim
                  dwfc(ig,m1) = CMPLX (gk(ig)*wfc(2,ig,offset+m1),      &
     &                  -1*gk(ig)*wfc(1,ig,offset+m1), kind=dp )
            end do
         end do
!
         do ibnd=1,n
            do m1=1,ldim
               temp(:)=real(conjg(dwfc(:,m1))*spsi(:,ibnd))
               dproj(ibnd,offset+m1)=2.d0*SUM(temp) 
               if (gstart==2) dproj(ibnd,offset+m1)=dproj(ibnd,offset+m1)-temp(1)
            end do
         end do
         call mp_sum( dproj, intra_bgrp_comm )
      end if
      do iv=1,nh(alpha_s)
         inl=ish(alpha_s)+(iv-1)*na(alpha_s)+alpha_a
         do i=1,n
            betapsi(iv,i)=bp(inl,i)
            dbetapsi(iv,i)=dbp(inl,i,ipol)
         end do
         do m=1,n_atomic_wfc
!                 do m1=1,2**Hubbard_l(is) + 1
            wfcbeta(m,iv)=becwfc(inl,m)
            wfcdbeta(m,iv)=wdb(inl,m,ipol)
         end do
      end do
      do ibnd=1,n
         do iv=1,nh(alpha_s)
            do jv=1,nh(alpha_s)
               do m=1,n_atomic_wfc
!                       do m1=1,2**Hubbard_l(is) + 1
                  dproj(ibnd,m) =                                       &
     &                        dproj(ibnd,m) + qq(iv,jv,alpha_s) *       &
     &                         ( wfcdbeta(m,iv)*betapsi(jv,ibnd) +      &
     &                           wfcbeta(m,iv)*dbetapsi(jv,ibnd) )
               end do
            end do
         end do
      end do
      deallocate(temp, gk)
      deallocate (betapsi)
      deallocate (dwfc)
      deallocate (dbetapsi)
      deallocate (wfcbeta)
      deallocate (wfcdbeta)
      return
      end subroutine dprojdtau
!
!
!-----------------------------------------------------------------------
      subroutine stepfn(A,sigma,x_value,g_value,step_value)
!-----------------------------------------------------------------------
!     This subroutine calculates the value of the gaussian and step
!     functions with a given x_value. A and sigma are given in the
!     input file. ... to be used in occupation_constraint...
!
      USE constants, ONLY : pi
      implicit none
      real(kind=8) A, sigma, x_value, g_value, step_value
      real(kind=8) x
      integer i
      step_value=0.0d0
      g_value=0.0d0
!
      do i=1,100000
         x=x_value + (i-100000)/100000.0d0*(x_value + 5.d0*sigma)
!
! Integrate from 5 sigma before the x_value
!
         g_value=A*dexp(-x*x/(2*sigma*sigma))/(sigma*dsqrt(2*pi))
!         write(6,*) 'step', step_value,'g',g_value
!         if (g_value.le.0.0) g_value=0.0
         if ((x_value+5*sigma).ge.0.0d0) then
         step_value=step_value+g_value/100000.0d0*(x_value+5.d0*sigma)
         end if
      end do
      return
      end subroutine stepfn
!
!-----------------------------------------------------------------------
      SUBROUTINE projwfc_hub( c, nx, eigr, betae, n, n_atomic_wfc,  &
     & wfc, becwfc, swfc, proj )
!-----------------------------------------------------------------------
      !
      ! Projection on atomic wavefunctions
      ! Atomic wavefunctions are not orthogonized
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: autoev
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_bgrp_comm
      USE mp,                 ONLY: mp_sum
      USE gvecw,              ONLY: ngw
      USE gvect, ONLY: gstart
      USE ions_base,          ONLY: nsp, na, nat
      USE uspp,               ONLY: nhsa => nkb
!
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: nx, n, n_atomic_wfc
      COMPLEX(DP), INTENT(IN) :: c( ngw, nx ), eigr(ngw,nat), betae(ngw,nhsa)
!
      COMPLEX(DP), INTENT(OUT):: wfc(ngw,n_atomic_wfc),    &
     & swfc( ngw, n_atomic_wfc )
      real(DP), intent(out):: becwfc(nhsa,n_atomic_wfc) !DEBUG
      REAL(DP),    ALLOCATABLE :: overlap(:,:), e(:), z(:,:)
      REAL(DP),    ALLOCATABLE :: temp(:)
      REAL(DP)                 :: somma, proj(n,n_atomic_wfc)
      INTEGER :: is, ia, nb, l, m, k, i
      !
      ! calculate number of atomic states
      !
      !
      IF ( n_atomic_wfc .EQ. 0 ) RETURN
      !
      !
      ! calculate wfc = atomic states
      !
      CALL atomic_wfc_northo( eigr, n_atomic_wfc, wfc )
      !
      ! calculate bec = <beta|wfc>
      !
      CALL nlsm1( n_atomic_wfc, 1, nsp, eigr, wfc, becwfc )
      !
      ! calculate swfc = S|wfc>
      !
      CALL s_wfc( n_atomic_wfc, becwfc, betae, wfc, swfc )
      !
      ! calculate proj = <c|S|wfc>
      !
      ALLOCATE(temp(ngw))
      DO m=1,n
         DO l=1,n_atomic_wfc
            temp(:)=DBLE(CONJG(c(:,m))*swfc(:,l)) !#@@@
            proj(m,l)=2.d0*SUM(temp)
            IF (gstart == 2) proj(m,l)=proj(m,l)-temp(1)
         END DO
      END DO
      DEALLOCATE(temp)
      CALL mp_sum( proj, intra_bgrp_comm )
!
      RETURN
      END SUBROUTINE projwfc_hub
!
!-----------------------------------------------------------------------
      SUBROUTINE atomic_wfc_northo( eigr, n_atomic_wfc, wfc )
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
! Atomic wavefunctions not orthogonalized
!
      USE kinds,              ONLY: DP
      USE gvecw,              ONLY: ngw
      USE gvect, ONLY: gstart, gg, g
      USE ions_base,          ONLY: nsp, na, nat
      USE cell_base,          ONLY: tpiba, omega !#@@@
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
!#@@@@
      USE constants,          ONLY: fpi
!#@@@@
!
      IMPLICIT NONE
      INTEGER,     INTENT(in) :: n_atomic_wfc
      COMPLEX(DP), INTENT(in) :: eigr( ngw, nat )
      COMPLEX(DP), INTENT(out):: wfc( ngw, n_atomic_wfc )
!
      INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
      REAL(DP), ALLOCATABLE ::  ylm(:,:), q(:), jl(:), vchi(:), chiq(:)

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' atomic_wfc_northo ', ' rgrid not allocated ', 1 )
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, g, gg, ylm)
      ndm = MAXVAL(rgrid(1:nsp)%mesh)
      !
      ALLOCATE(jl(ndm), vchi(ndm))
      ALLOCATE(q(ngw), chiq(ngw))
!
      DO i=1,ngw
         q(i) = SQRT(gg(i))*tpiba
      END DO
!
      natwfc=0
      isa   = 0
      DO is=1,nsp
         !
         !   radial fourier transform of the chi functions
         !   NOTA BENE: chi is r times the radial part of the atomic wavefunction
         !
         DO ia = 1 + isa, na(is) + isa
            DO nb = 1,upf(is)%nwfc
               l = upf(is)%lchi(nb)
               DO i=1,ngw
                  CALL sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
                  DO ir=1,rgrid(is)%mesh
                     vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
                  ENDDO
                  CALL simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
               ENDDO
               !
               !   multiply by angular part and structure factor
               !   NOTA BENE: the factor i^l MUST be present!!!
               !
               DO m = 1,2*l+1
                  lm = l**2 + m
                  !DO ia = 1 + isa, na(is) + isa
                  natwfc = natwfc + 1
                  wfc(:,natwfc) = (0.d0,1.d0)**l * eigr(:,ia)* ylm(:,lm)*chiq(:)
                  !ENDDO
               ENDDO
            ENDDO
         ENDDO
         isa = isa + na(is)
      ENDDO
!
      IF (natwfc.NE.n_atomic_wfc)                                       &
     &     CALL errore('atomic_wfc','unexpected error',natwfc)
!
!#@@@@
      do i = 1,n_atomic_wfc
        call dscal(2*ngw,fpi/sqrt(omega),wfc(1,i),1)
      end do
!#@@@@
      DEALLOCATE(q, chiq, vchi, jl, ylm)
!
      RETURN
      END SUBROUTINE atomic_wfc_northo
