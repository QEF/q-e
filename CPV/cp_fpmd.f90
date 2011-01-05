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

        SUBROUTINE newgb( a1, a2, a3, omega, alat )
!
!     re-generation of little box g-vectors
!
          USE kinds, ONLY: DP
          USE grid_dimensions, only: nr1, nr2, nr3
          USE smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
          USE small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
          USE constants, ONLY: pi

          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: a1( 3 ), a2( 3 ), a3( 3 ), omega, alat

          INTEGER :: i
          REAL(DP) :: alatb, b1b(3),b2b(3),b3b(3)

          IF ( nr1b == 0 .OR. nr2b == 0 .OR. nr3b == 0 ) return
          alatb  = alat / nr1*nr1b
          tpibab = 2.d0*pi / alatb
          do i=1,3
            a1b(i)=a1(i)/nr1*nr1b
            a2b(i)=a2(i)/nr2*nr2b
            a3b(i)=a3(i)/nr3*nr3b
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
