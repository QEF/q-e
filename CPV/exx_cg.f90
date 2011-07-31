      
!-------------------------------------------------------------------------------
!HPOTCG -- using Conjugate Gradient method to compute the Hartree Potential.
!-------------------------------------------------------------------------------

      subroutine hpotcg(rho,pot)

      USE kinds,                   ONLY  :  DP
      USE cp_main_variables,       ONLY  :  odtothd_in_sp, thdtood_in_sp
      USE cp_main_variables,       ONLY  :  thdtood
      USE cp_main_variables,       ONLY  :  n=>np_in_sp, n2=>np_in_sp2
      USE cp_main_variables,       ONLY  :  coeke, nord2
      USE wannier_base,            ONLY  : poisson_eps


      implicit real(DP) (a-h,o-z)
 
      real(DP)  rho(n), pot(n)
      real(DP), ALLOCATABLE  :: wk(:), wk2(:)
!
      integer  i, iou, ipar(16), lwk, mvstep
      real*8   fpar(16), dnrm2
      external dnrm2
!
      lwk = n*5
      ALLOCATE( wk(lwk) )
      ALLOCATE( wk2(n+n2) )

!     set up the parameter arrays
!
      ipar(1) = 0
      ipar(2) = 0
      ipar(3) = 1
      ipar(4) = lwk
      ipar(5) = 5
      ipar(6) = 300
      fpar(1) = 1.0D-6
      fpar(2) = poisson_eps
      fpar(11) = 0.0D0
!
      iou = 6
      mvstep = 0
 10   call cg(n,rho,pot,ipar,fpar,wk)
!
      if (ipar(1).eq.1) then
!        write (iou, *) ipar(7), fpar(5), mvstep
         call dcopy(n, wk(ipar(8)), 1, wk2, 1)
         call start_clock('lapmv')
         call lapmvs( wk2, wk(ipar(9)) )
         call stop_clock('lapmv')
         mvstep = mvstep + 1
         fpar(11) = fpar(11) + 74*n
         goto 10
      else if (ipar(1).le.0) then
         if (ipar(1).eq.0) then
            print *, 'Iterative sovler has satisfied convergence test.'
         else if (ipar(1).eq.-1) then
            print *, 'Iterative potver has iterated too many times.'
         else if (ipar(1).eq.-2) then
            print *, 'Iterative potver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4),  &
     &           ' elements.'
         else if (ipar(1).eq.-3) then
            print *, 'Iterative sovler is facing a break-down.'
            print *, 'ipar(12) =', ipar(12)
         else
            print *, 'Iterative potver terminated. code =', ipar(1)
         endif
      endif
!     write (iou, *) ipar(7), DBLE(fpar(6))
!     write (iou, *) '# ',
!    +     ipar(7), ' MATVECs   ', DBLE(fpar(11)), ' OPS'
      write (iou, *) '# retrun code = ', ipar(1),  '   cgstep = ', mvstep
!     write (iou, *) 'fpar follows'
      write (iou, *) (fpar(i),i=1,7)
!
!     check the error
      call lapmvs(pot, wk)
      do i = 1, n
         wk(i) = wk(i) - rho(i)
      enddo
      write (iou, *) '# the residual norm ', DBLE(dnrm2(n,wk,1)), DBLE(fpar(5))

      DEALLOCATE( wk, wk2 )

      return
      end
!-----end-of-hpotcg
!-----------------------------------------------------------------------
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
      end
!-----end-of-distdot
