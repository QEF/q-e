!
!-----------------------------------------------------------------------
      subroutine write_pseudo &
           (iunps,zed,xmin,dx,mesh,ndm,r,r2, &
           dft,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
           rhoc,vnl,phis,vpsloc,lls,ocs,etots,nwfs)
!-----------------------------------------------------------------------
!
      implicit none
      integer, parameter :: dp = kind(1.d0)
      integer :: ndm, mesh, nwfps, lmin,lmax,lloc,nlc,nnl,nwfs, &
              lls(nwfs)
      real(kind=dp) :: zed, zval, xmin,dx, cc(2),alpc(2),alc(6,0:3), &
             alps(3,0:3), phis(ndm,nwfs), ocs(nwfs), &
             r(ndm), r2(ndm), vnl(ndm,0:3), rhoc(ndm), erf, etots
      integer ios, mdum, i, l, k, n, ir, iunps, nb, ldum
      real(kind=dp) :: zdum, fourpi, vnloc,a_core, b_core,  &
             alfa_core, dum, xdum, dxdum, rdum, vpsloc(ndm)
      logical nlcc, bhstype, numeric
      parameter(fourpi=4.d0*3.141592653589793d0)
      character(len=3)  title_pseudo*70, cdum
      character(len=*) dft
!
!
      cdum='cc'
      nlc=0
      nnl=0
      bhstype=.false.
      if (dft.eq.'PW') then
         write( iunps, '(a)', err=300, iostat=ios ) 'slater-pz-ggx-ggc'
      else
         write( iunps, '(a)', err=300, iostat=ios ) dft
      endif

      write ( iunps, '('''''''',a2,'''''''',f8.4,3i5,l4,i5,l4,e17.9)', &
                         err=300, iostat=ios ) cdum, &
           zval, lmax, nlc, nnl, nlcc, &
           lloc, bhstype, etots

!
!   In numeric pseudopotentials both nlc and nnl are zero.
!
      numeric = nlc.le.0 .and. nnl.le.0

      if (.not.numeric) then
         write( iunps, *, err=300, iostat=ios ) &
             ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
         do l = 0, lmax
            write ( iunps, *, err=300, iostat=ios ) &
                ( alps(i,l),i=1,3 ), (alc(i,l),i=1,6)
         enddo
         if (nlcc) then
            write( iunps, *, err=300, iostat=ios ) a_core, &
                b_core, alfa_core
         endif
      endif

      write( iunps, '(3f15.10,2i5)', err=300, iostat=ios ) &
           zed, xmin, dx, mesh, nwfs

      if (numeric) then
!
!      pseudopotenziali in forma numerica
!
         do l = 0, lmax
               write( iunps, '(a)', err=300, iostat=ios )
               write( iunps, '(4e19.11)', err=300, iostat=ios ) &
                 (vnl(ir,l),ir=1,mesh)
         enddo
         if (lloc.eq.-1) then
            write( iunps, '(a)', err=300, iostat=ios )
            write( iunps, '(4e19.11)', err=300, iostat=ios ) &
                 (vpsloc(ir),ir=1,mesh)
         endif
         if(nlcc) then
            write( iunps, '(4e19.11)', err=300, iostat=ios ) &
                  ( rhoc(ir)/r2(ir)/fourpi, ir=1,mesh )
         endif
      endif

      do l=0,lmax
         do nb = 1, nwfs
            if (lls(nb).eq.l) then
               write( iunps, '(a)', err=300, iostat=ios )
               write( iunps, *, err=300, iostat=ios ) lls(nb), &
                                  (abs(ocs(nb))+ ocs(nb))*0.5d0 
               write( iunps, '(4e19.11)', err=300, iostat=ios ) &
                            (phis(ir,nb),ir=1,mesh)
            endif
         enddo
      enddo
300   call errore('ld1_readin','reading pseudo file',abs(ios))

      return
      end
