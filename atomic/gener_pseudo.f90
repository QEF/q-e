!
!--------------------------------------------------------------------------
      subroutine gener_pseudo
!--------------------------------------------------------------------------
!
!     This routine generate a pseudopotential in completely separated
!     form. It can be of NC type or of US type
!     Multiple projections are allowed.
!     Spin-orbit splitted pseudopotentials are also available.
!     NB: bmat indeces are as in the Vanderbilt paper PRB (1990)
!
!     The output of the routine are:
!
!     phis: the pseudo wavefunctions
!     betas: the nonlocal projectors
!     bmat:  the pseudopotential coefficients
!     qq:    the integrals of the q functions
!     qvan:  the augmentation functions
!     vpsloc: the local pseudopotential
!     chis:   auxiliary functions
!
!
!      
use ld1inc
implicit none

integer ::   &
      ik,    &  ! the point corresponding to rc
      ikus,  &  ! the point corresponding to rc ultrasoft
      ikloc, &  ! the point corresponding to rc locale
      ns,    &  ! counter on pseudo functions
      ns1,   &  ! counter on pseudo functions
      ib,jb, &  ! counter on beta functions
      nnode, &  ! the number of nodes of phi
      lam       ! the angualar momentum

real(kind=dp) ::    &
      xc(8),        & ! parameters of bessel functions
      gi(ndm,2),    & ! auxiliary to compute the integrals
      sum,          & ! the integral
      int_0_inf_dr, & ! the integral function
      work(nwfsx)     ! work space
     
real(kind=dp) :: &
      binv(nwfsx,nwfsx)  ! the inverse B matrix

integer :: &
      m, n, l, n1, n2, nwf0, nst, ikl, imax, iwork(nwfsx), &
      is, nbf, nc, ios

character(len=5) :: ind

logical :: &
      lbes4     ! 4 Bessel functions expansion

real(kind=dp), parameter :: thresh=1.d-12

if (pseudotype.eq.1.or.pseudotype.eq.2) then
   write(6, &
      '(/,5x,15(''-''),'' Generating NC pseudopotential '',15(''-''),/)')
elseif (pseudotype.eq.3) then
   write(6, &
      '(/,5x,15(''-''),'' Generating US pseudopotential '',21(''-''),/)')
endif
if (pseudotype.eq.1.and.rel.eq.2) call errore('gener_pseudo', &
                  'not programmed' ,1)
!
!   compute the local potential from the all-electron potential
!
call pseudovloc
!
!   if nlcc is true compute here the core charge
!
if (nlcc) call set_rho_core
!
!   set the appropriate energies and the correspondence all-electron
!   pseudo
!
do n=1,nwfs
   if (enls(n).eq.0.d0) enls(n)=enl(nstoae(n))
enddo
!
!   compute the pseudowavefunctions by expansion in spherical
!   bessel function before r_c
!
do ns=1,nbeta
   lbes4=.false.
   lam=lls(ns)
   nst=(lam+1)*2
   nwf0=nstoae(ns)
!    
!  compute the ik closer to r_cut, r_cutus, rcloc
!
   ik=0
   ikus=0
   ikloc=0
   do n=1,mesh
      if (r(n).lt.rcut(ns)) ik=n
      if (r(n).lt.rcutus(ns)) ikus=n
      if (r(n).lt.rcloc) ikloc=n
   enddo
   if (mod(ik,2).eq.0) ik=ik+1
   if (mod(ikus,2).eq.0) ikus=ikus+1
   if (mod(ikloc,2).eq.0) ikloc=ikloc+1
   if (ikus.gt.mesh) call errore('gener_pseudo','ik is wrong ',1)
   if (pseudotype.eq.3) then
      ikk(ns)=max(ikus+10,ikloc+5)
   elseif (pseudotype.eq.1.or.pseudotype.eq.2) then
      ikk(ns)=max(ik+10,ikloc+5)
   else
      call errore('gener_pseudo','pseudotype not programmed',1)
   endif 
!
!  compute the phi functions
!
   nnode=0
   call compute_phi(lam,ik,nwf0,ns,xc,1,nnode,ocs(ns))
   if (nnode.ne.0) call errore('gener_pseudo','too many nodes',1)
!
!   US only on the components where ikus <> ik
! 
   do n=1,mesh
      psipsus(n,ns)=phis(n,ns) 
   enddo   
   if (ikus.ne.ik) then
       call compute_phius(lam,ikus,ns,xc,1)
       lbes4=.true.
   endif
   call compute_chi(lam,ns,xc,lbes4)
!
!    check that the chi are zero beyond ikk
   nst=0
   do n=1,mesh
      gi(n,1)=0.d0
   enddo
!   do n=ikk(ns)+1,min(ikk(ns)+20,mesh)
   do n=ikk(ns)+1,mesh
      gi(n,1)=chis(n,ns)**2
   enddo
   do n=min(ikk(ns)+20,mesh),mesh
        chis(n,ns)=0.d0
   enddo
   sum=int_0_inf_dr(gi,r,r2,dx,mesh,nst)
   if (sum.gt.2.d-6) then
       write(6, &
      &  '(5x,''ns='',i4,'' l='',i4, '' sum='',f15.9,'' r(ikk) '',f15.9)') &
                 ns,lam,sum,r(ikk(ns))
      call errore('gener_pseudo ','chi too large beyon r_c',-1)
      do n=ikk(ns),mesh  
         write(6,*) r(n),gi(n,1)
      enddo
      stop
   endif
enddo

!      do n=1,mesh
!         write(6,'(5e15.7)') r(n),psipsus(n,1),chis(n,1),
!     +                            psipsus(n,2),chis(n,2)
!      enddo
!      stop

!
!    for each angular momentum take the same integration point
!
do ns=1,nbeta
   do ns1=1,nbeta
      if (lls(ns).eq.lls(ns1).and.ikk(ns1).gt.ikk(ns)) &
            ikk(ns)=ikk(ns1)
   enddo
enddo
!
!     if pseudotype is 1 construct the potential
!
if (pseudotype.eq.1) then
   vnl=0.d0
   do ns=1,nbeta
      lam=lls(ns)
      do n=1,ikk(ns)
         vnl(n,lam)=chis(n,ns)/phis(n,ns)
      enddo 
   enddo
endif
!
!     construct B_{ij}
!
bmat=0.d0
do ns=1,nbeta
  do ns1=1,nbeta
    if (lls(ns).eq.lls(ns1).and.abs(jjs(ns)-jjs(ns1)).lt.1.d-7) then
      nst=(lls(ns)+1)*2
      ikl=ikk(ns1)
      do n=1,mesh
         gi(n,1)=phis(n,ns)*chis(n,ns1)
      enddo
      bmat(ns,ns1)=int_0_inf_dr(gi,r,r2,dx,ikl,nst)
    endif
  enddo
enddo

if (pseudotype.eq.2) then
!
!     symmetrize the B matrix
!
   do ns=1,nbeta
      do ns1=1,ns
         binv(ns,ns1)=0.5d0*(bmat(ns,ns1)+bmat(ns1,ns))
         binv(ns1,ns)=binv(ns,ns1)
      enddo
   enddo
   do ns=1,nbeta
      do ns1=1,nbeta
         bmat(ns,ns1)=binv(ns,ns1)
      enddo
   enddo
elseif (pseudotype.eq.3) then
!
!    find the inverse
!
   do ns=1,nbeta
      do ns1=1,nbeta
         binv(ns,ns1)=bmat(ns,ns1)
      enddo
   enddo
elseif (pseudotype.eq.1) then
   goto 500
else
   call errore('gener_pseudo','pseudotype not programmed',1)
endif
!
!   compute the inverse of the matrix B_{ij}^-1
!
call invmat1(binv,binv,nwfsx,nbeta)
!
!   compute the beta functions
!
betas=0.d0
do ns=1,nbeta
   do ns1=1,nbeta
      do n=1,mesh
         betas(n,ns)=betas(n,ns)+ binv(ns1,ns)*chis(n,ns1)
      enddo
   enddo
enddo

qq=0.d0
if (pseudotype.eq.3) then
!
!    compute the Q functions
!
  do ns=1,nbeta
    do ns1=1,ns
      ikl=max(ikk(ns),ikk(ns1))
      do n=1, ikl
         qvan(n,ns,ns1) = psipsus(n,ns) * psipsus(n,ns1) &
                        - phis(n,ns) * phis(n,ns1)
         gi(n,1)=qvan(n,ns,ns1)
      enddo
      do n=ikl+1,mesh
         qvan(n,ns,ns1)=0.d0
      enddo
!
!     and puts its integral in qq
!
      if (lls(ns).eq.lls(ns1).and.abs(jjs(ns)-jjs(ns1)).lt.1.d-8) then
         nst=(lls(ns)+1)*2
         qq(ns,ns1)=int_0_inf_dr(gi,r,r2,dx,ikk(ns),nst)
      endif
!
!     set the bmat with the eigenvalue part
!
      bmat(ns,ns1)=bmat(ns,ns1)+enls(ns1)*qq(ns,ns1)
!
!    Use symmetry of the n,ns1 indeces to set qvan and qq and bmat
!
      if (ns.ne.ns1) then
         do n=1,mesh
            qvan(n,ns1,ns)=qvan(n,ns,ns1)
         enddo
         qq(ns1,ns)=qq(ns,ns1)
         bmat(ns1,ns)=bmat(ns1,ns)+enls(ns)*qq(ns1,ns)
      endif
    enddo
  enddo
  write(6,'(/5x,'' The bmat matrix'')')
  do ns1=1,nbeta
    write(6,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
  enddo
  write(6,'(/5x,'' The qq matrix'')')
  do ns1=1,nbeta
     write(6,'(6f12.5)') (qq(ns1,ns),ns=1,nbeta)
  enddo
endif

if (pseudotype.ne.1) then
   do ib=1,nbeta
      do jb=1,nbeta
         ddd(ib,jb,1)=bmat(ib,jb)
      enddo
   enddo
endif
!
!    descreening the local potential and the D coefficients
!
500   call descreening
!
!    In pseudotype 1 the local potential is added to all channel
!
if (pseudotype.eq.1) then
   do l=0,3
      do n=1,mesh
         vnl(n,l)=vnl(n,l)+vpsloc(n)
      enddo
   enddo
   vpsloc=0.d0
endif
!
!     print the main functions on files
!
if (file_wavefunctionsps .ne. ' ') then
  open(unit=19,file=file_wavefunctionsps, status='unknown', iostat=ios,err=300)
300      call errore('gener_pseudo','opening file '//file_wavefunctionsps,&
                     abs(ios))

   do n=1,mesh
      write(19,'(i5,7e13.5)') n,r(n), (phis(n,ns), ns=1,nwfs)
   enddo
   close(19)
endif
if (file_beta .ne. ' ') then
  open(unit=19,file=file_beta, status='unknown', iostat=ios, err=400)
400    call errore('gener_pseudo','opening file '//file_beta,abs(ios))
   do n=1,mesh
      write(19,'(8f12.6)') r(n), (betas(n,ns), ns=1,nbeta)
   enddo
   close(19)
endif
if (file_chi .ne. ' ') then
  open(unit=19,file=file_chi, status='unknown', iostat=ios, err=600)
600    call errore('gener_pseudo','opening file '//file_chi,abs(ios))
   do n=1,mesh
      write(19,'(8f12.6)') r(n), (chis(n,ns), ns=1,nbeta)
   enddo
   close(19)
endif
if (file_qvan .ne. ' ') then
  do ns1=1,nbeta
     ind=' '
     if (ns1.lt.10) then
        write(ind,'(".",i1)') ns1
     elseif (ns1.lt.100) then
        write(ind,'(".",i2)') ns1
     else
        write(ind,'(".",i3)') ns1
     endif
     open(unit=19,file=TRIM(file_qvan)//TRIM(ind), status='unknown', &
                  iostat=ios, err=700)
700    call errore('gener_pseudo','opening file '//file_qvan,abs(ios))
     do n=1,mesh
        write(19,'(8f12.6)') r(n), (qvan(n,ns,ns1), ns=1,ns1)
     enddo
     close(19)
  enddo
endif

write(6, &
     & '(/,5x,12(''-''),'' End of pseudopotential generation '',20(''-''),/)')

return
end
