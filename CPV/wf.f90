#include "f_defs.h"
!-----------------------------------------------------------------------
      subroutine wf(clwf,c,bec,eigr,eigrb,taub,irb,b1,b2,b3,Uall,becdr,what1,wfc,jw,ibrav)
!-----------------------------------------------------------------------
!
!     this routine calculates overlap matrices
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
  use constants
  use ions_base, only : nsp, na, nas=>nax
  use parameters, only : natx, nsx
  !use gvec
  use gvecs
  use cvan
  use cell_base, only : omega, a1, a2, a3, alat
  use elct
  use gvecb
  use gvecw, only : ngw, ng0
  use small_box
!  use parms
  use smooth_grid_dimensions , nnrs => nnrsx
  use control_flags
  use qgb_mod
  use wfparm
  use wfparm2
  use grid_dimensions, only : nr1, nr2, nr3
  use smallbox_grid_dimensions ,nnrb => nnrbx
  use uspp_param, only: nh, nhm
  use uspp, only : nhsa=> nkb
  use parallel_include
#ifdef __PARA
      use para_mod
#endif
  implicit none
  !
  logical what1
!  integer, intent(in) :: nas
  integer, intent(in) :: irb(3,natx,nsx),jw, ibrav
!  integer, intent(in) :: irb(3,nax,nsx),jw, ibrav
!  integer :: irb(3,nax,nsx),jw, ibrav
  real(kind=8), intent(inout) :: bec(nhsa,n), becdr(nhsa,n,3)
  real(kind=8), intent(in) :: b1(3),b2(3),b3(3),taub(3,nas)
  real(kind=8) :: becwf(nhsa,n) , becdrwf(nhsa,n), temp3(nhsa,n)
  complex(kind=8), intent(inout) :: c(ngw,nx)
  complex(kind=8), intent(in) :: eigr(ngw,nas,nsp),eigrb(ngb,nas,nsp)
  complex(kind=8) :: cwf(ngw,nx), bec2(n), bec3(n), bec2up(nupdwn(1))
  complex(kind=8) :: bec2dw(nupdwn(2)), bec3up(nupdwn(1)), bec3dw(nupdwn(2))
  complex(kind=8),allocatable :: c_m(:,:),c_p(:,:),c_psp(:,:),c2(:,:)
  complex(kind=8),allocatable :: c_msp(:,:)
  real(kind=8), intent(inout) :: Uall(n,n)
  integer, allocatable :: tagz(:)
  !
  integer :: inl, jnl, iss, isa, is, ia, ijv, i, j, k, l,          &
     ig, ierr, ti, tj, tk, iv, jv, inw, iqv, ibig1, ibig2, ibig3,           &
        ir1, ir2, ir3, ir, b5(3), b6(3), clwf, m,ib,jb, total, nstat, jj
  integer :: ngpww, irb3
  real(kind=8) :: t1, t2, t3, taup(3), pi2
  real(kind=8), allocatable, dimension(:, :) :: Uspin
  complex(kind=8) :: ci, ct1, ct2, ct3, qvt, alpha, beta1
  complex(kind=8), allocatable, dimension(:, :) ::  X, X1,Xsp,X2,X3
  complex(kind=8), allocatable, dimension(:, :, :) :: O, Ospin, Oa
  real(kind=8), parameter :: autoaf=0.529177d0
  real(kind=8) alen,blen,clen

  complex(kind=8), allocatable, dimension(:) :: qv
!
!
  integer,allocatable :: f3(:), f4(:)
  real(kind=8), allocatable :: gr(:,:), mt(:), W(:,:), wr(:)
  integer :: adjust,ini, ierr1,nnn
  complex(kind=8) :: U2(n,n)
  complex(kind=8) :: boxdotgridcplx
  external boxdotgridcplx

!
#ifdef __PARA
      integer proc, ntot, ncol, mc, ngpwpp(nproc)
      integer ncol1,nz1, nz_1 
      integer nmin(3), nmax(3), n1,n2,nzx,nz,nz_
      integer nmin1(3), nmax1(3)
      integer root, rdispls(nproc), recvcount(nproc), sendcount(nproc), sdispls(nproc)
      integer rdispls2(nproc), recvcount2(nproc), sendcount2(nproc), sdispls2(nproc)
      integer   rdispls1(nproc), recvcount1(nproc), sendcount1(nproc), sdispls1(nproc)
      complex*16, allocatable:: psitot(:,:), psitot_pl(:,:)
      complex*16, allocatable:: psitot1(:), psitot_p(:)
      complex*16, allocatable:: psitot_mi(:,:)
      complex*16, allocatable:: psitot_m(:)
      integer , allocatable :: ns(:)
#else 
#endif
   


      integer igx,igy,igz
      real(kind=8) wfcx, wfcy, wfcz
      real(kind=8) wfc(3,n)

      real(kind=8) te(6)
!
  ci=(0.d0,1.d0)
  alpha=(1.d0,0.d0)
  beta1=(0.d0,0.d0)
  pi2=2.d0*pi
   if(iprsta.gt.4) then 
   write(6,*) "Now Entering wf..."
   end if

!
!
!set up the weights and the G vectors for wannie-function calculation

   te=0.d0
!
  select case (ibrav)

  case (0)
  !free cell used for cpr
     nw=6
     allocate(wfg(nw,3), weight(nw))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     call tric_wts(b1,b2,b3,alat,weight) !**** assigns weights for the triclinic/cpr case
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=1
     wfg(5,2)=1
     wfg(5,3)=1
     wfg(6,1)=1
     wfg(6,3)=1
     go to 99

  case (1)
  !cubic P [sc]
     nw=3
     allocate(tagz(nw))
     allocate(wfg(3,3), weight(3))
     weight=1.d0
     wfg=0
     tagz=1
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     tagz(3)=0
     go to 99

  case(2)
  !cubic F [fcc]
     nw=4
     allocate (tagz(nw))
     allocate(wfg(4,3), weight(4))
     weight=0.25d0
     wfg=0
     tagz=1
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,:)=-1
     tagz(3)=0
     go to 99

  case(3)
  !cubic I [bcc]
     nw=6
     allocate(wfg(6,3), weight(6))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     weight=0.25d0
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=1
     wfg(5,2)=1
     wfg(5,3)=1
     wfg(6,1)=-1
     wfg(6,3)=1
     go to 99

  case(4)
  !hexagonal and trigonal P
     nw=4
     allocate(wfg(4,3), weight(4))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     weight(1)=0.5d0
     weight(2)=0.5d0
     weight(3)=1.d0/(b3(3)**2)
     weight(4)=0.5d0
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=-1
     go to 99

  case(5)
  !trigonal R
     nw=6
     allocate(wfg(6,3), weight(6))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     t1=1.d0-1.d0/(b2(2)**2*1.5d0)
     weight(1)=1.d0
     weight(2)=1.d0+2.d0*t1
     weight(3)=1.d0
     weight(4)=t1
     weight(5)=-t1
     weight(6)=-t1
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,3)=1
     wfg(5,2)=-1
     wfg(5,3)=1
     wfg(6,1)=1
     wfg(6,2)=-1
     go to 99

  case(6)
  !tetragonal P[st]
     nw=3
     allocate(wfg(3,3), weight(3))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     weight(1)=1.d0
     weight(2)=1.d0
     weight(3)=1.d0/(b3(3)**2)
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     go to 99

  case(7)
  !tetragonal I [bct]
     nw=6
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0  
     allocate(wfg(6,3), weight(6))
     t1=0.25d0/(b2(3)**2)
     weight(1)=0.5d0-t1
     weight(2)=t1
     weight(3)=t1
     weight(4)=t1
     weight(5)=weight(1)
     weight(6)=t1
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,3)=1
     wfg(5,2)=1
     wfg(5,3)=-1
     wfg(6,1)=1
     wfg(6,2)=1
     go to 99

  case(8)
  !orthorhombic P
     nw=3
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     allocate(wfg(3,3), weight(3))
     weight(1)=1.d0
     weight(2)=1.d0/(b2(2)**2)
     weight(3)=1.d0/(b3(3)**2)
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     go to 99

  case(9)
  !one face centered orthorhombic C
     if (b1(2).eq.1) then
        nw=3
        allocate(wfg(3,3), weight(3))
        allocate(tagz(nw))
        wfg=0
        weight(1)=0.5d0
        weight(2)=0.5d0
     else
        if (b1(2).gt.1) then
           write(6, *) "Please make celldm(2) not less than 1"
#ifdef __PARA
           call mpi_finalize(i)
#endif
           stop
        else
           nw=4
           allocate(wfg(4,3), weight(4))
           allocate(tagz(nw))
           wfg=0
           weight(1)=0.5d0
           weight(2)=0.5d0
           weight(4)=(1.d0/(b1(2)**2)-1.d0)/4.d0
           wfg(4,1)=1
           wfg(4,2)=1
        end if
     end if
     weight(3)=1.d0/(b3(3)**2)
     tagz=1
     tagz(3)=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     go to 99

  case(10)
  !all face centered orthorhombic F
     if (b1(2).eq.-1.and.b1(3).eq.1) then
        write(6, *) "Please change ibrav to 2"
#ifdef __PARA
        call mpi_finalize(i)
#endif
        stop
     end if
     if ((b1(1).gt.b1(3)).or.(b1(1)+b1(2).lt.0)) then
        write(6, *) "Please make celldm(2) >= 1 >= celldm(3)"
#ifdef __PARA
        call mpi_finalize(i)
#endif
        stop
     end if
     if (b1(3).eq.1) then
        nw=5
        allocate(wfg(5,3), weight(5))
   allocate(tagz(nw))
     else
        nw=6
        allocate(wfg(6,3), weight(6))
   allocate(tagz(nw))
     end if
     weight=0.25d0/(b1(3)**2)
     wfg=0
     tagz=1
     tagz(3)=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,:)=1
     weight(5)=0.25d0/(b1(2)**2)-weight(1)
     wfg(5,2)=1
     wfg(5,3)=1
     if (b1(3).ne.1) then
        weight(6)=0.25d0-weight(1)
        wfg(6,1)=1
        wfg(6,2)=1
     end if
     go to 99

  case(11)
  !body centered orthohombic I
     if ((b1(3).eq.1).and.(b2(2).eq.1)) then
        write(6, *) "Please change ibrav to 3"
#ifdef __PARA
        call mpi_finalize(i)
#endif
        stop
     end if
     if ((b1(3).eq.1).or.(b2(2).eq.1)) then
        write(6, *) "Please change ibrav to 7"
#ifdef __PARA
        call mpi_finalize(i)
#endif
        stop
     end if
     nw=6
     allocate(wfg(6,3), weight(6))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     t1=0.25d0
     t2=t1/(b2(2)**2)
     t3=t1/(b1(3)**2)
     weight(1)=t1-t2+t3
     weight(2)=t1+t2-t3
     weight(3)=-t1+t2+t3
     weight(4)=weight(3)
     weight(5)=weight(1)
     weight(6)=weight(2)
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=1
     wfg(5,2)=1
     wfg(5,3)=1
     wfg(6,1)=-1
     wfg(6,3)=1
     go to 99
     
  case(12)
  !monoclinic P
     nw=4
     allocate(wfg(4,3), weight(4))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     t1=-b1(2)/b2(2)
     k=nint(t1)
     if ((k.eq.0).and.(t1.ge.0)) k=1
     if ((k.eq.0).and.(t1.le.0)) k=-1
     weight(4)=t1/dfloat(k)
     weight(1)=1.d0-weight(4)
     t2=t1/sqrt(1.d0-1.d0/(1.d0+b1(2)**2))
     weight(2)=t2**2-t1*k
     weight(3)=1.d0/(b3(3)**2)
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=k
     go to 99

  case(13)
  !one face centered monoclinic C
     nw=6
     allocate(wfg(6,3), weight(6))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     t1=-b1(2)/b3(2)
     k=nint(t1)
     if ((k.eq.0).and.(t1.ge.0)) k=1
     if ((k.eq.0).and.(t1.le.0)) k=-1
     t2=1.d0/(4.d0*b1(1)**2)
     t3=1.d0/(b3(2)**2)
     weight(1)=2.d0*t2*(1.d0-t1)
     weight(2)=weight(1)
     weight(3)=t3+4.d0*t2*t1*(t1-k)
     weight(4)=1.d0-2.d0*t2
     weight(5)=2.d0*t1*t2/dfloat(k)
     weight(6)=weight(5)
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=-1
     wfg(5,1)=1
     wfg(5,3)=k
     wfg(6,2)=1
     wfg(6,3)=k
     go to 99

  case default
  !free cell used for cpr
     nw=6
     allocate(wfg(nw,3), weight(nw))
     allocate(tagz(nw))
     tagz=1
     tagz(3)=0
     call tric_wts(b1,b2,b3,alat,weight) !**** assigns weights for the triclinic/cpr case
     wfg=0
     wfg(1,1)=1
     wfg(2,2)=1
     wfg(3,3)=1
     wfg(4,1)=1
     wfg(4,2)=1
     wfg(5,2)=1
     wfg(5,3)=1
     wfg(6,1)=1
     wfg(6,3)=1
     go to 99

  !ibrav=14: triclinic P

!     write(6, *) "not available to calculate wannier-function for this ibrav"
!#ifdef __PARA
!     call mpi_finalize(i)
!#endif
!     stop
  end select

!
!     set up matrix O
!
99 allocate(O(nw, n, n), X(n,n),Oa(nw, n, n))
   if(nspin.eq.2.and.nvb.gt.0) then
      allocate(X2(nupdwn(1),nupdwn(1)))
      allocate(X3(nupdwn(2),nupdwn(2)))
   end if


!      do inw=1,nw
!          do i=1,ngw
!                write(6,*) indexplus(i,inw),  tagp(i,inw)
!                write(6,*) indexminus(i,inw),  tag(i,inw)
!          end do
!        end do



#ifdef __PARA

   allocate (ns(nproc))
   do i=1,nproc
      ns(i)=0
   end do

   if(n.eq.nproc) then
     do i=1,n
        ns(i)=1
     end do
   else
     i=0
1       do j=1,nproc   
          ns(j)=ns(j)+1
          i=i+1
          if(i.ge.n) go to 2
       end do
          if(i.lt.n) go to 1 
   end if
2        if(iprsta.gt.4) then
     do j=1,nproc
        write(6,*) ns(j)
     end do
   end if

       total = 0   
   do proc=1,nproc
      ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
           total=total+ngpwpp(proc)
!           nstat=ns(proc)
        if(iprsta.gt.4) then
           write(6,*) "I am proceessor", proc, "and i have ",ns(me)," states."
        end if
   end do
   nstat=ns(me)

   allocate(psitot(total,nstat))
   allocate(psitot1(total*nstat))
   allocate(psitot_pl(total,nstat))
   allocate(psitot_p(total*nstat))
   allocate(psitot_mi(total,nstat))
   allocate(psitot_m(total*nstat))

   allocate(c_p(ngw,nx))
   allocate(c_m(ngw,nx))
   if(iprsta.gt.4) then
     write(6,*) "All allocations done"
   end if
   

   do proc=1,nproc
        sendcount(proc)=ngpwpp(me)*ns(proc)
        recvcount(proc)=ngpwpp(proc)*ns(me)
   end do
   sdispls(1)=0
   rdispls(1)=0
   
   do proc=2,nproc
      sdispls(proc)=sdispls(proc-1)+sendcount(proc-1)
      rdispls(proc)=rdispls(proc-1)+recvcount(proc-1)
   end do
!
!   Step 1. Communicate to all Procs so that each proc has all
!   G-vectors and some states instead of all states and some
!   G-vectors. This information is stored in the 1-d array 
!   psitot1.
!
   
   call MPI_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('WF','ierr<>0',ierr)
   call MPI_alltoallv(c, sendcount, sdispls, MPI_DOUBLE_COMPLEX,             &
   &             psitot1, recvcount, rdispls, MPI_DOUBLE_COMPLEX,       &
   &            MPI_COMM_WORLD,ierr)
   if (ierr.ne.0) call errore('WF','alltoallv 1',ierr)
   if(iprsta.gt.4) then
     write(6,*) "Step 1. Communicate to all Procs ... Done, wf"
   end if
#endif   
   if(clwf.eq.5) then
#ifdef __PARA
        call write_psi(c,jw)
        call MPI_finalize(ierr)
   write(6,*) "State written", jw
   STOP
   end if
#else
   do i=1,ngw
      write(22,*) c(i,jw)
   end do
   write(6,*) "State written", jw
   STOP
   end if
#endif
   
#ifdef __PARA
!
!   Step 2. Convert the 1-d array psitot1 into a 2-d array consistent with the
!   original notation c(ngw,n). Psitot contains ntot = SUM_Procs(ngw) G-vecs
!   and nstat states instead of all n states
!
   
   ngpww=0
   do proc=1,nproc
      do i=1,ns(me)
      do j=1,ngpwpp(proc)
           psitot(j+ngpww,i)=psitot1(rdispls(proc)+j+(i-1)*ngpwpp(proc))
         end do
      end do
    ngpww=ngpww+ngpwpp(proc)
   end do
   if(iprsta.gt.4) then
     write(6,*) "Step 2. Convert the 1-d array psitot1 into a 2-d array... Done, wf"
   end if

!
!   Step 3. do the translation of the 2-d array to get the transtalted
!   arrays psitot_pl and psittot_mi, corresponding to G+G' and -G+G'
!   

  do inw=1,nw   
!
!   Intermediate Check. If the translation is only along the z-direction
!   no interprocessor communication and data rearrangement is required 
!   because each processor contains all the G- components in the z-dir.
!
    if(tagz(inw).eq.0) then
        do i=1,n
          do ig=1,ngw
                if(indexplusz(ig).eq.-1) then
                   c_p(ig,i)=(0.d0,0.d0)
                else
                   c_p(ig,i)=c(indexplusz(ig),i)
                end if
                if(indexminusz(ig).eq.-1) then
                   c_m(ig,i)=(0.d0,0.d0)
                else
                     c_m(ig,i)=conjg(c(indexminusz(ig),i))
                end if
           end do
        end do
    else
      do i=1,ns(me)
   do ig=1,total
      if(indexplus(ig,inw).eq.-1) then
         psitot_pl(ig,i)=(0.d0,0.d0)
        else   
        if(tagp(ig,inw).eq.1) then
         psitot_pl(ig,i)=conjg(psitot(indexplus(ig,inw),i))
        else
         psitot_pl(ig,i)=psitot(indexplus(ig,inw),i)
        end if
      end if
      if(indexminus(ig,inw).eq.-1) then
         psitot_mi(ig,i)=(0.d0,0.d0)
        else
                  if(tag(ig,inw).eq.1) then
           psitot_mi(ig,i)=conjg(psitot(indexminus(ig,inw),i))
        else
           psitot_mi(ig,i)=psitot(indexminus(ig,inw),i)
        end if
      end if
      end do
   end do
   if(iprsta.gt.4) then
    write(6,*) "Step 3. do the translation of the 2-d array...Done, wf"
   end if
!
!   Step 4. Convert the 2-d arrays psitot_p and psitot_m into 1-d
!   arrays
!
   ngpww=0
   do proc=1,nproc
          do i=1,ns(me)
                do j=1,ngpwpp(proc)
                     psitot_p(rdispls(proc)+j+(i-1)*ngpwpp(proc))=psitot_pl(j+ngpww,i)
                     psitot_m(rdispls(proc)+j+(i-1)*ngpwpp(proc))=psitot_mi(j+ngpww,i)
                end do
           end do
        ngpww=ngpww+ngpwpp(proc)
        end do
   if(iprsta.gt.4) then
     write(6,*) "Convert the 2-d arrays psitot_p and psitot_m into 1-d arrays...Done, wf"
   end if
!
!   Step 5. Redistribute among processors. The result is stored in 2-d
!   arrays c_p and c_m consistent with the notation c(ngw,n), such that
!   c_p(j,i) contains the coefficient for c(j,i) corresponding to G+G'
!       and c_m(j,i) contains the coefficient for c(j,i) corresponding to -G+G'
!
   c_p = 0.0d0 ! call ZERO(2*ngw*nx,c_p)
        call MPI_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('WF','ierr<>0',ierr)
        call MPI_alltoallv(psitot_p, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
   &                       c_p, sendcount , sdispls, MPI_DOUBLE_COMPLEX,              &
   &                       MPI_COMM_WORLD,ierr)
   if (ierr.ne.0) call errore('WF','alltoallv 2',ierr)

   c_m = 0.0d0 ! call ZERO(2*ngw*nx,c_m)
        call MPI_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('WF','ierr<>0',ierr)
        call MPI_alltoallv(psitot_m, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
   &                       c_m, sendcount, sdispls, MPI_DOUBLE_COMPLEX,               &
   &                       MPI_COMM_WORLD,ierr)
   if (ierr.ne.0) call errore('WF','alltoallv 3',ierr)
        if(iprsta.gt.4) then
          write(6,*) "Step 5. Redistribute among processors...Done, wf"
        end if
    end if
#else
   allocate(c_p(ngw,nx))
   allocate(c_m(ngw,nx))
   do inw=1,nw
   if(tagz(inw).eq.0) then
        do i=1,n
          do ig=1,ngw
                if(indexplusz(ig).eq.-1) then
                   c_p(ig,i)=(0.d0,0.d0)
                else
                   c_p(ig,i)=c(indexplusz(ig),i)
                end if
                if(indexminusz(ig).eq.-1) then
                   c_m(ig,i)=(0.d0,0.d0)
                else
                   c_m(ig,i)=conjg(c(indexminusz(ig),i))
                end if
           end do
        end do
    else
      do i=1,n
        do ig=1,ngw
                if(indexplus(ig,inw).eq.-1) then
                   c_p(ig,i)=(0.d0,0.d0)
                else
                  if(tagp(ig,inw).eq.1) then
                   c_p(ig,i)=conjg(c(indexplus(ig,inw),i))
                  else
                   c_p(ig,i)=c(indexplus(ig,inw),i)
                  end if
                end if
                if(indexminus(ig,inw).eq.-1) then
                   c_m(ig,i)=(0.d0,0.d0)
                else
                  if(tag(ig,inw).eq.1) then
                     c_m(ig,i)=conjg(c(indexminus(ig,inw),i))
                  else
                     c_m(ig,i)=c(indexminus(ig,inw),i)
                  end if
                end if
           end do
        end do
    end if
#endif

!
!   Step 6. Calculate Overlaps
!
!   Augmentation Part first

     allocate( qv( nnrb ) )

     X=(0.d0, 0.d0)
     !
     !
     do is = 1, nvb
        do ia =1, na(is)
           ijv = 0
           do iv = 1, nh(is)
              inl = ish(is) + (iv-1)*na(is) + ia
              ijv = ijv + 1
              qv( 1 : nnrb ) = 0.0d0 
              do ig=1,ngb
                 qv(npb(ig))=eigrb(ig,ia,is)*qgb(ig,ijv,is)
                 qv(nmb(ig))=conjg(eigrb(ig,ia,is)*qgb(ig,ijv,is))
              end do
#ifdef __PARA
     irb3=irb(3,ia,is)
#endif
              call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
              iqv=1
              qvt=(0.d0,0.d0)
              qvt=boxdotgridcplx(irb(1,ia,is),qv,expo(1,inw))

#ifdef __PARA
              call reduce(2, qvt)
#endif
!
              if (nspin.eq.1) then
                 bec2(1:n)=(0.d0,0.d0)
                 bec2(1:n)=bec(inl,1:n)*alpha
                 call ZSYRK('U','T',n,1,qvt,bec2,1,alpha,X,n)
              else
                 X2=(0.d0,0.d0)
                 X3=(0.d0,0.d0)
                 bec2up(1:nupdwn(1))=(0.d0,0.d0)
                 bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))
                 call ZSYRK('U','T',nupdwn(1),1,qvt,bec2up,1,alpha,X2,nupdwn(1))
                 bec2dw(1:nupdwn(2))=(0.d0,0.d0)
                 bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):n)
                 call ZSYRK('U','T',nupdwn(2),1,qvt,bec2dw,1,alpha,X3,nupdwn(2))
                 do i = 1, nupdwn(1)
                   do j=i, nupdwn(1)
                     X(i,j)=X(i,j)+X2(i,j)
                   end do
                 end do
                 do i = 1,nupdwn(2)
                    do j=i,nupdwn(2)
                      X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                    end do
                 end do
              end if
              do jv = iv + 1, nh(is)
                 jnl = ish(is) + (jv-1)*na(is) + ia
                 ijv = ijv + 1
                 qv( 1:nnrb ) = 0.0d0 ! call zero(2*nnrb,qv)
                 do ig=1,ngb
                    qv(npb(ig))=eigrb(ig,ia,is)*qgb(ig,ijv,is)
                    qv(nmb(ig))=conjg(eigrb(ig,ia,is)*qgb(ig,ijv,is))
                 end do
                 call ivfftbold(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
                 iqv=1
                 qvt=0.d0
                 qvt=boxdotgridcplx(irb(1,ia,is),qv,expo(1,inw))
#ifdef __PARA
                 call reduce(2, qvt)
#endif
!
                 if (nspin.eq.1) then
                    bec2(1:n)=(0.d0,0.d0)
                    bec3(1:n)=(0.d0,0.d0)
                    bec2(1:n)=bec(inl,1:n)*alpha
                    bec3(1:n)=bec(jnl,1:n)*alpha
                    call ZSYR2K('U','T',n,1,qvt,bec2,1,bec3,1,alpha,X,n)
                 else
                    X2=(0.d0,0.d0)
                    X3=(0.d0,0.d0)
                    bec2up(1:nupdwn(1))=(0.d0,0.d0)
                    bec3up(1:nupdwn(1))=(0.d0,0.d0)
                    bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))*alpha
                    bec3up(1:nupdwn(1))=bec(jnl,1:nupdwn(1))*alpha
                    call ZSYR2K('U','T',nupdwn(1),1,qvt,bec2up,1,bec3up,1,alpha,X2,nupdwn(1))
                    bec2dw(1:nupdwn(2))=(0.d0,0.d0)
                    bec3dw(1:nupdwn(2))=(0.d0,0.d0)
                    bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):n)*alpha
                    bec3dw(1:nupdwn(2))=bec(jnl,iupdwn(2):n)*alpha
                    call ZSYR2K('U','T',nupdwn(2),1,qvt,bec2dw,1,bec3dw,1,alpha,X3,nupdwn(2))
                    do i = 1, nupdwn(1)
                      do j=i, nupdwn(1)
                        X(i,j)=X(i,j)+X2(i,j)
                      end do
                    end do
                    do i = 1,nupdwn(2)
                       do j=i,nupdwn(2)
                         X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                       end do
                    end do
                 end if
              end do
           end do
        end do
     end do
     t1=omega/dfloat(nr1*nr2*nr3)
     X=X*t1
     do i=1, n
        do j=i+1, n
           X(j, i)=X(i, j)
        end do
     end do
     Oa(inw, :, :)=X(:, :)
   if(iprsta.gt.4) then
     write(6,*) "Augmentation Part Done"
   end if

   deallocate( qv )

!   Then Soft Part
   if(nspin.eq.1) then
!   Spin Unpolarized calculation
               X=0.d0   
    if(ng0.eq.2) then
         c_m(1,:)=0.d0
    end if
!           cwf(:,:)=beta1
!           cwf(:,:)=c(:,:)
      call ZGEMM('c','N',n,n,ngw,alpha,c,ngw,c_p,ngw,alpha,X,n)
      call ZGEMM('T','N',n,n,ngw,alpha,c,ngw,c_m,ngw,alpha,X,n)
#ifdef __PARA
         call reduce (2*n*n,X)
#endif
              O(inw,:,:)=Oa(inw,:,:)+X(:,:)
   if(iprsta.gt.4) then
     write(6,*) "Soft Part Done"
   end if
     
   else
!   Spin Polarized case
!   Up Spin First
     allocate(Xsp(n,nupdwn(1)))
     allocate(c_psp(ngw,nupdwn(1)))
     allocate(c_msp(ngw,nupdwn(1)))
               Xsp=0.d0
          c_psp=0.d0 
          c_msp=0.d0
    do i=1,nupdwn(1)
       c_psp(:,i)=c_p(:,i)
       c_msp(:,i)=c_m(:,i)
    end do
           if(ng0.eq.2) then
              c_msp(1,:)=0.d0
           end if
!           cwf(:,:)=beta1
!           cwf(:,:)=c(:,:,1,1)
           call ZGEMM('c','N',n,nupdwn(1),ngw,alpha,c,ngw,c_psp,ngw,alpha,Xsp,n)
           call ZGEMM('T','N',n,nupdwn(1),ngw,alpha,c,ngw,c_msp,ngw,alpha,Xsp,n)
#ifdef __PARA
              call reduce (2*n*nupdwn(1),Xsp)
#endif
              do i=1,nupdwn(1)
            do j=1,n
                   X(j,i)=Xsp(j,i)
       end do
              end do
   deallocate(Xsp,c_psp,c_msp)
!    Then Down Spin
   allocate(Xsp(n,iupdwn(2):n))
        allocate(c_psp(ngw,iupdwn(2):n))
        allocate(c_msp(ngw,iupdwn(2):n))
               Xsp=0.d0
               c_psp=0.d0
               c_msp=0.d0
         do i=iupdwn(2),n
            c_psp(:,i)=c_p(:,i)
            c_msp(:,i)=c_m(:,i)
    end do
           if(ng0.eq.2) then
              c_msp(1,:)=0.d0
           end if
!           cwf(:,:)=beta1
!           cwf(:,:)=c(:,:,1,1)
           call ZGEMM('c','N',n,nupdwn(2),ngw,alpha,c,ngw,c_psp,ngw,alpha,Xsp,n)
           call ZGEMM('T','N',n,nupdwn(2),ngw,alpha,c,ngw,c_msp,ngw,alpha,Xsp,n)
#ifdef __PARA
              call reduce (2*n*nupdwn(2),Xsp)
#endif
              do i=iupdwn(2),n
                do j=1,n
                   X(j,i)=Xsp(j,i)
                end do
              end do
   deallocate(Xsp,c_psp,c_msp)
   O(inw,:,:)=Oa(inw,:,:)+X(:,:)
   end if
     end do
#ifdef __PARA
   deallocate(ns)
#endif

 if(clwf.eq.2) then
!    output the overlap matrix to fort.38
#ifdef __PARA
  if(me.eq.1) then
#endif
  rewind 38
  write(38, '(i5, 2i2, i3, f9.5)') n, nw, nspin, ibrav, alat
  if (nspin.eq.2) then
     write(38, '(i5)') nupdwn(1)
  end if
  write(38, *) a1
  write(38, *) a2
  write(38, *) a3
  write(38, *) b1
  write(38, *) b2
  write(38, *) b3
  do inw=1, nw
     write(38, *) wfg(inw, :), weight(inw)
  end do
  do inw=1, nw
     do i=1, n
        do j=1, n
           write(38, *) O(inw, i, j)
        end do
     end do
  end do
  do i=1, n
     do j=1, n
        write(38, *) Uall(i, j)
     end do
  end do
    close(38)
#ifdef __PARA
  end if
#endif
#ifdef __PARA
  call Mpi_finalize(ierr)
#endif
  STOP
  end if
   
   if(clwf.eq.3.or.clwf.eq.4) then
   if(nspin.eq.1) then
   if(.not.what1) then
        if(wfsd) then
      call wfsteep(n,O,Uall,b1,b2,b3)
        else
      call ddyn(n,O,Uall,b1,b2,b3)
        end if
   end if
   if(iprsta.gt.4) then
     write(6,*) "Out from DDYN"
   end if
   else
     allocate(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uspin(i, j)=Uall(i, j)
           Ospin(:, i, j)=O(:, i, j)
        end do
     end do
     if(.not.what1) then
     if(wfsd) then
        call wfsteep(nupdwn(1), Ospin, Uspin,b1,b2,b3)
     else 
        call ddyn(nupdwn(1), Ospin, Uspin,b1,b2,b3)
     end if
     end if
     do i=1, nupdwn(1)
        do j=1, nupdwn(1)
           Uall(i, j)=Uspin(i, j)
       O(:,i,j)  =Ospin(:,i,j)
        end do
     end do
     deallocate(Uspin, Ospin)
     allocate(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
           Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
        end do
     end do
     if(.not.what1) then
     if(wfsd) then
         call wfsteep(nupdwn(2), Ospin, Uspin,b1,b2,b3)
     else
         call ddyn(nupdwn(2), Ospin, Uspin,b1,b2,b3)
     end if
     end if
     do i=1, nupdwn(2)
        do j=1, nupdwn(2)
           Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
      O(:,i+nupdwn(1),j+nupdwn(1))=Ospin(:,i,j)
        end do
     end do
     deallocate(Uspin, Ospin)
    end if
   end if   

  !       Update C and bec
   cwf=beta1
!        cwf(:,:)=c(:,:,1,1)
   becwf=0.0
        U2=Uall*alpha
           call ZGEMM('N','N',ngw,n,n,alpha,c,ngw,U2,n,beta1,cwf,ngw)
!           call ZGEMM('N','N',ngw,n,n,alpha,cwf,ngw,U2,n,beta1,cwf,ngw)
           call DGEMM('N','N',nhsa,n,n,alpha,bec,nhsa,Uall,n,beta1,becwf,nhsa)
        U2=beta1
     if(iprsta.gt.4) then
      write(6,*) "Updating Wafefunctions and Bec"
     end if


!          do inw=1, 3
!           do i=1, nhsa
!             do j=1, n
!               temp3(i,j)=becdr(i,j,inw)
!               becdrwf(i,j)=0.d0
!                 do nnn=1,n
!                      becdrwf(i,j)=becdrwf(i,j)+temp3(i,nnn)*Uall(nnn,j)
!                  end do
!             end do
!           end do
!            do i=1,nhsa
!              do j=1,n
!                becdr(i,j,inw)=becdrwf(i,j)
!              end do
!            end do
!          end do

         c(:,:)=cwf(:,:)
         bec(:,:)=becwf(:,:)

!        do i=1,n
!           do j=1,ngw
!                c(j,i)=cwf(j,i)
!           end do
!           do k=1,nhsa
!                bec(k,i)=becwf(k,i)
!           end do
!        end do

         if(iprsta.gt.4) then
           write(6,*) "Wafefunctions and Bec Updated"
          end if

!
! calculate wannier-function centers
!
  allocate(wr(nw), W(nw, nw),gr(nw,3),f3(nw),f4(nw),mt(nw))
  do inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  end do
!
! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
! to check the correctness of choices on G vectors
!
  do i=1, nw
     do j=1, nw
        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
     end do
  end do
!  write(24, *) "wannier function centers: (unit:\AA)"
  do i=1, n
     mt=-aimag(log(O(:,i,i)))/pi2
     wfc(1, i)=SUM(mt*weight*gr(:,1))
     wfc(2, i)=SUM(mt*weight*gr(:,2))
     wfc(3, i)=SUM(mt*weight*gr(:,3))
     do inw=1, nw
        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
     end do
     mt=wr
     f3=0
     adjust=0
!
!   balance the phase factor if necessary
!
#ifdef VARIABLECELL 
#else
!     do while(SUM((mt-f3)**2).gt.0.01d0)
!        f4=f3
!        f3=nint(mt-mt(1))
!        if (adjust.gt.200) f3=f3-1
!        if (adjust.gt.100.and.adjust.le.200) f3=f3+1
!        mt=wr+matmul(W, f3)
!    if(iprsta.gt.4) then
!        write(6,*) "mt:", mt
!        write(6,*) "f3:", f3
!    end if
!        adjust=adjust+1
!        if (adjust.gt.300) stop "unable to balance the phase!"
!     end do
#endif
     wfc(1,i)=(wfc(1,i)+SUM(mt*weight*gr(:,1)))*alat
     wfc(2,i)=(wfc(2,i)+SUM(mt*weight*gr(:,2)))*alat
     wfc(3,i)=(wfc(3,i)+SUM(mt*weight*gr(:,3)))*alat
  end do

!  if (ibrav.eq.1.or.ibrav.eq.6.or.ibrav.eq.8) then
!     do i=1, n
!        if (wfc(1, i).lt.0) wfc(1, i)=wfc(1, i)+a1(1)
!        if (wfc(2, i).lt.0) wfc(2, i)=wfc(2, i)+a2(2)
!        if (wfc(3, i).lt.0) wfc(3, i)=wfc(3, i)+a3(3)
!     end do
!  end if
!  if(what1) then
!  alen=a1(1)+a2(1)+a3(1)
!  blen=a1(2)+a2(2)+a3(2)
!  clen=a1(3)+a2(3)+a3(3)
!  do i=1, n
!        do while (wfc(1,i).lt.0.)
!           wfc(1,i)=wfc(1,i)+alen
!        end do
!        do while (wfc(2,i).lt.0.)
!           wfc(2,i)=wfc(2,i)+blen
!        end do
!        do while (wfc(3,i).lt.0.)
!           wfc(3,i)=wfc(3,i)+clen
!        end do
!        do while (wfc(1,i).gt.alen)
!           wfc(1,i)=wfc(1,i)-alen
!        end do
!        do while (wfc(2,i).gt.blen)
!           wfc(2,i)=wfc(2,i)-blen
!        end do
!        do while (wfc(3,i).gt.clen)
!           wfc(3,i)=wfc(3,i)-clen
!        end do
!      end do
! end if
#ifdef __PARA
        if(me.eq.1) then
#endif
  if(.not.what1) then
  do i=1,n
     write(26, '(3f11.6)') wfc(:,i)*autoaf
  end do
  end if
#ifdef __PARA
        end if
#endif
   if(nspin.eq.2.and.nvb.gt.0) then
      deallocate(X2,X3)
   end if
   deallocate(wr, W,mt,f3,f4,gr)
     if(iprsta.gt.4) then
   write(6,*) "deallocated wr, w, f3, f4, gr"
     end if
#ifdef __PARA
        deallocate (psitot)
        deallocate (psitot1)
    if(iprsta.gt.4) then
        write(6,*) "deallocated psitot, psitot1"
     end if
        deallocate (psitot_pl)
        deallocate (psitot_p)
    if(iprsta.gt.4) then
        write(6,*) "deallocated psitot_pl,psuitot_p"
     end if
        deallocate (psitot_mi)
        deallocate (psitot_m)
    if(iprsta.gt.4) then
        write(6,*) "deallocated psitot_mi,psitot_m"
     end if
        deallocate(c_p)
        deallocate(c_m)
    if(iprsta.gt.4) then
        write(6,*) "deallocated c_p, c_m"
     end if
#else
        deallocate(c_p,c_m)
    if(iprsta.gt.4) then
        write(6,*) "deallocated c_p,c_m"
     end if
#endif
   deallocate(X,O,Oa,wfg,weight,tagz)
    if(iprsta.gt.4) then
        write(6,*) "deallocated X,O,Oa,wfg,weight,tagz"
     end if
!  deallocate(Oac,Oa,Osc,Os,mr1,mc1)
   if(iprsta.gt.4) then
       write(6,*) "Now Leaving subroutine wf..."
   end if

  return
end subroutine wf
!------------------------------------------------------------------------
   subroutine ddyn(m,Omat,Umat,b1,b2,b3)
!  This part of the subroutine wf has been added by Manu. It performes
!  Damped Dynamics on the A matrix to get the Unitary transformation to
!  obtain the wannier function at time(t+delta). It also updates the
!  quantities bec and becdr
!  
!                                                      MANU
!                                                      November 26, 2001
!--------------------------------------------------------------------------
 
  use wfparm2
  use wfparm
  use cell_base
  use constants
  use elct
  use control_flags, only: iprsta
  use parallel_include
#ifdef __PARA
  use para_mod
#endif
  implicit none

  integer :: f3(nw), f4(nw), i,j,inw
  integer ,intent(in) :: m
  real(kind=8), intent(in) :: b1(3),b2(3),b3(3)
  real(kind=8), intent(inout) :: Umat(m,m)
  complex(kind=8), intent(inout) :: Omat(nw,m,m)
  complex(kind=8) :: U2(m,m),U3(m,m)
  integer :: adjust,ini, ierr1,nnn
  real(kind=8), allocatable, dimension(:) :: wr
  real(kind=8), allocatable, dimension(:,:) :: W
  real(kind=8) :: t0, fric,U(m,m), t2
  real(kind=8) :: A(m,m),oldt0,Wm(m,m),U1(m,m)
  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(4*m)
!  real(kind=8) :: Aminus(m,m), Aplus(m,m),f2(4*m)
  real(kind=8) :: temp(m,m)
  complex(kind=8) :: d(m,m), alpha, beta1, ci
  complex(kind=8) :: f1(2*m-1), wp(m*(m+1)/2),z(m,m)
  complex(kind=8), allocatable, dimension(:, :) :: X1
  complex(kind=8), allocatable, dimension(:, :, :) :: Oc
  real(kind=8) , allocatable , dimension(:) :: mt
  real(kind=8), parameter :: autoaf=0.529177d0
  real(kind=8) :: spread, sp, pi2
  real(kind=8) :: wfc(3,n), gr(nw,3)

  alpha=(1.d0,0.d0)
  beta1=(0.d0,0.d0)
  ci   =(0.d0,1.d0)
  pi2  =2.d0*pi

  allocate(mt(nw))
  allocate(X1(m,m))
  allocate(Oc(nw,m,m))

   fric=friction
  allocate (W(m,m),wr(m))

   Umat=0.d0
   do i=1,m
       Umat(i,i)=1.d0
   end do
   
   U2=Umat*alpha

!
! update Oc using the initial guess of Uspin
!
  do inw=1, nw
    X1(:, :)=Omat(inw, :, :)
     U3=beta1
!    call ZGEMUL(U2, m, 'T', X1, m, 'N', U3, m, m,m,m) 
     call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
!    call ZGEMUL(U3, m, 'N', U2, m, 'N', X1, m, m,m,m) 
     call ZGEMM ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  end do

   U2=beta1
   U3=beta1

    oldt0=0.d0
    A=0.d0
    Aminus=A
    temp=Aminus


!   START ITERATIONS HERE

  do ini=1, nsteps

    t0=0.d0     !use t0 to store the value of omega
    do inw=1, nw
       do i=1, m
          t0=t0+real(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       end do
    end do

        if(ABS(t0-oldt0).lt.tolw) then
#ifdef __PARA
   if(me.eq.1) then
#endif
      write(27,*) "MLWF Generated at Step",ini
#ifdef __PARA
   end if
#endif
   if(iprsta.gt.4) then
      write(6,*) "MLWF Generated at Step",ini
   end if
      go to 241
   end if

        if(adapt) then
   if(oldt0.lt.t0) then
       fric=fric/2.
       A=Aminus
       Aminus=temp
   end if
   end if

!   calculate d(omega)/dA and store result in W
!   this is the force for the damped dynamics
!

    W=0.d0
    do inw=1, nw
       t2=weight(inw)
       do i=1,m
          do j=1,m
             W(i,j)=W(i,j)+t2*real(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
                  -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          end do
       end do
    end do
   

!   the verlet scheme to calculate A(t+dt)
   
   Aplus=0.d0

   do i=1,m
     do j=i+1,m
    Aplus(i,j)=Aplus(i,j)+(2*dt/(2*dt+fric))*(2*A(i,j)               &
         -Aminus(i,j)+(dt*dt/q)*W(i,j)) + (fric/(2*dt+fric))*Aminus(i,j)
     enddo
   enddo

   Aplus=Aplus-transpose(Aplus)
   Aplus=(Aplus-A)

    do i=1, m
       do j=i,m 
        wp(i + (j-1)*j/2) = cmplx(0.0, Aplus(i,j))
       end do
    end do

#if ! defined __AIX
    call zhpev('V','U',m,wp,wr,z,m,f1,f2,ierr1)
#else
    call zhpev(21, wp, wr, z, m, m, f2, 4*m)
    ierr1 = 0
#endif

    if (ierr1.ne.0) then 
   write(6,*) "failed to diagonalize W!"
    stop
    end if

    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*wr(i)*dt)
    end do      !d=exp(d)

!   U=z*exp(d)*z+
!   
     U3=beta1
!    call ZGEMUL(z, m, 'N', d, m, 'N', U3, m, m,m,m)
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)  
     U2=beta1
!    call ZGEMUL(U3, m, 'N', z, m, 'c', U2, m, m,m,m)
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
    U=real(U2)
    U2=beta1
    U3=beta1

   temp=Aminus
   Aminus=A
   A=Aplus


!   update Umat
!
!    call DGEMUL(Umat, m, 'N', U, m, 'N', U1, m, m,m,m) 
      U1=beta1
     call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
   
    Umat=U1 

!   update Oc
!
    U2=Umat*alpha
    U3=beta1
  do inw=1, nw
    X1(:, :)=Omat(inw, :, :)
!    call ZGEMUL(U2, m, 'T', X1, m, 'N', U3, m, m,m,m)
     call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
    X1=beta1
!    call ZGEMUL(U3, m, 'N', U2, m, 'N', X1, m, m,m,m)
     call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
    Oc(inw, :, :)=X1(:, :)
  end do
    U2=beta1
    U3=beta1

   if(ABS(t0-oldt0).ge.tolw.and.ini.ge.nsteps) then
#ifdef __PARA
   if(me.eq.1) then
#endif 
      write(27,*) "MLWF Not generated after",ini,"Steps." 
#ifdef __PARA
   end if
#endif
   if(iprsta.gt.4) then
      write(6,*) "MLWF Not generated after",ini,"Steps." 
   end if
           go to 241
        end if

    oldt0=t0

   end do

241  deallocate(wr, W)
  spread=0.0
    do i=1, m
       mt=1.d0-real(Oc(:,i,i)*conjg(Oc(:,i,i)))
       sp= (alat*autoaf/pi2)**2*SUM(mt*weight)
#ifdef __PARA
       if(me.eq.1) then
#endif
       write(25, '(f10.7)') sp
#ifdef __PARA
       end if
#endif
     if(sp.lt.0.d0) then
      write(6,*) "Something wrong WF Spread negative. The Program will Stop."
#ifdef __PARA
      call MPI_FINALIZE(ierr1)
#endif
      STOP
     end if
     spread=spread+sp
    end do
    spread=spread/m

#ifdef __PARA
   if(me.eq.1) then
#endif
    write(24, '(f10.7)') spread
    write(27,*) "Average spread = ", spread
#ifdef __PARA
      end if
#endif
    Omat=Oc
   if(iprsta.gt.4) then
           write(6,*) "Average spread = ", spread
   end if

!
! calculate wannier-function centers
!
!  allocate(wr(nw), W(nw, nw))
!  do inw=1, nw
!     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
!  end do
!
! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
! to check the correctness of choices on G vectors
!
!  do i=1, nw
!     do j=1, nw
!        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
!     end do
!  end do
!  write(24, *) "wannier function centers: (unit:\AA)"
!  do i=1, m
!     mt=-aimag(log(Oc(:,i,i)))/pi2
!     wfc(1, i)=SUM(mt*weight*gr(:,1))
!     wfc(2, i)=SUM(mt*weight*gr(:,2))
!     wfc(3, i)=SUM(mt*weight*gr(:,3))
!     do inw=1, nw
!        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
!     end do
!     mt=wr
!     f3=0
!     adjust=0
!
!   balance the phase factor if necessary
!
!     do while(SUM((mt-f3)**2).gt.0.01d0)
!        f4=f3
!        f3=nint(mt-mt(1))
!        if (adjust.gt.200) f3=f3-1
!        if (adjust.gt.100.and.adjust.le.200) f3=f3+1
!        mt=wr+matmul(W, f3)
!!        write(6,*) "mt:", mt
!!        write(6,*) "f3:", f3
!        adjust=adjust+1
!        if (adjust.gt.300) stop "unable to balance the phase!"
!     end do
!     wfc(1,i)=(wfc(1,i)+SUM(mt*weight*gr(:,1)))*alat
!     wfc(2,i)=(wfc(2,i)+SUM(mt*weight*gr(:,2)))*alat
!     wfc(3,i)=(wfc(3,i)+SUM(mt*weight*gr(:,3)))*alat
!  end do
!
!  if (ibrav.eq.1.or.ibrav.eq.6.or.ibrav.eq.8) then
!     do i=1, m
!        if (wfc(1, i).lt.0) wfc(1, i)=wfc(1, i)+a1(1)
!        if (wfc(2, i).lt.0) wfc(2, i)=wfc(2, i)+a2(2)
!        if (wfc(3, i).lt.0) wfc(3, i)=wfc(3, i)+a3(3)
!     end do
!  end if
!#ifdef __PARA
!   if(me.eq.1) then
!#endif
!  if(.not.what1) then
!  do i=1, m
!     write(26, '(3f11.6)') wfc(:,i)*autoaf
!  end do
!  end if
!#ifdef __PARA
!   end if
!#endif
!   deallocate(wr, W)
   deallocate (mt,X1,Oc)
   if(iprsta.gt.4) then
   write(6,*) "Leaving DDYN"
   end if
    return
   end subroutine ddyn
!-----------------------------------------------------------------------
 subroutine wfunc_init(clwf,b1,b2,b3,ibrav)
!-----------------------------------------------------------------------

   use gvec, only: gx, mill_l
   use gvecw, only : ngw, ng0
   use elct
   use wfparm
!   use cell_base
        use cvan
   use parallel_include     
#ifdef __PARA
   use para_mod
#endif
   implicit none   
        real(kind=8), intent(in) :: b1(3),b2(3),b3(3)
#ifdef __PARA
   integer :: ntot, proc, ierr, root, i,j,inw,ngppp(nproc)
   integer :: ii,ig,recvcount(nproc), sendcount(nproc),displs(nproc)
#else
   integer :: ierr, i,j,inw, ntot
   integer :: ii,ig
#endif
        real (kind=8), allocatable:: bigg(:,:)
        integer, allocatable :: bign(:,:)
   integer :: igcount,nw1,jj,nw2, in, kk, ibrav
   integer, allocatable :: i_1(:), j_1(:), k_1(:)
   real(kind=8) :: ti, tj, tk, t1, vt, err1, err2, err3
   integer :: ti1,tj1,tk1, clwf

# ifdef __PARA
   if(n.lt.nproc) then
      write(6,*) "Number of Processors (",nproc,") is greater than the number of states (",n,")."
      write(6,*) "The Program will Stop."

   call MPI_FINALIZE(ierr)
   STOP
   end if
#endif
   allocate(gnx(3,ngw))
   allocate(gnn(3,ngw))
#ifdef __PARA
   root=0
#endif
   vt=1.0d-4
   j=0
   do i=1,ngw
             gnx(1,i)=gx(i,1)
             gnx(2,i)=gx(i,2)
             gnx(3,i)=gx(i,3)
        gnn(1,i)=mill_l(1,i)
        gnn(2,i)=mill_l(2,i)
        gnn(3,i)=mill_l(3,i)
   end do
#ifdef __PARA
   ntot=0
   do i=1,nproc
         ngppp(i)=(dfftp%nwl(i)+1)/2
   end do
   
   do proc=1,nproc
      recvcount(proc)=ngppp(proc)*3
      if(proc.eq.1) then
         displs(proc)=0
      else
         displs(proc)=displs(proc-1)+recvcount(proc-1)
      end if
      ntot=ntot+recvcount(proc)/3
   end do
   
   if(me.eq.1) then
      allocate(bigg(3,ntot))
      allocate(bign(3,ntot))
   end if
#else
   ntot=ngw
   allocate(bigg(3,ntot))
   allocate(bign(3,ntot))
   bigg(1:3,1:ntot)=gnx(1:3,1:ntot)
   bign(1:3,1:ntot)=gnn(1:3,1:ntot)
!     do j=1,ntot
!      write(50,'(6(2x,i6))') gnn(:,j), bign(:,j)
!     end do
#endif
   
   select case(ibrav)
        case(0)
!       free cell for cpr

           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminus(ntot,nw1))
           allocate(indexminusz(ngw))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=1     ! 1
           k_1(4)=0     ! 0

           i_1(5)=0     ! 0
           j_1(5)=1     ! 1
           k_1(5)=1     ! 1

           i_1(6)=1     ! 1
           j_1(6)=0     ! 0
           k_1(6)=1     ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(1)
!   cubic P [sc]

           nw1=3
           write(6,*) "Translations to be done", nw1
      allocate(indexplus(ntot,nw1))
      allocate(indexminus(ntot,nw1))
      allocate(tag(ntot,nw1))
      allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
      allocate(i_1(nw1))
      allocate(j_1(nw1))
      allocate(k_1(nw1))

           i_1(1)=1   ! 1
           j_1(1)=0   ! 0 
           k_1(1)=0   ! 0

           i_1(2)=0   ! 0 
           j_1(2)=1   ! 1
           k_1(2)=0   ! 0

           i_1(3)=0   ! 0
           j_1(3)=0   ! 0
           k_1(3)=1   ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

   go to 99

   case(2)
!   cubic F [FCC]
           nw1=4
           write(6,*) "Translations to be done", nw1
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))



           i_1(1)=1   ! 1
           j_1(1)=0   ! 0
           k_1(1)=0   ! 0

           i_1(2)=0   ! 0
           j_1(2)=1   ! 1
           k_1(2)=0   ! 0

           i_1(3)=0   ! 0
           j_1(3)=0   ! 0
           k_1(3)=1   ! 1 
   
      i_1(4)=-1    ! -1
      j_1(4)=-1    ! -1
      k_1(4)=-1    ! -1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(3)
!       cubic I [bcc]
           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))



           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

      i_1(3)=0   ! 0
      j_1(3)=0   ! 0
      k_1(3)=1   ! 1
   
      i_1(4)=1   ! 1
      j_1(4)=1   ! 1
      k_1(4)=0   ! 0   

      i_1(5)=1   ! 1
      j_1(5)=1   ! 1
      k_1(5)=1   ! 1
   
      i_1(6)=-1   ! -1
      j_1(6)=0   ! 0
      k_1(6)=1   ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

   go to 99

   case(4)
!   hexagonal and trigonal P

           nw1=4
           write(6,*) "Translations to be done", nw1
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))



           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=-1    ! -1
           k_1(4)=0     ! 0

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(5)
!   trigonal R

          nw1=6   
           write(6,*) "Translations to be done", nw1
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))



           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0
   
      i_1(3)=0   ! 0
      j_1(3)=0   ! 0
      k_1(3)=1   ! 1

      i_1(4)=1   ! 1
      j_1(4)=0   ! 0
      k_1(4)=1   ! 1
   
      i_1(5)=0   ! 0
      j_1(5)=-1   ! -1
      k_1(5)=1   ! 1

      i_1(6)=1   ! 1
      j_1(6)=-1   ! -1
      k_1(6)=0   ! 0

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(6)
! tetragonal P [st]

           nw1=3
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))

           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
      k_1(3)=1   ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(7)
! tetragonal I [bct]

           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

      i_1(4)=1   ! 1
      j_1(4)=0   ! 0
      k_1(4)=1   ! 1

      i_1(5)=0   ! 0
      j_1(5)=1   ! 1
      k_1(5)=-1   ! -1

      i_1(6)=1   ! 1
      j_1(6)=1   ! 1
      k_1(6)=0   ! 0

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(8)
! Orthorhombic P
           nw1=3
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
      indexplus=0
      indexminus=0
      tag=0
      tagp=0
      indexplusz=0
      indexminusz=0

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(9)
! One face centered Orthorhombic C
   if(b1(2).eq.1) then 

      nw1=3
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
   
   else
     if(b1(2).gt.1) then
        write (6,*) "Please make celldm(2) not less than 1"
#ifdef __PARA
   call mpi_finalize(i)
#endif
        STOP
     else
           nw1=4
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))

           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

      i_1(4)=1
      j_1(4)=1
      k_1(4)=0
     end if
   end if
    
           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0
   
   go to 99
   
   case(10)
! all face centered orthorhombic F
   
   if(b1(2).eq.-1.and.b1(3).eq.1) then
      write(6,*) "Please change ibrav to 2"
#ifdef __PARA
   call mpi_finalize(i)
#endif
      STOP
   end if
   if((b1(1).gt.b1(3)).or.(b1(1)+b1(2).lt.0)) then
      write(6,*) "Please make celldm(2)>=1>=celldm(3)"
#ifdef __PARA
   call mpi_finalize(i)
#endif
      STOP
   end if
   if(b1(3).eq.1) then
      nw1=5
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))

           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
   else
      nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))

           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))
   
      i_1(6)=1   ! 1
      j_1(6)=1   ! 1
      k_1(6)=0   ! 0
   
   end if

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1
   
      i_1(4)=1   ! 1
      j_1(4)=1    ! 1
      k_1(4)=1   ! 1

      i_1(5)=0   ! 0
      j_1(5)=1   ! 1
      k_1(5)=0   ! 1
   
           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

   go to 99

   case(11) 
! Body centered orthorhombic I
   
   if((b1(3).eq.1).and.(b2(2).eq.1)) then
      write(6,*) "Please change ibrav to 3"
#ifdef __PARA
   call mpi_finalize(i)
#endif
   STOP
   end if
   if((b1(3).eq.1).or.(b2(2).eq.1)) then
       write(6,*) "Please change ibrav to 7"
#ifdef __PARA
   call mpi_finalize(i)
#endif
   STOP
   end if

           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))

           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=1     ! 1
           k_1(4)=0     ! 0

           i_1(5)=1     ! 1
           j_1(5)=0     ! 0
           k_1(5)=1     ! 1

           i_1(6)=-1    ! -1
           j_1(6)=0     ! 0
           k_1(6)=1     ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case(12)
! monoclinic P

           nw1=4
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

      t1=-b1(2)/b2(2)
      kk=nint(t1)
      if((kk.eq.0).and.(t1.ge.0)) kk=1
      if((kk.eq.0).and.(t1.le.0)) kk=-1

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=kk    ! kk
           k_1(4)=0     ! 0

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

   go to 99

   case(13)
! one face centered monoclinic C

           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexminus(ntot,nw1))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminusz(ngw))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           t1=-b1(2)/b3(2)
           kk=nint(t1)
           if((kk.eq.0).and.(t1.ge.0)) kk=1
           if((kk.eq.0).and.(t1.le.0)) kk=-1

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=-1    ! -1
           k_1(4)=0     ! 0

      i_1(5)=1   ! 1
      j_1(5)=0   ! 0
      k_1(5)=kk   ! kk

      i_1(6)=0   ! 0
      j_1(6)=1   ! 1
      k_1(6)=kk   ! kk

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0

        go to 99

   case default
!       free cell for cpr

           nw1=6
           write(6,*) "Translations to be done", nw1
           allocate(indexplus(ntot,nw1))
           allocate(indexplusz(ngw))
           allocate(indexminus(ntot,nw1))
           allocate(indexminusz(ngw))
           allocate(tag(ntot,nw1))
           allocate(tagp(ntot,nw1))
           allocate(i_1(nw1))
           allocate(j_1(nw1))
           allocate(k_1(nw1))

           i_1(1)=1     ! 1
           j_1(1)=0     ! 0
           k_1(1)=0     ! 0

           i_1(2)=0     ! 0
           j_1(2)=1     ! 1
           k_1(2)=0     ! 0

           i_1(3)=0     ! 0
           j_1(3)=0     ! 0
           k_1(3)=1     ! 1

           i_1(4)=1     ! 1
           j_1(4)=1     ! 1
           k_1(4)=0     ! 0

           i_1(5)=0     ! 0
           j_1(5)=1     ! 1
           k_1(5)=1     ! 1

           i_1(6)=1     ! 1
           j_1(6)=0     ! 0
           k_1(6)=1     ! 1

           indexplus(:,3)=0
           indexminus(:,3)=0
           tag(:,3)=0
           tagp(:,3)=0
           go to 99

! ibrav 14 : Triclinic P
!   write(6,*) "Hey!! I'm not superman... Implement triclinic yourself"
!#ifdef __PARA
!   call mpi_finalize(i)
!#endif
!   STOP
   end select
   
99   write(6,*) "ibrav selected:", ibrav
        if(nvb.gt.0) call small_box_wf(i_1, j_1, k_1, nw1)
#ifdef __PARA
   call mpi_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_barrier' , ierr)
   call mpi_gatherv(gnx, recvcount(me), MPI_REAL8,              &
                    bigg, recvcount, displs, MPI_REAL8,         &
          root, MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_gatherv' , ierr)


        call mpi_barrier(MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) call errore('wfunc_init','mpi_barrier' , ierr)
        call mpi_gatherv(gnn, recvcount(me), MPI_INTEGER,              &
                         bign, recvcount, displs, MPI_INTEGER,         &
                         root, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) call errore('wfunc_init','mpi_gatherv' , ierr)
#endif
#ifdef __PARA
   if(me.eq.1) then
#endif
   if(clwf.eq.5) then
#ifdef __PARA
     do ii=1,ntot
      write(21,*) bigg(:,ii)
     end do
#else
     do ii=1,ngw
       write(21,*) gx(ii,1), gx(ii,2), gx(ii,3)
     end do
#endif
     close(21)
   end if
#ifdef __PARA 
   end if
#endif

   do inw=1,nw1
    if(i_1(inw).eq.0.and.j_1(inw).eq.0) then
          do ig=1,ngw
            if(ng0.eq.2) then
             indexminusz(1)=-1
            end if
           ti=(gnn(1,ig)+i_1(inw))*b1(1)+(gnn(2,ig)+j_1(inw))*b2(1)+(gnn(3,ig)+k_1(inw))*b3(1)
           tj=(gnn(1,ig)+i_1(inw))*b1(2)+(gnn(2,ig)+j_1(inw))*b2(2)+(gnn(3,ig)+k_1(inw))*b3(2)
           tk=(gnn(1,ig)+i_1(inw))*b1(3)+(gnn(2,ig)+j_1(inw))*b2(3)+(gnn(3,ig)+k_1(inw))*b3(3)
           do ii=1,ngw
        err1=ABS(gnx(1,ii)-ti)
        err2=ABS(gnx(2,ii)-tj)
        err3=ABS(gnx(3,ii)-tk)
!             if(gnx(1,ii).eq.ti.and.gnx(2,ii).eq.tj.and.gnx(3,ii).eq.tk) then
             if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                indexplusz(ig)=ii
!               write (6,*) "Found +", ig,ii,inw
                go to 224
             else
             end if
           end do
           indexplusz(ig)=-1
!                write (6,*) "Not Found +", ig,-1,inw
224        ti=(-gnn(1,ig)+i_1(inw))*b1(1)+(-gnn(2,ig)+j_1(inw))*b2(1)+(-gnn(3,ig)+k_1(inw))*b3(1)
           tj=(-gnn(1,ig)+i_1(inw))*b1(2)+(-gnn(2,ig)+j_1(inw))*b2(2)+(-gnn(3,ig)+k_1(inw))*b3(2)
           tk=(-gnn(1,ig)+i_1(inw))*b1(3)+(-gnn(2,ig)+j_1(inw))*b2(3)+(-gnn(3,ig)+k_1(inw))*b3(3)
           ti1=-gnn(1,ig)+i_1(inw)
           tj1=-gnn(2,ig)+j_1(inw)
           tk1=-gnn(3,ig)+k_1(inw)
           if(ti1.lt.0.or.(ti1.eq.0.and.(tj1.lt.0.or.(tj1.eq.0.and.tk1.lt.0)))) then
                do ii=1,ngw
                   err1=ABS(gnx(1,ii)+ti)
                   err2=ABS(gnx(2,ii)+tj)
                   err3=ABS(gnx(3,ii)+tk)
!                   if(gnx(1,ii).eq.-ti.and.gnx(2,ii).eq.-tj.and.gnx(3,ii).eq.-tk) then
                    if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                     indexminusz(ig)=ii
!                     tag(ig,inw)=1
!                     write (6,*) "Found -", ig,ii,inw
                     go to 223
                    else
                    end if
                end do
                indexminusz(ig)=-1
!                tag(ig,inw)=1
!                write (6,*) "Not Found -", ig,-1,inw
           else
              do ii=1,ngw
                   err1=ABS(gnx(1,ii)-ti)
                   err2=ABS(gnx(2,ii)-tj)
                   err3=ABS(gnx(3,ii)-tk)
!                 if(gnx(1,ii).eq.ti.and.gnx(2,ii).eq.tj.and.gnx(3,ii).eq.tk) then
                  if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                   indexminusz(ig)=ii
!                   tag(ig,inw)=-1
!                   write (6,*) "Found -", ig,ii,inw
                   go to 223
                 else
                 end if
              end do
              indexminusz(ig)=-1
!              tag(ig,inw)=-1
!              write (6,*) "Not Found -", ig,-1,inw
          end if
223     continue
   end do
        write(6,*) "Translation", inw, "for", ngw, "G vectors"
      else
#ifdef __PARA
    if(me.eq.1) then   
#endif
     do ig=1,ntot
       if(ng0.eq.2) then
             indexminus(1,inw)=-1
           end if
           ti=(bign(1,ig)+i_1(inw))*b1(1)+(bign(2,ig)+j_1(inw))*b2(1)+(bign(3,ig)+k_1(inw))*b3(1)
           tj=(bign(1,ig)+i_1(inw))*b1(2)+(bign(2,ig)+j_1(inw))*b2(2)+(bign(3,ig)+k_1(inw))*b3(2)
           tk=(bign(1,ig)+i_1(inw))*b1(3)+(bign(2,ig)+j_1(inw))*b2(3)+(bign(3,ig)+k_1(inw))*b3(3)
      ti1=bign(1,ig)+i_1(inw)
      tj1=bign(2,ig)+j_1(inw)
      tk1=bign(3,ig)+k_1(inw)
          if(ti1.lt.0.or.(ti1.eq.0.and.(tj1.lt.0.or.(tj1.eq.0.and.tk1.lt.0)))) then
           do ii=1,ntot
             err1=ABS(bigg(1,ii)+ti)
             err2=ABS(bigg(2,ii)+tj)
             err3=ABS(bigg(3,ii)+tk)
              if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
!             if(bigg(1,ii).eq.-ti.and.bigg(2,ii).eq.-tj.and.bigg(3,ii).eq.-tk) then
                indexplus(ig,inw)=ii
      tagp(ig,inw)=1
!                write (6,*) "Found +", ig,ii,inw 
                go to 214
             else
             end if
           end do
           indexplus(ig,inw)=-1
      tagp(ig,inw)=1
!          write (6,*) "Not Found +", ig,-1,inw 
     else
           do ii=1,ntot
             err1=ABS(bigg(1,ii)-ti)
             err2=ABS(bigg(2,ii)-tj)
             err3=ABS(bigg(3,ii)-tk)
              if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
!             if(bigg(1,ii).eq.ti.and.bigg(2,ii).eq.tj.and.bigg(3,ii).eq.tk) then
                indexplus(ig,inw)=ii
                tagp(ig,inw)=-1
!                write (6,*) "Found +", ig,ii,inw
                go to 214
             else
             end if
           end do
           indexplus(ig,inw)=-1
      tagp(ig,inw)=-1
!          write (6,*) "Not Found +", ig,-1,inw
     end if
214        ti=(-bign(1,ig)+i_1(inw))*b1(1)+(-bign(2,ig)+j_1(inw))*b2(1)+(-bign(3,ig)+k_1(inw))*b3(1)
           tj=(-bign(1,ig)+i_1(inw))*b1(2)+(-bign(2,ig)+j_1(inw))*b2(2)+(-bign(3,ig)+k_1(inw))*b3(2)
           tk=(-bign(1,ig)+i_1(inw))*b1(3)+(-bign(2,ig)+j_1(inw))*b2(3)+(-bign(3,ig)+k_1(inw))*b3(3)
      ti1=-bign(1,ig)+i_1(inw)
      tj1=-bign(2,ig)+j_1(inw)
      tk1=-bign(3,ig)+k_1(inw)
           if(ti1.lt.0.or.(ti1.eq.0.and.(tj1.lt.0.or.(tj1.eq.0.and.tk1.lt.0)))) then
                do ii=1,ntot
                    err1=ABS(bigg(1,ii)+ti)
                    err2=ABS(bigg(2,ii)+tj)
                    err3=ABS(bigg(3,ii)+tk)
                    if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
!                   if(bigg(1,ii).eq.-ti.and.bigg(2,ii).eq.-tj.and.bigg(3,ii).eq.-tk) then
                     indexminus(ig,inw)=ii
                tag(ig,inw)=1
!                     write (6,*) "Found -", ig,ii,inw 
                     go to 213
                    else
                    end if
                end do
                indexminus(ig,inw)=-1
           tag(ig,inw)=1
!                write (6,*) "Not Found -", ig,-1,inw 
           else 
              do ii=1,ntot
                 err1=ABS(bigg(1,ii)-ti)
                 err2=ABS(bigg(2,ii)-tj)
                 err3=ABS(bigg(3,ii)-tk)
                 if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
!                 if(bigg(1,ii).eq.ti.and.bigg(2,ii).eq.tj.and.bigg(3,ii).eq.tk) then
                   indexminus(ig,inw)=ii
              tag(ig,inw)=-1
!                   write (6,*) "Found -", ig,ii,inw 
                   go to 213
                 else
                 end if
              end do
              indexminus(ig,inw)=-1
         tag(ig,inw)=-1
!              write (6,*) "Not Found -", ig,-1,inw 
          end if
213     continue
        end do
   write(6,*) "Translation", inw, "for", ntot, "G vectors"
#ifdef __PARA
       end if
#endif
      end if
     end do

!   do inw=1,nw1
!     do i=1,ntot
!      write(36,222) indexplus(i,inw),  tagp(i,inw)
!      write(37,222) indexminus(i,inw),  tag(i,inw)
!     end do
!   end do
!222 format (1x,i10,1x,i10)

#ifdef __PARA

   call mpi_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_barrier 2' , ierr)

   call mpi_bcast(indexplus,nw1*ntot  , MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_bcast +' , ierr)


   call mpi_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_barrier 3' , ierr)

   call mpi_bcast(indexminus,nw1*ntot , MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_bcast -' , ierr)


   call mpi_barrier(MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_barrier 4' , ierr)

   call mpi_bcast(tag,nw1*ntot , MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
   if (ierr.ne.0) call errore('wfunc_init','mpi_bcast -' , ierr)

        call mpi_barrier(MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) call errore('wfunc_init','mpi_barrier 4' , ierr)

        call mpi_bcast(tagp,nw1*ntot , MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) call errore('wfunc_init','mpi_bcast -' , ierr)

   if (me.eq.1) then
#endif
      deallocate(bigg)
      deallocate(bign)
#ifdef __PARA
   end if
#endif
      deallocate(i_1,j_1,k_1)
   
   return
   end subroutine wfunc_init   
!--------------------------------------------------------------
   subroutine grid_map
!--------------------------------------------------------------
   use efcalc
   use cell_base
!   use parms
   use smooth_grid_dimensions, nnrs=>nnrsx
   use parallel_include
#ifdef __PARA
   use para_mod
#endif
   implicit none

   integer ir1, ir2, ir3, ibig3
   allocate(xdist(nnrs))
   allocate(ydist(nnrs))
   allocate(zdist(nnrs))

   do ir3=1,nr3s
#ifdef __PARA
   ibig3 = ir3 - dffts%ipp( me )
   if(ibig3.gt.0.and.ibig3.le.dffts%npp(me)) then
#else
   ibig3=ir3
#endif
      do ir2=1,nr2s
       do ir1=1,nr1s
         xdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
      &                  ((ir1-1)/dfloat(nr1sx))
         ydist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                   &
      &                  ((ir2-1)/dfloat(nr2sx))
         zdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
      &                  ((ir3-1)/dfloat(nr3sx))
!         
      end do
      end do
#ifdef __PARA
   end if
#endif
   end do
   return
   end subroutine grid_map

!--------------------------------------------------------------
subroutine tric_wts(rp1,rp2,rp3,alat,wts)
!--------------------------------------------------------------

!***** This subroutine computes the weights to be used for ****
!***** R.P. translations in the WF calculation in the case ****
!***** of ibrav=0 or ibrav=14 *********************************

   implicit none
   real(kind=8),intent(in) :: rp1(3), rp2(3), rp3(3)
   real(kind=8),intent(in) :: alat
   real(kind=8),intent(out) :: wts(6) 
   real(kind=8) :: pi, tpiba, tpiba2
   real(kind=8) :: b1x, b2x, b3x, b1y, b2y, b3y, b1z, b2z, b3z
        integer :: i

!        WRITE(6,*) 'ENTERED TRIC_WTS !'
!        WRITE(6,*) 'b3(3) :',rp3(3)
!        WRITE(6,*) 'alat  :',alat
   pi = 2.0d0*asin(1.0)  
   tpiba = 2.0d0*pi/alat
   tpiba2 = tpiba**2

   b1x = rp1(1)*tpiba
   b2x = rp2(1)*tpiba
   b3x = rp3(1)*tpiba
   b1y = rp1(2)*tpiba
   b2y = rp2(2)*tpiba
   b3y = rp3(2)*tpiba
   b1z = rp1(3)*tpiba
   b2z = rp2(3)*tpiba
   b3z = rp3(3)*tpiba
!        WRITE(6,*) 'COMPUTING WEIGHTS NOW ...'

   wts(1) = tpiba2*(-b1z*b2x*b2z*b3x + b2y**2*b3x**2 + b1z*b2z*b3x**2 + & 
              b2z**2*b3x**2 - b1z*b2y*b2z*b3y - 2.d0*b2x*b2y*b3x*b3y + & 
              b2x**2*b3y**2 + b1z*b2z*b3y**2 + b2z**2*b3y**2 + & 
              b1z*b2x**2*b3z + b1z*b2y**2*b3z - b1z*b2x*b3x*b3z - & 
              2.d0*b2x*b2z*b3x*b3z - b1z*b2y*b3y*b3z - &
              2.d0*b2y*b2z*b3y*b3z + b2x**2*b3z**2 + b2y**2*b3z**2 + &
              b1x*(b2y**2*b3x + b2z**2*b3x - b2y*(b2x + b3x)*b3y - &
              b2z*(b2x + b3x)*b3z +  b2x*(b3y**2 + b3z**2)) + &
              b1y*(b2x**2*b3y - b2x*b3x*(b2y + b3y) + &
              b2z*b3y*(b2z - b3z) + b2y*(b3x**2 - b2z*b3z +  b3z**2)))/ &
              ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y &
               + b1y*b2x*b3z - b1x*b2y*b3z)**2)

   wts(2) = tpiba2*(b1z**2*(b2x*b3x + b3x**2 + b3y*(b2y + b3y)) + &
              b1y**2*(b2x*b3x + b3x**2 + b3z*(b2z + b3z)) - &
              b1z*(-b2z*(b3x**2 + b3y**2) + (b2x*b3x + b2y*b3y)*b3z + &
              b1x*(b2z*b3x + (b2x + 2.d0* b3x)*b3z)) - &
              b1y*(b1x*(b2y*b3x + (b2x + 2.d0*b3x)*b3y) + &
              b3y*(b1z*b2z + b2x*b3x + 2.d0*b1z*b3z + b2z*b3z) - &
              b2y*(b3x**2 - b1z* b3z + b3z**2)) + &
              b1x*(-b2y*b3x*b3y + b2x*b3y**2 - b2z*b3x*b3z + b2x*b3z**2 + &
              b1x*(b2y*b3y + b3y**2 + b3z*(b2z + b3z))))/ &                        
              ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
              b1y*b2x*b3z - b1x*b2y*b3z)**2)

   wts(3) = tpiba2*(b1z**2*(b2x**2 + b2x*b3x + b2y*(b2y + b3y)) - &
              b1y*(2.d0*b1z*b2y*b2z + b2x*b2y*b3x - b2x**2*b3y + &
              b1z*b2z*b3y - b2z**2*b3y + b1x*(2.d0*b2x*b2y + b2y*b3x + b2x*b3y) + &
              b1z*b2y*b3z + b2y*b2z*b3z) + b1y**2*(b2x**2 + b2x*b3x + b2z*(b2z + b3z)) - &
              b1z*(b2x*b2z*b3x + b2y*b2z*b3y - b2x**2*b3z - b2y**2*b3z + &
              b1x*(2.d0*b2x*b2z + b2z*b3x + b2x*b3z)) + &
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

!        WRITE(6,*) 'WEIGHTS ARE :'
!        WRITE(6,"(F10.6)") (wts(i),i=1,6)

   end subroutine tric_wts
!-----------------------------------------------------------------
        subroutine small_box_wf(i_1,j_1,k_1,nw1)
!-----------------------------------------------------------------
        use efcalc
        use cell_base
!        use parms
        use smooth_grid_dimensions, nnrs=>nnrsx
        use constants
        use wfparm
        use grid_dimensions, only : nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr => nnrx
        use parallel_include
#ifdef __PARA
        use para_mod
#endif
        implicit none

        integer ir1, ir2, ir3, ibig3 , inw
        real(kind=8) x
        integer , intent(in) :: nw1, i_1(nw1), j_1(nw1), k_1(nw1)
        allocate(expo(nnr,nw1))

        do inw=1,nw1

        write(6,*) inw ,":", i_1(inw), j_1(inw), k_1(inw)

        do ir3=1,nr3
#ifdef __PARA
        ibig3 = ir3 - dfftp%ipp( me )
        if(ibig3.gt.0.and.ibig3.le.dfftp%npp(me)) then
#else
        ibig3=ir3
#endif
           do ir2=1,nr2
                do ir1=1,nr1
                   x =  (((ir1-1)/dfloat(nr1x))*i_1(inw) +                          &
      &                  ((ir2-1)/dfloat(nr2x))*j_1(inw) +             &
      &                  ((ir3-1)/dfloat(nr3x))*k_1(inw))*0.5d0*fpi
                   expo(ir1+(ir2-1)*nr1x+(ibig3-1)*nr1x*nr2x,inw) =  cmplx(cos(x), -sin(x))
                end do
           end do
#ifdef __PARA
        end if
#endif
        end do
        end do
        return
        end subroutine small_box_wf
!-----------------------------------------------------------------------
      complex(kind=8) function boxdotgridcplx(irb,qv,vr)
!-----------------------------------------------------------------------
!
! Calculate \sum_i qv(r_i)*vr(r_i)  with r_i on box grid
! array qv(r) is defined on box grid, array vr(r)on dense grid
! irb   : position of the box in the dense grid
! Parallel execution: remember to sum the contributions from other nodes
!
!      use ion_parameters
      use grid_dimensions, nnr => nnrx
      use smallbox_grid_dimensions ,nnrb => nnrbx
      use cell_base
      use small_box
#ifdef __PARA
      use para_mod
#endif
      implicit none
      integer, intent(in):: irb(3)
      complex(kind=8), intent(in):: qv(nnrb), vr(nnr)
!
      integer ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
!
!      if(nfft.le.0.or.nfft.gt.2) call errore('box2grid','wrong data',nfft)

      boxdotgridcplx=(0.d0, 0.d0)

      do ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+mod(ibig3-1,nr3)
#ifdef __PARA
         ibig3 = ibig3 - dfftp%ipp( me )
         if (ibig3.gt.0.and.ibig3.le.dfftp%npp(me)) then
#endif
            do ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+mod(ibig2-1,nr2)
               do ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+mod(ibig1-1,nr1)
                  ibig=ibig1 + (ibig2-1)*nr1x + (ibig3-1)*nr1x*nr2x
                  ir  =ir1 + (ir2-1)*nr1bx + (ir3-1)*nr1bx*nr2bx
                  boxdotgridcplx = boxdotgridcplx + qv(ir)*vr(ibig)
               end do
            end do
#ifdef __PARA
         endif
#endif
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine write_rho_g(rhog)
!-----------------------------------------------------------------------

      use gvecp, only: ng=>ngm
      use gvec, only: gx, mill_l
      use elct
      use parallel_include
#ifdef __PARA 
      use para_mod
#endif

      implicit none

      real(kind=8), allocatable:: gnx(:,:), bigg(:,:)
      complex(kind=8) ,intent(in) :: rhog(ng,nspin) 
      complex(kind=8),allocatable :: bigrho(:)
      complex(kind=8) :: rhotmp_g(ng)
      integer ntot, i, j
#ifdef __PARA
      integer proc, ierr, root, ngdens(nproc),recvcount(nproc), displs(nproc)
      integer recvcount2(nproc), displs2(nproc)
#endif
      character (len=6) :: name
      character (len=15) :: name2

      allocate(gnx(3,ng))


#ifdef __PARA
      root=0
#endif

    do i=1,ng
      gnx(1,i)=gx(i,1)
      gnx(2,i)=gx(i,2)
      gnx(3,i)=gx(i,3)
    end do

#ifdef __PARA
    ntot=0
    do i=1,nproc
        ngdens(i)=(dfftp%ngl(i)+1)/2
    end do

    do proc=1,nproc
       recvcount(proc)=ngdens(proc)*3
       recvcount2(proc)=ngdens(proc)
       if(proc.eq.1) then
           displs(proc)=0
           displs2(proc)=0
       else
           displs(proc)=displs(proc-1)+recvcount(proc-1)
           displs2(proc)=displs2(proc-1)+recvcount2(proc-1)
       end if
       ntot=ntot+recvcount(proc)/3
    end do

    if(me.eq.1) then
       allocate(bigg(3,ntot))
    end if

    call mpi_barrier(MPI_COMM_WORLD,ierr)
    if(ierr.ne.0) call errore('write_rho_g0','mpi_barrier', ierr)
    call mpi_gatherv(gnx,recvcount(me),MPI_REAL8,                        &
                     bigg,recvcount,displs,MPI_REAL8,                    &
                     root,MPI_COMM_WORLD, ierr)
    if(ierr.ne.0) call errore('write_rho_g0','mpi_gatherv', ierr)
    do i=1,nspin

      rhotmp_g(1:ng)=rhog(1:ng,i)

      if(me.eq.1) then 
        allocate (bigrho(ntot))
      end if

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(ierr.ne.0) call errore('write_rho_g1','mpi_barrier', ierr)
      call mpi_gatherv(rhotmp_g,recvcount2(me),MPI_DOUBLE_COMPLEX,                &
                     bigrho,recvcount2,displs2,MPI_DOUBLE_COMPLEX,              &
                     root,MPI_COMM_WORLD, ierr)
      if(ierr.ne.0) call errore('write_rho_g1','mpi_gatherv', ierr)

    
      if(me.eq.1) then
      if(i.eq.1) name2="CH_DEN_G_PARA.1"
      if(i.eq.2) name2="CH_DEN_G_PARA.2"
       OPEN(unit=57, file=name2) 
       do j=1,ntot
          write(57,*) bigrho(j)
       end do
      close(57)
      deallocate(bigrho)
      end if

    write(6,*) "Charge density written to ", name2

    end do

    if(me.eq.1) then
      name="G_PARA"
       OPEN(unit=56, file=name) 
       do i=1,ntot
          write(56,*) bigg(:,i)
       end do
      close(56)
      deallocate(bigg)
    end if
    write(6,*) "G-vectors written to G_PARA"
#else
    ntot=ng
    allocate(bigg(3,ntot))
    bigg(1:3,1:ntot)=gnx(1:3,1:ng)
    do i=1,nspin
      allocate(bigrho(ntot))
      bigrho(1:ng)=rhog(1:ng,i)

      if(i.eq.1) name2="CH_DEN_G_SERL.1"
      if(i.eq.2) name2="CH_DEN_G_SERL.2"

       OPEN(unit=57, file=name2) 
       do j=1,ntot
          write(57,*) bigrho(j)
       end do
      close(57)
      deallocate(bigrho)

      write(6,*) "Charge density written to", name2

    end do

      name="G_SERL"
       OPEN(unit=56, file=name) 
       do i=1,ntot
          write(56,*) bigg(:,i)
       end do
      close(56)
   deallocate(bigg)
   write(6,*) "G-vectors written to G_SERL"
#endif
     deallocate(gnx)
     return
     end subroutine write_rho_g

!-----------------------------------------------------------------------
      subroutine macroscopic_average(rhog,tau0,e_tuned)
!-----------------------------------------------------------------------
      use gvec
      use elct
      use tune
      use cell_base
!      use ion_parameters
      use ions_base
      use constants
      use parallel_include
#ifdef __PARA
      use para_mod
#endif
      implicit none

      real(kind=8), allocatable:: gnx(:,:), bigg(:,:)
      complex(kind=8) ,intent(in) :: rhog(ng,nspin)
      complex(kind=8),allocatable :: bigrho(:)
      complex(kind=8) :: rhotmp_g(ng)
      complex(kind=8),parameter :: zero=(0.d0,0.d0), ci=(0.d0,1.d0)
      integer ntot, i, j, ngz, l, isa
      integer ,allocatable :: g_red(:,:)
#ifdef __PARA
      integer proc, ierr, root, ngdens(nproc),recvcount(nproc), displs(nproc)
      integer recvcount2(nproc), displs2(nproc)
#endif
       real(kind=8) zlen,vtot, pos(3,nax,nsx), a_direct(3,3),a_trans(3,3), e_slp, e_int
       real(kind=8), intent(out) :: e_tuned(3)
       real(kind=8), intent(in) :: tau0(3,nax)
       real(kind=8),allocatable :: v_mr(:), dz(:), gz(:), g_1(:,:), vbar(:), cd(:), v_final(:)
       real(kind=8), allocatable:: cdion(:), cdel(:), v_line(:), dist(:)
       complex(kind=8),allocatable :: rho_ion(:),v_1(:),vmac(:),rho_tot(:),rhogz(:), bigrhog(:)

       allocate(gnx(3,ng))


#ifdef __PARA
      root=0
#endif

    do i=1,ng
      gnx(1,i)=gx(i,1)
      gnx(2,i)=gx(i,2)
      gnx(3,i)=gx(i,3)
    end do

#ifdef __PARA
    ntot=0
    do i=1,nproc
        ngdens(i)=(dfftp%ngl(i)+1)/2
    end do

    do proc=1,nproc
       recvcount(proc)=ngdens(proc)*3
       recvcount2(proc)=ngdens(proc)
       if(proc.eq.1) then
           displs(proc)=0
           displs2(proc)=0
       else
           displs(proc)=displs(proc-1)+recvcount(proc-1)
           displs2(proc)=displs2(proc-1)+recvcount2(proc-1)
       end if
       ntot=ntot+recvcount(proc)/3
    end do

!    if(me.eq.1) then
       allocate(bigg(3,ntot))
       allocate(g_1(3,2*ntot-1))
       allocate(g_red(3,2*ntot-1))
!    end if

    call mpi_barrier(MPI_COMM_WORLD,ierr)
    if(ierr.ne.0) call errore('macroscopic_average','mpi_barrier', ierr)
    call mpi_gatherv(gnx,recvcount(me),MPI_REAL8,                        &
                     bigg,recvcount,displs,MPI_REAL8,                    &
                     root,MPI_COMM_WORLD, ierr)
    if(ierr.ne.0) call errore('macroscopic_avergae','mpi_gatherv', ierr)

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) call errore('macroscopic_average','mpi_barrier 2' , ierr)

    call mpi_bcast(bigg,3*ntot  , MPI_REAL8, root, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) call errore('macroscopic_average','mpi_bcast_bigg' , ierr)

      rhotmp_g(1:ng)=rhog(1:ng,1)

!      if(me.eq.1) then
        allocate (bigrho(ntot))
        allocate (bigrhog(2*ntot-1))
!      end if

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(ierr.ne.0) call errore('macroscopic_average','mpi_barrier', ierr)
      call mpi_gatherv(rhotmp_g,recvcount2(me),MPI_DOUBLE_COMPLEX,       &
                     bigrho,recvcount2,displs2,MPI_DOUBLE_COMPLEX,       &
                     root,MPI_COMM_WORLD, ierr)
      if(ierr.ne.0) call errore('macroscopic_avergae','mpi_gatherv', ierr)

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('macroscopic_average','mpi_barrier 3' , ierr)

      call mpi_bcast(bigrho,ntot  , MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('macroscopic_average','mpi_bcast_bigrho' , ierr)



#else
    ntot=ng
    allocate(bigg(3,ntot))
    allocate(g_1(3,2*ntot-1))
    allocate(g_red(3,2*ntot-1))
    bigg(1:3,1:ntot)=gnx(1:3,1:ng)
    allocate(bigrho(ntot))
    allocate(bigrhog(2*ntot-1))
    bigrho(1:ng)=rhog(1:ng,1)
#endif
     
!       e_tuned=0.d0
!       return
              
       allocate(v_mr(npts))
       allocate(v_final(npts))
       allocate(dz(npts))
       allocate(vbar(npts))
       allocate(cd(npts))
       allocate(cdel(npts))
       allocate(cdion(npts))

       !-- needed for non-orthogonal cells

       a_direct(1,1:3)=a1(1:3)
       a_direct(2,1:3)=a2(1:3)
       a_direct(3,1:3)=a3(1:3)

       a_trans=TRANSPOSE(a_direct)
  
       !--- end 

!       e_tuned=0.d0
!       return

       !--- Construct rho(-g) from rho(g). rgo(-g)=rho*(g)
 
          bigrhog(1:ntot)=bigrho(1:ntot)
          g_1(:,1:ntot)=bigg(:,1:ntot)
       do i=2,ntot
          bigrhog(ntot+i-1)=conjg(bigrho(i))
          g_1(:,ntot+i-1)=-bigg(:,i)
       end do

       !--- needed fot non-orthogonal cells
       
       do i=1,2*ntot-1
         g_red(:,i)=NINT(MATMUL(a_trans(:,:),g_1(:,i))*tpiba/(2.d0*pi))
       end do

       !--- end
!       e_tuned=0.d0
!       return
      !--- define the direction of the line

        xdir=1
        ydir=2

        if ((zdir).eq.1) xdir=3
        if ((zdir).eq.2) ydir=3

        if(zdir.eq.1) zlen=DSQRT(a1(1)**2+a1(2)**2+a1(3)**2)
        if(zdir.eq.2) zlen=DSQRT(a2(1)**2+a2(2)**2+a2(3)**2)
        if(zdir.eq.3) zlen=DSQRT(a3(1)**2+a3(2)**2+a3(3)**2)


       !--- We need the potentiail only along zdir, so pick the appropriate G-vectors with Gxdir=Gydir=0

       ngz=0
       do i=1,2*ntot-1
         if((g_red(xdir,i).eq.0).and.(g_red(ydir,i).eq.0)) ngz=ngz+1
       end do

!       write(6,*) 'ngz=',ngz

       allocate(gz(ngz))
       allocate(rhogz(ngz))
       allocate(rho_ion(ngz))
       allocate(rho_tot(ngz))
       allocate(vmac(ngz))
       allocate(v_1(ngz))
       
!       write(6,*) 'allocated'

       !--- The G-vectors are output in units of 2*pi/a, so convert them to the correct values

       j=0
       do i=1,2*ntot-1
         if((g_red(xdir,i).eq.0).and.(g_red(ydir,i).eq.0)) then
           j=j+1
           gz(j)=g_1(zdir,i)*tpiba
           rhogz(j)=bigrhog(i)
         end if
       end do

!       write(6,*) 'mapped'

!       e_tuned=0.d0
!       return

       isa = 0
       do i=1,nsp
         do j=1,na(i)
            isa = isa + 1
            pos(:,j,i)=tau0(:,isa)
         end do
       end do

       !--- Construct the ionic Charge density in G-space

       rho_ion=zero
       do j=1,ngz
         do i=1,nsp
           do l=1,na(i)
              rho_ion(j)=rho_ion(j)+zv(i)*exp(-ci*gz(j)*pos(zdir,l,i))*exp(-gz(j)**2/(4.d0*alpha))
           end do
         end do
       end do

!       write(6,*) 'rho_ion'

       rho_ion=rho_ion/omega

       !--- Construct the total Charge density in G-space

       rho_tot=rho_ion-rhogz

       !--- Construct the electrostatic potential and macroscopic average in G-space

       v_1(1)=zero
       vmac(1)=zero
       v_1(2:ngz)=4*pi*rho_tot(2:ngz)/gz(2:ngz)**2
       vmac(2:)=v_1(2:)*sin(gz(2:)*b)/(gz(2:)*b)


       !--- Calculate planewise average in R-space and FFT V(Gz) ---> V(z) ... well not really FFT but FT

       vbar=0.d0
       v_mr=0.d0
       cdel=0.d0
       cdion=0.d0
       cd=0.d0
       do j=1,npts
         dz(j)=(j-1)*zlen/(npts*1.d0)
         do i=1,ngz
           vbar(j)=vbar(j)-Real(exp(ci*gz(i)*dz(j))*v_1(i))
           v_mr(j)=v_mr(j)-Real(exp(ci*gz(i)*dz(j))*vmac(i))
           cdel(j)=cdel(j)-Real(exp(ci*gz(i)*dz(j))*rhogz(i))
           cdion(j)=cdion(j)+Real(exp(ci*gz(i)*dz(j))*rho_ion(i))
           cd(j)=cd(j)+Real(exp(ci*gz(i)*dz(j))*rho_tot(i))
         end do
!           write(6,*) vbar(j), v_mr(j), cdel(j), cdion(j)
       end do
       if (shift) then
           vtot=(v_mr(start)+v_mr(start-1))/2.d0
           v_final(1:npts-start+1)=v_mr(start:npts)-vtot
           v_final(npts-start+2:npts)=v_mr(1:start-1)-vtot
       else
           vtot=(v_mr(1)+v_mr(npts))/2.d0
           v_final(1:npts)=v_mr(1:npts)-vtot
       end if
!       write(6,*) vtot
       e_tuned=0.d0

       allocate(v_line(1:av1-av0+1))
       allocate(dist(1:av1-av0+1))


       v_line(1:av1-av0+1)=v_final(av0:av1)
       dist(1:av1-av0+1) =dz(av0:av1)

!       call least_square(av1-av0+1, dist, v_line, e_slp, e_int)

!       e_tuned(zdir)=-e_slp

       e_tuned(zdir)=-(v_final(av1)-v_final(av0))/((av1-av0)*zlen/(npts*1.d0))
!       write(6,*) 'e_tuned=',e_tuned
       
#ifdef __PARA
!       if(me.eq.1) then
        deallocate(bigg,g_1,bigrho,bigrhog,g_red)     
!       end if
#else
        deallocate(bigg,g_1,bigrho,bigrhog,g_red)     
#endif

     deallocate(gnx,v_mr,v_final,dz,vbar,cd,cdel,cdion)
     deallocate(v_line, dist)
     deallocate(gz,rhogz,rho_ion,rho_tot,vmac,v_1)

!     write(6,*) 'deallocated'
 
     return
     end subroutine macroscopic_average
!--------------------------------------------------------------------
     subroutine least_square(npts,x,y,slope,intercept)
!--------------------------------------------------------------------

     implicit none
     integer,intent(in):: npts
     integer i
     real(kind=8),intent(in) :: x(npts),y(npts)
     real(kind=8), intent(out):: slope, intercept 
     real(kind=8) :: sumx,sumy,sumx2,sumxy,sumsqx
     real(kind=8) :: xav,yav

     sumxy=0.d0
     sumx =0.d0
     sumy =0.d0
     sumx2=0.d0
    do i=1,npts
     sumxy=sumxy+x(i)*y(i)
     sumx =sumx +x(i)
     sumy =sumy +y(i)
     sumx2=sumx2+x(i)*x(i)
    end do
     sumsqx=sumx**2
     xav=sumx/dfloat(npts)
     yav=sumy/dfloat(npts)

     slope=(npts*sumxy - sumx*sumy)/(npts*sumx2 - sumsqx)

     intercept=yav-slope*xav

     return

     end subroutine least_square
!-----------------------------------------------------------------------
      subroutine wfsteep(m, Omat, Umat,b1,b2,b3)
!-----------------------------------------------------------------------
    use wfparm
    use wfparm2
    use control_flags
    use cell_base
    use parallel_include
#ifdef __PARA
      use para_mod
#endif
  implicit none

!    (m,m) is the size of the matrix Ospin.
!    Ospin is input overlap matrix.
!    Uspin is the output unitary transformation.
!             Rough guess for Uspin can be carried in.
!
!     conjugated gradient to search maximization
!
  real(kind=8), parameter :: autoaf=0.529177d0
  integer, intent(in) :: m
  real(kind=8), intent(in) :: b1(3),b2(3),b3(3)
  complex(kind=8), intent(inout) :: Omat(nw, m, m)
  real(kind=8), intent(inout) :: Umat(m,m)
!
  integer :: i, j, k, l, ig, ierr, ti, tj, tk, inw, ir, adjust
!  integer :: f3(nw), f4(nw),isteep , ierr1
  integer :: f3(nw), f4(nw), ierr1
  real(kind=8) :: slope, slope2, t1, t2, t3, pi2, mt(nw),t21,temp1,maxdt
  real(kind=8) :: U(m,m), wfc(3, m), Wm(m,m), schd(m,m), f2(4*m), gr(nw, 3)
  real(kind=8) :: Uspin2(m,m),temp2,wfdtold,oldt1,t01, d3(m,m), d4(m,m), U1(m,m)
  real(kind=8) :: spread, sp
  real(kind=8), allocatable, dimension(:) :: wr
  real(kind=8), allocatable, dimension(:,:) :: W
  complex(kind=8) :: ci, ct1, ct2, ct3, z(m, m), X(m, m), d(m,m), d2(m,m)
  complex(kind=8) :: f1(2*m-1), wp(m*(m+1)/2), Oc(nw, m, m), alpha, beta1
  complex(kind=8) ::  Oc2(nw, m, m),wp1(m*(m+1)/2), X1(m,m), U2(m,m), U3(m,m)

!
  ci=(0.d0,1.d0)
  alpha=(1.0d0, 0.0d0)
  beta1=(0.0d0, 0.0d0)
  pi2=2.d0*3.14159265358979d0
!
  allocate(W(m,m), wr(m))


  Umat=0.d0
  do i=1,m
    Umat(i,i)=1.d0
  end do
  Oc=beta1
  Oc2=beta1
  X1=beta1
  U2=Umat*alpha

!
! update Oc using the initial guess of Uspin
!
  do inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=beta1
     call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
     X1=beta1
     call ZGEMM ('N','N', m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
     Oc(inw, :, :)=X1(:, :)
  end do

     U2=beta1
     U3=beta1

  W=0.d0
  schd=0.d0
  oldt1=0.d0
  wfdtold=0.d0

  do k=1, nit
    t01=0.d0     !use t1 to store the value of omiga
    do inw=1, nw
       do i=1, m
          t01=t01+real(conjg(Oc(inw, i, i))*Oc(inw, i, i))
       end do
    end do

!    write(6,*) t01

    if(ABS(oldt1-t01).lt.tolw) then 
#ifdef __PARA
        if(me.eq.1) then
#endif
           write(27,*) "MLWF Generated at Step",k
#ifdef __PARA
        end if
#endif
        if(iprsta.gt.4) then
           write(6,*) "MLWF Generated at Step",k
        end if
       go to 40
      end if
    
!    oldt1=t01
  
!   calculate d(omiga)/dW and store result in W
!   W should be a real symmetric matrix for gamma-point calculation
!
    Wm=W
    W=0.d0
    do inw=1, nw
       t2=weight(inw)
       do i=1,m
          do j=i+1,m
             W(i,j)=W(i,j)+t2*real(Oc(inw,i,j)*conjg(Oc(inw,i,i)        &
              -Oc(inw,j,j))+conjg(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
          end do
       end do
    end do
    W=W-transpose(W)
  
!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(W**2)
  
!   calculate slope2=d2(omiga)/d(lamda)2
    slope2=0.d0
    do ti=1, m
       do tj=1, m
          do tk=1, m
             t2=0.d0
             do inw=1, nw
                t2=t2+real(Oc(inw,tj,tk)*conjg(Oc(inw,tj,tj)+Oc(inw,tk,tk) &
                          -2.d0*Oc(inw,ti,ti))-4.d0*Oc(inw,ti,tk)          &
                          *conjg(Oc(inw,ti,tj)))*weight(inw)
             end do
             slope2=slope2+W(tk,ti)*W(ti,tj)*2.d0*t2
          end do
       end do
     end do
    slope2=2.d0*slope2
  
!   use parabola approximation. Defined by 1 point and 2 slopes
    if (slope2.lt.0) wfdt=-slope/2.d0/slope2
    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
  
    if (k.lt.nsd) then
       schd=W    !use steepest-descent technique

!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(schd**2)

!       schd=schd*maxwfdt
    do i=1, m
       do j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0, schd(i,j))
       end do
    end do

#if ! defined __AIX
    call zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
#else
    call zhpev(21, wp1, wr, z, m, m, f2, 4*m)
    ierr1 = 0
#endif

    if (ierr.ne.0) stop 'failed to diagonalize W!'

    else
!
!     construct conjugated gradient
!        d(i)=-g(i)+belta(i)*d(i-1)
!        belta^FR(i)=g(i)t*g(i)/(g(i-1)t*g(i-1))
!        belta^PR(i)=g(i)t*(g(i)-g(i-1))/(g(i-1)t*g(i-1))
!
        call DGEMM ('T','N', m,m,m,alpha,Wm,m,Wm,m,beta1,d3,m)

       
       t1=0.d0
       do i=1, m
          t1=t1+d3(i, i)
       end do
       if (t1.ne.0) then
          d4=(W-Wm)
          call DGEMM ('T','N', m,m,m,alpha,W,m,d4,m,beta1,d3,m)
          t2=0.d0
          do i=1, m
             t2=t2+d3(i, i)
          end do
          t3=t2/t1
          schd=W+schd*t3
       else
          schd=W
       end if
!
!   calculate the new d(Lambda) for the new Search Direction
!   added by Manu. September 19, 2001
!
!   calculate slope=d(omiga)/d(lamda)
    slope=SUM(schd**2)
!------------------------------------------------------------------------
!   schd=schd*maxwfdt
    do i=1, m
       do j=i, m
        wp1(i + (j-1)*j/2) = cmplx(0.0, schd(i,j))
       end do
    end do

#if ! defined __AIX
    call zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
#else
    call zhpev(21, wp1, wr, z, m, m, f2, 4*m)
    ierr1 = 0
#endif
    if (ierr.ne.0) stop 'failed to diagonalize W!'

      maxdt=maxwfdt

11    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*(maxwfdt)*wr(i))
    end do      !d=exp(d)
 
!   U=z*exp(d)*z+
     U3=beta1
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=real(U2)
     U2=beta1
     U3=beta1
!
!   update Uspin
    U1=beta1
    call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1
 
!    Uspin2=matmul(Uspin, U2)
!
!   update Oc
!
     U2=Umat*alpha
     U3=beta1
     do inw=1, nw
      X1(:,:)=Omat(inw,:,:)
      call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
      X1=beta1
      call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
      Oc2(inw, :,:)=X(:,:)
     end do
     U2=beta1 
     U3=beta1
!
    t21=0.d0     !use t21 to store the value of omiga
    do inw=1, nw
       do i=1, m
          t21=t21+real(conjg(Oc2(inw, i, i))*Oc2(inw, i, i))
       end do
    end do
 
      temp1=-((t01-t21)+slope*maxwfdt)/(maxwfdt**2)
      temp2=slope
      wfdt=-temp2/(2*temp1)

        if (wfdt.gt.maxwfdt.or.wfdt.lt.0.d0) then
        maxwfdt=2*maxwfdt
        go to 11
        end if

        maxwfdt=maxdt
!
!
!   use parabola approximation. Defined by 2 point and 1 slopes
!    if (slope2.lt.0) wfdt=-slope/2.d0/slope2
!    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
!
!    write(6, '(e12.5E2,1x,e11.5E2,1x,f6.2)') slope2, slope, wfdt
!-------------------------------------------------------------------------
!
!      schd is the new searching direction
!
    end if
  
    d=0.d0
    do i=1, m
       d(i, i)=exp(ci*wfdt*wr(i))
    end do          !d=exp(d)

 
!   U=z*exp(d)*z+
!
     U3=beta1
     call ZGEMM ('N', 'N', m,m,m,alpha,z,m,d,m,beta1,U3,m)
     U2=beta1
     call ZGEMM ('N','c', m,m,m,alpha,U3,m,z,m,beta1,U2,m)
     U=real(U2)
     U2=beta1
     U3=beta1

!   update Uspin
!
    U1=beta1
    call DGEMM ('N', 'N', m,m,m,alpha,Umat,m,U,m,beta1,U1,m)
    Umat=U1

!   update Oc
!
       U2=Umat*alpha
       U3=beta1
     do inw=1, nw
       X1(:, :)=Omat(inw, :, :)
       call ZGEMM ('T', 'N', m,m,m,alpha,U2,m,X1,m,beta1,U3,m)
       X1=beta1
       call ZGEMM ('N','N',m,m,m,alpha,U3,m,U2,m,beta1,X1,m)
       Oc(inw, :, :)=X1(:, :)
     end do
    U2=beta1
    U3=beta1
          if(ABS(t01-oldt1).ge.tolw.and.k.ge.nit) then
#ifdef __PARA
        if(me.eq.1) then
#endif
           write(27,*) "MLWF Not generated after",k,"Steps."
#ifdef __PARA
        end if
#endif
        if(iprsta.gt.4) then
           write(6,*) "MLWF Not generated after",k,"Steps."
        end if
           go to 40
        end if
    oldt1=t01
  end do

40  deallocate(W, wr)

!
! calculate the spread
!
!  write(24, *) "spread: (unit \AA^2)"
  do i=1, m
     mt=1.d0-real(Oc(:,i,i)*conjg(Oc(:,i,i)))
     sp = (alat*autoaf/pi2)**2*SUM(mt*weight)
#ifdef __PARA
       if(me.eq.1) then
#endif
       write(25, '(f10.7)') sp
#ifdef __PARA
       end if
#endif
      if(sp.lt.0.d0) then
      write(6,*) "Something wrong WF Spread negative. The Program will Stop."
#ifdef __PARA
      call MPI_FINALIZE(ierr1)
#endif
      STOP
     end if
     spread=spread+sp
  end do
  spread=spread/dfloat(m)

#ifdef __PARA
        if(me.eq.1) then
#endif
    write(24, '(f10.7)') spread
    write(27,*) "Average spread = ", spread
#ifdef __PARA
        end if
#endif

!
! calculate wannier-function centers
!
!allocate(wr(nw), W(nw, nw))
! do inw=1, nw
!     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
!  end do
!
! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
! to check the correctness of choices on G vectors
!
!  do i=1, nw
!     do j=1, nw
!        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
!        write(6, *) i,j,W(i,j)
!     end do
!  end do
!! write(24, *) "wannier function centers: (unit:\AA)"
!  do i=1, m
!     mt=-aimag(log(Oc(:,i,i)))/pi2
!     wfc(1, i)=SUM(mt*weight*gr(:,1))
!     wfc(2, i)=SUM(mt*weight*gr(:,2))
!     wfc(3, i)=SUM(mt*weight*gr(:,3))
!     do inw=1, nw
!        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
!     end do
!     mt=wr
!     wfc(1, i)=(wfc(1,i)+SUM(mt*weight*gr(:,1)))*alat
!     wfc(2, i)=(wfc(2,i)+SUM(mt*weight*gr(:,2)))*alat
!     wfc(3, i)=(wfc(3,i)+SUM(mt*weight*gr(:,3)))*alat
!  end do
!
!
!  do i=1, m
!     write(26, '(3f11.6)') wfc(:,i)*autoaf
!  end do
!  deallocate(wr, W)
!   deallocate (mt,X1,Oc)
        if(iprsta.gt.4) then
        write(6,*) "Leaving WFSTEEP"
        end if
  return
  end subroutine wfsteep
!
!-------------------------------------------------------------------------
      subroutine dforce_field (bec,deeq,betae,i,c,ca,df,da,v,v1)
!-----------------------------------------------------------------------
!computes: the generalized force df=cmplx(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=cmplx(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
      use control_flags, only: iprint, tbuff
      use gvec
      use gvecs
      use gvecw, only: ngw
      use cvan
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use elct
      use constants, only: pi, fpi
      use ions_base, only: nsp, na, nat
      use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
      use uspp_param, only: nh, nhm
      use uspp, only : nhsa=> nkb, dvan
!
      implicit none
!
      complex(kind=8) betae(ngw,nhsa), c(ngw), ca(ngw), df(ngw), da(ngw)
      real(kind=8) bec(nhsa,n), deeq(nhm,nhm,nat,nspin), v(nnrsx,nspin), v1(nnrsx,nspin)
      integer i
! local variables
      integer iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      real(kind=8) fi, fip, dd
      complex(kind=8) fp,fm,ci
      real(kind=8) af(nhsa), aa(nhsa) ! automatic arrays
      complex(kind=8)  dtemp(ngw)    !
      complex(kind=8), allocatable :: psi(:)
!
!
      call start_clock( 'dforce_field' )
      allocate( psi( nnrsx ) )
!
!     important: if n is odd => c(*,n+1)=0.
! 
      if (mod(n,2).ne.0.and.i.eq.n) then
         do ig=1,ngw
            ca(ig)=(0.,0.)
         end do
      endif
!
      ci=(0.0,1.0)
!
      if (.not.tbuff) then
!
         psi( 1:nnrsx ) = 0.0d0 ! call zero(2*nnrsx,psi)
         do ig=1,ngw
            psi(nms(ig))=conjg(c(ig)-ci*ca(ig))
            psi(nps(ig))=c(ig)+ci*ca(ig)
         end do
!
         call ivfftw(psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!     
      else
!
!     read psi from buffer 21
!
#if defined(__CRAYY)
         buffer in(21,0) (psi(1),psi(nnrsx))
         ios = unit(21)
#else
         read(21,iostat=ios) psi
#endif
         if(ios.ne.0) call errore                                        &
     &       (' dforce',' error in reading unit 21',ios)
!
      endif
! 
      iss1=ispin(i)
!
! the following avoids a potential out-of-bounds error
!
      if (i.ne.n) then
         iss2=ispin(i+1)
      else
         iss2=iss1
      end if
!
      do ir=1,nnrsx
         psi(ir)=cmplx(v(ir,iss1)* real(psi(ir)),                       &
     &                 v1(ir,iss2)*aimag(psi(ir)) )
      end do
!
      call fwfftw(psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
!     note : the factor 0.5 appears 
!       in the kinetic energy because it is defined as 0.5*g**2
!       in the potential part because of the logics
!
      fi =-  f(i)*0.5
      fip=-f(i+1)*0.5
      do ig=1,ngw
         fp= psi(nps(ig)) + psi(nms(ig))
         fm= psi(nps(ig)) - psi(nms(ig))
         df(ig)= fi*(tpiba2*ggp(ig)* c(ig)+cmplx(real(fp), aimag(fm)))
         da(ig)=fip*(tpiba2*ggp(ig)*ca(ig)+cmplx(aimag(fp),-real(fm)))
      end do
!
!     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
! 
      if(nhsa.gt.0)then
         do inl=1,nhsa
            af(inl)=0.
            aa(inl)=0.
         end do
!
         do is=1,nsp
            do iv=1,nh(is)
               do jv=1,nh(is)
                  isa=0
                  do ism=1,is-1
                     isa=isa+na(ism)
                  end do
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     isa=isa+1
                     dd = deeq(iv,jv,isa,iss1)+dvan(iv,jv,is)
                     af(inl)=af(inl)-  f(i)*dd*bec(jnl,  i)
                     dd = deeq(iv,jv,isa,iss2)+dvan(iv,jv,is)
                     if (i.ne.n) aa(inl)=aa(inl)-f(i+1)*dd*bec(jnl,i+1)
                  end do
               end do
            end do
         end do
!
         do ig=1,ngw
            dtemp(ig)=(0.,0.)
         end do
         call MXMA                                                      &
     &        (betae,1,2*ngw,af,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         do ig=1,ngw
            df(ig)=df(ig)+dtemp(ig)
         end do
!
         do ig=1,ngw
            dtemp(ig)=(0.,0.)
         end do
         call MXMA                                                      &
     &        (betae,1,2*ngw,aa,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         do ig=1,ngw
            da(ig)=da(ig)+dtemp(ig)
         end do
      endif

      deallocate( psi )
!
      call stop_clock( 'dforce_field' )
!
      return
      end
!----------------------------------------------------------------------
#ifdef __PARA
!-----------------------------------------------------------------------
      subroutine write_psi(c,jw)
!----------------------------------------------------------------------
! for calwf 5             - M.S
! collect wavefunctions on first node and write to file
!
      use gvec
      use gvecs
      use elct
      use para_mod
      use smooth_grid_dimensions , nnrs => nnrsx
      use gvecw , only : ngw
      use parallel_include
!
      implicit none

      integer unit, jw
      complex*16 c(ngw,nx)
      complex*16, allocatable :: psis(:)
!
      integer i, ii, ig, proc, ierr, ntot, ncol, mc,ngpwpp(nproc)
      integer nmin(3), nmax(3), n1,n2,nzx,nz,nz_
      integer root, displs(nproc), recvcount(nproc)
      complex*16, allocatable:: psitot(:), psiwr(:,:,:)
!
! nmin, nmax are the bounds on (i,j,k) indexes of wavefunction G-vectors
!
      call nrbounds(ngw,nr1s,nr2s,nr3s,mill_l,nmin,nmax)
!
! nzx is the maximum length of a column along z
!
      nzx=nmax(3)-nmin(3)+1
!
      root = 0
! root is the first node
      ntot = 0


        do proc=1,nproc
          ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
          ntot=ntot+ngpwpp(proc)
        end do

      do proc=1,nproc

         recvcount(proc) = ngpwpp(proc)
!
!!
! recvcount(proc) = size of data received from processor proc
!                   (number of columns times length of each column)
!
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
!
! displs(proc) is the position of data received from processor proc
!
!         ntot = ntot + recvcount(proc)
!
! ntot = total size of gathered data
!
      end do
!
! allocate the needed work spaces
!
      allocate( psis( nnrs ) ) 
      psis( 1:nnrs ) = 0.0d0 ! call zero(2*nnrs,psis)
      if (me.eq.1) then
         allocate(psitot(ntot))
         allocate(psiwr(nmin(3):nmax(3),nmin(1):nmax(1),nmin(2):nmax(2)))
!         write(unit) n, nmin, nmax
      end if
!
! fill array psis with c_i(G) (as packed columns along z)
!
!      do i=1,n
         do ig=1,ngw
!
! ncol+1 is the index of the column
!
            ncol=(nps(ig)-1)/nr3sx
!
! nz_ is the z component in FFT style (refolded between 1 and nr3s)
!
            nz_ =nps(ig)-ncol*nr3sx
!
! nz is the z component in "natural" style (between nmin(3) and nmax(3))
!
            nz  =nz_-1
            if (nz.ge.nr3s/2) nz=nz-nr3s

! ncpw(me) columns along z are stored in contiguous order on each node
!
!            psis(nz-nmin(3)+1+ncol*nzx)=c(ig,jw)
             psis(ig)=c(ig,jw)
         end do
!
! gather all psis arrays on the first node, in psitot
!
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_barrier 2',ierr)
         call mpi_gatherv (psis, recvcount(me),     MPI_DOUBLE_COMPLEX, &
    &                      psitot,recvcount, displs,MPI_DOUBLE_COMPLEX, &
    &                 root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('write_wfc','mpi_gatherv',ierr)
!
! write the node-number-independent array
!
        if(me.eq.1) then
        do i=1,ntot
            write(22,*) psitot(i)
        end do
            write(6,*) "State Written", jw
        end if
!
      if (me.eq.1) then
         deallocate(psiwr)
         deallocate(psitot)
      end if

      deallocate( psis ) 
!
      return
!
      end subroutine write_psi
!
#endif
!-----------------------------------------------------------------------
   subroutine rhoiofr (nfi,c,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin,ndwwf)
!-----------------------------------------------------------------------
!     the normalized electron density rhor in real space
!     the kinetic energy ekin
!     subroutine uses complex fft so it computes two ft's
!     simultaneously
!
!     rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
!     < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                   2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     e_v = sum_i,ij rho_i,ij d^ion_is,ji
!
      use control_flags, only: iprint, tbuff, iprsta, thdyn, tpre, trhor
      use ions_base , nas => nax
      use parameters, only: natx, nsx
      use gvec
      use gvecs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use cvan
      !use parm
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use cell_base, only: omega
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use elct
      use constants, only: pi, fpi
      use pseu
      use wfparm, only : iwf
!
      use cdvan
      use dener
      use io_global, only: stdout
      use uspp_param, only: nh, nhm
      use uspp, only : nhsa=> nkb
!
      implicit none
      real(kind=8) bec(nhsa,n), rhovan(nhm*(nhm+1)/2,nat,nspin)
      real(kind=8) rhovanaux(nhm,nhm,nat,nspin)
      real(kind=8) rhor(nnr,nspin), rhos(nnrsx,nspin)
      real(kind=8) enl, ekin
      complex(kind=8) eigrb(ngb,nas,nsp), c(ngw,nx), rhog(ng,nspin)
      integer irb(3,natx,nsx), nfi, ndwwf
! local variables
      integer iss, isup, isdw, iss1, iss2, ios, i, ir, ig
      integer is,iv,jv,isa,isn, jnl, j, k, inl, ism, ia
      real(kind=8) rsumr(2), rsumg(2), sa1, sa2, sums(2)
      real(kind=8) rnegsum, rmin, rmax, rsum
      real(kind=8) enkin, ennl
      complex(kind=8) ci,fp,fm
      complex(kind=8), allocatable :: psi(:), psis(:)
      external ennl, enkin
!
!
      call start_clock(' rhoiofr ')
      allocate( psi( nnr ) )
      allocate( psis( nnrsx ) ) 
      ci=(0.0,1.0)
      do iss=1,nspin
         rhor( 1:nnr, iss ) = 0.0d0 ! call zero(nnr,rhor(1,iss))
         rhos( 1:nnrsx, iss ) = 0.0d0 ! call zero(nnrsx,rhos(1,iss))
         rhog( 1:ng, iss ) = 0.0d0 ! call zero(2*ng,rhog(1,iss))
      end do
!
!     ==================================================================
!     calculation of kinetic energy ekin
!     ==================================================================
      ekin=enkin(c)
      if(tpre) call denkin(c,dekin)
!
!     ==================================================================
!     calculation of non-local energy
!     ==================================================================
!      enl=ennl(rhovan,bec)
       do is=1,nsp
           do iv=1, nh(is)
              do jv=iv,nh(is)
                 isa=0
                 do ism=1,is-1
                    isa=isa+na(ism)
                 end do
                 do ia=1,na(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    jnl=ish(is)+(jv-1)*na(is)+ia
                    isa=isa+1
                    sums(1)=f(iwf)*bec(inl,iwf)*bec(jnl,iwf)
                    rhovanaux(iv,jv,isa,1) = sums(1)
                    rhovanaux(jv,iv,isa,1) = sums(1)
                 end do
             end do
         end do
      end do

      k=1
      do i=1,nhm
        do j=i,nhm
           rhovan(k,:,:)=rhovanaux(j,i,:,:)
           k=k+1
        end do
      end do
   
      if(tpre) call dennl(bec,denl)
!    
!    warning! trhor and thdyn are not compatible yet!   
!
      if(trhor.and.(.not.thdyn))then
!     ==================================================================
!     charge density is read from unit 47
!     ==================================================================
#ifdef __PARA
         call read_rho(47,nspin,rhor)
#else
         read(47) ((rhor(ir,iss),ir=1,nnr),iss=1,nspin)
#endif
         rewind 47
!
         if(nspin.eq.1)then
            iss=1
            do ir=1,nnr
               psi(ir)=cmplx(rhor(ir,iss),0.)
            end do
            call fwfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ig=1,ng
               rhog(ig,iss)=psi(np(ig))
            end do
         else
            isup=1
            isdw=2
            do ir=1,nnr
               psi(ir)=cmplx(rhor(ir,isup),rhor(ir,isdw))
            end do
            call fwfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ig=1,ng
               fp=psi(np(ig))+psi(nm(ig))
               fm=psi(np(ig))-psi(nm(ig))
               rhog(ig,isup)=0.5*cmplx( real(fp),aimag(fm))
               rhog(ig,isdw)=0.5*cmplx(aimag(fp),-real(fm))
            end do
         endif
!
      else
!     ==================================================================
!     self-consistent charge
!     ==================================================================
!
!     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
! 
!         if (mod(n,2).ne.0) then
!            do ig=1,ngw
!               c(ig,n+1)=(0.,0.)
!            end do
!         endif
!
!         do i=1,n,2
            i=iwf
            psis( 1:nnrsx ) = 0.0d0 ! call zero(2*nnrsx,psis)
            do ig=1,ngw
!               c(ig,i+1)=(0.,0.)
               psis(nms(ig))=conjg(c(ig,i))
               psis(nps(ig))=c(ig,i)
            end do
!
            call ivfftw(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
!     wavefunctions in unit 21
!
#if defined(__CRAYY)
            if(tbuff) buffer out(21,0) (psis(1),psis(nnrsx))
#else
            if(tbuff) write(21,iostat=ios) psis
#endif
!            iss1=ispin(i)
            iss1=1
            sa1=f(i)/omega
!            if (i.ne.n) then
!              iss2=ispin(i+1)
!              sa2=f(i+1)/omega
!            else
               iss2=iss1  ! carlo
               sa2=0.0
!            end if
            do ir=1,nnrsx
               rhos(ir,iss1)=rhos(ir,iss1) + sa1*( real(psis(ir)))**2
               rhos(ir,iss2)=rhos(ir,iss2) + sa2*(aimag(psis(ir)))**2
            end do
!
!       buffer 21
!     
            if(tbuff) then
#if defined(__CRAYY)
               ios=unit(21)
#endif
               if(ios.ne.0) call errore                                  &
     &              (' rhoofr',' error in writing unit 21',ios)
            endif
!
!         end do
!
         if(tbuff) rewind 21
!
!     smooth charge in g-space is put into rhog(ig)
!
         if(nspin.eq.1)then
            iss=1
            do ir=1,nnrsx
               psis(ir)=cmplx(rhos(ir,iss),0.)
            end do
            call fwffts(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            do ig=1,ngs
               rhog(ig,iss)=psis(nps(ig))
            end do
         else
            isup=1
            isdw=2
             do ir=1,nnrsx
               psis(ir)=cmplx(rhos(ir,isup),rhos(ir,isdw))
            end do
            call fwffts(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            do ig=1,ngs
               fp= psis(nps(ig)) + psis(nms(ig))
               fm= psis(nps(ig)) - psis(nms(ig))
               rhog(ig,isup)=0.5*cmplx( real(fp),aimag(fm))
               rhog(ig,isdw)=0.5*cmplx(aimag(fp),-real(fm))
            end do
         endif
!
         if(nspin.eq.1) then
!     ==================================================================
!     case nspin=1
!     ------------------------------------------------------------------
            iss=1
            psi( 1:nnr ) = 0.0d0 ! call zero(2*nnr,psi)
            do ig=1,ngs
               psi(nm(ig))=conjg(rhog(ig,iss))
               psi(np(ig))=      rhog(ig,iss)
            end do
            call invfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ir=1,nnr
               rhor(ir,iss)=real(psi(ir))
            end do
         else 
!     ==================================================================
!     case nspin=2
!     ------------------------------------------------------------------
            isup=1
            isdw=2
            psi( 1:nnr ) = 0.0d0 ! call zero(2*nnr,psi)
            do ig=1,ngs
               psi(nm(ig))=conjg(rhog(ig,isup))+ci*conjg(rhog(ig,isdw))
               psi(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            end do
            call invfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ir=1,nnr
               rhor(ir,isup)= real(psi(ir))
               rhor(ir,isdw)=aimag(psi(ir))
            end do
         endif
!
!         if(iprsta.ge.3)then
          write( stdout,*) 'Smooth part of charge density :'
            do iss=1,nspin
               rsumg(iss)=omega*real(rhog(1,iss))
               rsumr(iss)=SUM(rhor(1:nnr,iss))*omega/dfloat(nr1*nr2*nr3)
            end do
#ifdef __PARA
            if (ng0.ne.2) then
! in the parallel case, only one processor has G=0 ! 
               do iss=1,nspin
                  rsumg(iss)=0.0
               end do
            end if
            call reduce(nspin,rsumg)
            call reduce(nspin,rsumr)
#endif
            if (nspin.eq.1) then
               WRITE( stdout,1) rsumg(1),rsumr(1)
            else
               WRITE( stdout,2) (rsumg(iss),iss=1,nspin),(rsumr(iss),iss=1,nspin)
            endif
!         endif
!     ==================================================================
!
!     add vanderbilt contribution to the charge density
!
!     drhov called before rhov because input rho must be the smooth part
!
         if (tpre) call drhov(irb,eigrb,rhovan,rhog,rhor)
!
         call rhov(irb,eigrb,rhovan,rhog,rhor)
      endif
!     ======================================endif for trhor=============
        rewind ndwwf
#ifdef __PARA
        call write_rho(ndwwf,nspin,rhor)
#else
        write(ndwwf,'(f12.7)')     &
          ((rhor(ir,iss),ir=1,nnr),iss=1,nspin)
#endif


!
!     here to check the integral of the charge density
!
!
      if(iprsta.ge.2) then
         call checkrho(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
         rnegsum=rnegsum*omega/dfloat(nr1*nr2*nr3)
         rsum=rsum*omega/dfloat(nr1*nr2*nr3)
         WRITE( stdout,'(a,4(1x,f12.6))')                                     &
     &     ' rhoofr: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
      end if
!
      if(nfi.eq.0.or.mod(nfi-1,iprint).eq.0) then

         write( stdout, * ) 
         write( stdout, * ) 'Smooth part + Augmentatio Part: '
         do iss=1,nspin
            rsumg(iss)=omega*real(rhog(1,iss))
            rsumr(iss)=SUM(rhor(1:nnr,iss))*omega/dfloat(nr1*nr2*nr3)
         end do
#ifdef __PARA
         if (ng0.ne.2) then
! in the parallel case, only one processor has G=0 ! 
            do iss=1,nspin
               rsumg(iss)=0.0
            end do
         end if
         call reduce(nspin,rsumg)
         call reduce(nspin,rsumr)
#endif
         if (nspin.eq.1) then
            WRITE( stdout,1) rsumg(1),rsumr(1)
         else
            if(iprsta.ge.3)                                             &
     &          WRITE( stdout,2)  rsumg(1),rsumg(2),rsumr(1),rsumr(2)
            WRITE( stdout,1) rsumg(1)+rsumg(2),rsumr(1)+rsumr(2)
         endif
      endif
!
    2 format(//' subroutine rhoofr: total integrated electronic',       &
     &     ' density'/' in g-space =',f10.6,2x,f10.6,4x,                &
     &     ' in r-space =',f10.6,2x,f10.6)
    1 format(//' subroutine rhoofr: total integrated electronic',       &
     &     ' density'/' in g-space =',f10.6,4x,                         &
     &     ' in r-space =',f10.6)
!

      write( stdout, * ) 
      write( stdout , * ) 'State Written : ' ,iwf, 'to unit',ndwwf
      write( stdout, * ) 

      deallocate( psi  )
      deallocate( psis ) 

      call stop_clock(' rhoiofr ')
!
      return
      end
!
!-----------------------------------------------------------------------
