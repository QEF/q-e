!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine read_pseudo_ncpp (file_pseudo,zed,grid,ndmx,&
                        dft,lmax,lloc,zval,nlcc,rhoc,vnl,vpsloc,rel)
  !-----------------------------------------------------------------------
  !
  use kinds, only : DP
  use constants, only : fpi
  use radial_grids, only: radial_grid_type, do_mesh
  implicit none
  !
  ! I/O variables
  !
  type (radial_grid_type), intent(out):: grid
  !
  integer  ::    &
       ndmx, &    ! input: the mesh dimensions
       rel, &    ! input: rel=2 for spin-orbit pseudopotential
       mesh,&    ! output: the number of mesh points
       lmax,&    ! output: the maximum angular momentum
       lloc      ! output: the local potential

  real(DP) ::       &
       zed,            & ! input: the atomic charge
       zval,           & ! output: the valence charge
!       xmin,dx,        & ! output: the mesh 
!       rmax,           & ! output: the maximum mesh value
!       r(ndmx),r2(ndmx), & ! output: the mesh
!       rab(ndmx),       & ! output: derivative of the mesh
!       sqr(ndmx),       & ! output: the square root of the mesh
       vnl(ndmx,0:3,2), & ! output: the potential in numerical form
       vpsloc(ndmx),    & ! output: the local pseudopotential
       rhoc(ndmx)         ! output: the core charge
 
  real(DP) :: xmin, dx, rmax ! auxiliary mesh data

  logical :: &
       nlcc    ! output: if true the pseudopotential has nlcc

  character ::   &
       file_pseudo*20, &    ! input: the file with the pseudopotential
       dft*20              ! output: the type of xc

  integer :: &
       ios, i, l, k, ir, iunps, nbeta,nlc,nnl   

  real(DP) :: &
       vnloc, a_core, b_core, alfa_core, &
       cc(2),alpc(2),alc(6,0:3),alps(3,0:3)
  real(DP), external :: qe_erf 

  logical :: &
       bhstype, numeric

  character(len=3)  cdum

  iunps=2
  open(unit=iunps,file=file_pseudo,status='old',form='formatted', &
       err=100, iostat=ios)
100 call errore('read_pseudo_ncpp','open error on file '//file_pseudo,ios)
  !
  !     reads the starting lines
  !
  read( iunps, '(a)', end=300, err=300, iostat=ios ) dft
  if (dft(1:2).eq.'**') dft='LDA'
  if (dft(1:17).eq.'slater-pz-ggx-ggc') dft='PW'

  read ( iunps, *, err=300, iostat=ios ) cdum,  &
       zval, lmax, nlc, nnl, nlcc,  &
       lloc, bhstype

  if ( nlc.gt.2 .or. nnl.gt.3)  &
       call errore( 'read_pseudo_ncpp','Wrong nlc or nnl', 1)
  if ( nlc .lt.0 .or. nnl .lt. 0 )  &
       call errore( 'read_pseudo_ncpp','nlc or nnl < 0 ? ', 1 )
  if ( zval.le.0.0_dp )  &
       call errore( 'read_pseudo_ncpp','Wrong zval ', 1 )

  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric = nlc.le.0 .and. nnl.le.0
  if (lloc.eq.-1000) lloc=lmax

  if (.not.numeric) then
     read( iunps, *, err=300, iostat=ios )  &
          ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
     if ( abs(cc(1)+cc(2)-1.0_dp).gt.1.0e-6_dp) call errore  &
          ('read_pseudo_ncpp','wrong pseudopotential coefficients',1)
     do l = 0, lmax
        read ( iunps, *, err=300, iostat=ios ) &
             ( alps(i,l),i=1,3 ), (alc(i,l),i=1,6)
     enddo
     if (nlcc) then
        read( iunps, *, err=300, iostat=ios ) a_core,  &
             b_core, alfa_core
        if (alfa_core.le.0.0_dp)  &
             call errore('read_pseudo_ncpp','nlcc but alfa=0',1)
     endif
     if (cc(2).ne.0.0_dp.and.alpc(2).ne.0.0_dp.and.bhstype) then
        call bachel(alps,alc,1,lmax)
     endif

  endif
  !
  !     read the mesh parameters
  !
  read( iunps, *, err=300, iostat=ios ) zed, xmin, dx, mesh, nbeta
  rmax=exp( (mesh-1)*dx+xmin )/zed
  !
  !    and generate the mesh: this overwrites the mesh defined in the
  !    input parameters
  !
  call do_mesh(rmax,zed,xmin,dx,0,grid)
  if (mesh.ne.grid%mesh) &
       call errore('read_pseudo_ncpp','something wrong in mesh',1)
  !
  !    outside this routine all pseudo are numeric: construct vnl and
  !    core charge
  !    
  if (.not.numeric) then
     !
     !  obsolescent format with analytical coefficients
     !
     read( iunps, '(a)', err=300, iostat=ios ) cdum
     if (nlcc) then 
        do ir=1, mesh
           rhoc(ir)=(a_core+b_core*grid%r2(ir))*exp(-alfa_core*grid%r2(ir)) &
                *grid%r2(ir)*fpi
        enddo
     else
        rhoc=0.0_dp
     endif
     do l=0,lmax
        do ir=1,mesh
           vnloc = 0.0_dp
           do k=1,3
              vnloc = vnloc + exp(-alps(k,l)*grid%r2(ir))*  &
                   ( alc(k,l) + grid%r2(ir)*alc(k+3,l) )
           enddo
           !
           !  NB: the factor 2 converts from hartree to rydberg
           !      spin-orbit not implemented for analytical PP
           !
           vnl(ir,l,1) = 2.0_dp*vnloc
        enddo
     enddo
     do ir=1,mesh
        vpsloc(ir) = -2.0_dp*zval/grid%r(ir)* &
                   ( cc(1)*qe_erf(grid%r(ir)*sqrt(alpc(1))) &
                   + cc(2)*qe_erf(grid%r(ir)*sqrt(alpc(2))) )
     end do
  endif

  if (numeric) then
     !
     !      pseudopotentials in numerical form
     !
     do l = 0, lmax
        read( iunps, '(a)', err=300, iostat=ios ) cdum
        read( iunps, *, err=300, iostat=ios )  &
             (vnl(ir,l,1),ir=1,mesh)
        if (rel ==2 .and. l > 0) then
           !
           read( iunps, '(a)', err=300, iostat=ios ) cdum
           read( iunps, *, err=300, iostat=ios )  &
                (vnl(ir,l,2),ir=1,mesh)
        endif
     enddo
     !
     !  the local part is subtracted from the non-local part
     !  and stored in vpsloc - just for consistency
     !
     if (lloc == -1) then
        read( iunps, '(a)', err=300, iostat=ios )
        read( iunps, *, err=300, iostat=ios ) &
             (vpsloc(ir),ir=1,mesh)
     else if (rel < 2) then
        vpsloc(1:mesh) = vnl(1:mesh,lloc,1)
     else
        vpsloc(1:mesh) = (l*vnl(1:mesh,lloc,1)+(l+1)*vnl(1:mesh,lloc,2)) / &
             (2*l+1)
     end if
     do l=0,lmax
        vnl(1:mesh,l,1) = vnl(1:mesh,l,1) - vpsloc(1:mesh) 
        if (rel == 2) vnl(1:mesh,l,2) = vnl(1:mesh,l,2) - vpsloc(1:mesh) 
     end do
     !
     if(nlcc) then
        read( iunps, *, err=300, iostat=ios )  &
             ( rhoc(ir), ir=1,mesh )
        do ir=1, mesh
           rhoc(ir)=rhoc(ir)*grid%r2(ir)*fpi
        enddo
     else
        rhoc=0.0_dp
     endif

  endif
300 call errore('read_pseudo_ncpp','reading pseudofile',abs(ios))
  !
  !   all the components of the nonlocal potential beyond lmax are taken
  !   equal to the local part
  !
  do l=lmax+1,3
     vnl(1:mesh,l,1)=vpsloc(1:mesh)
     if (rel ==2 .and. l > 0) vnl(1:mesh,l,2)=vpsloc(1:mesh)
  enddo
  close(iunps)
  return
end subroutine read_pseudo_ncpp
