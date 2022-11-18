! Copyright (C) 2017 Quantum ESPRESSO group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 

MODULE dftd3_qe

  USE dftd3_sizes
  USE dftd3_common
  USE dftd3_core
  USE dftd3_api

  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: dftd3, dftd3_in, dftd3_xc, dftd3_pbc_gdisp, dftd3_printout, dftd3_clean
  PUBLIC :: dftd3_pbc_gdisp_new, dftd3_pbc_hdisp, print_dftd3_hessian 
  SAVE
  !
  type(dftd3_calc) :: dftd3
  type(dftd3_input):: dftd3_in
   
  CONTAINS

    !---------------------------------------------------------------------------
    SUBROUTINE print_dftd3_hessian( mat, n, label )
      USE kinds,            ONLY: DP
      USE io_global,        ONLY: stdout
    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      COMPLEX(DP), INTENT(IN) :: mat(3,n,3,n)
      CHARACTER(LEN=*), INTENT(IN) :: label
      !
      INTEGER :: i, iat, ixyz, j, jat, jxyz
      COMPLEX(DP), ALLOCATABLE :: buffer(:)
      CHARACTER(LEN=256):: formt, filout, string
      !
      WRITE(formt,'(A,I9,A)') '(', 2*3*n, 'f24.16)'
      WRITE(filout, '(A,A,A)')   'dynamical.',TRIM(label),'.dat'
      !
      WRITE( stdout, '(/,5x,2A)') 'Writing Hessian on file ', TRIM(filout)
      !
      ALLOCATE( buffer(3*n) )
      !
      OPEN (unit = 1, file = TRIM(filout), status = 'unknown')
      !
      WRITE(1,'(A)') 'Hessian matrix of the Grimme-D3 dispersion term'
      WRITE(1,'(A,5I9,3x,L)') 'System: '
      !
      DO i = 1, 3*n         ! 1 2 3 4 5 6 7 8 9 ... 3*nat
        iat  = (i+2)/3        ! 1 1 1 2 2 2 3 3 3 ... nat 
        ixyz = i - 3* (iat-1) ! 1 2 3 1 2 3 1 2 3 ... 3 
        !
        DO j = 1, 3*n
          jat  = (j+2)/3 
          jxyz = j - 3* (jat-1) 
          buffer(j) = mat(ixyz,iat,jxyz,jat)
        END DO 
        !  
        WRITE(1, formt ) buffer(1:3*n)
        !
      END DO 
      !
      CLOSE (1)
      !
      DEALLOCATE(buffer)
      !
      RETURN
      !
    END SUBROUTINE print_dftd3_hessian   

    !> Clean memory after a dftd3 calculator is run.
    !!
    SUBROUTINE dftd3_clean( )

      IF (ALLOCATED(dftd3%c6ab)) DEALLOCATE(dftd3%c6ab)
      IF (ALLOCATED(dftd3%mxc) ) DEALLOCATE(dftd3%mxc)
      IF (ALLOCATED(dftd3%r0ab)) DEALLOCATE(dftd3%r0ab)

    END SUBROUTINE dftd3_clean
    
    !> Convert XC labels from QE to those used by DFT-D3
    FUNCTION dftd3_xc ( dft )
      CHARACTER(LEN=*), INTENT(in) :: dft
      CHARACTER(LEN=256) :: dftd3_xc
      CHARACTER(LEN=1), EXTERNAL :: lowercase
      integer :: i
       
      dftd3_xc = ''
      DO i=1,LEN_TRIM(dft)
         dftd3_xc(i:i) = lowercase(dft(i:i))
      END DO
      IF( TRIM(dftd3_xc)=="bp")      dftd3_xc="b-p"
      IF( TRIM(dftd3_xc)=="blyp")    dftd3_xc="b-lyp"
      IF( TRIM(dftd3_xc)=="b3lyp")   dftd3_xc="b3-lyp"
      IF( TRIM(dftd3_xc)=="hse")     dftd3_xc="hse06"
      IF( TRIM(dftd3_xc)=="pw86pbe") dftd3_xc="rpw86-pbe"
      IF( TRIM(dftd3_xc)=="olyp")    dftd3_xc="o-lyp"

    END FUNCTION dftd3_xc
  
  !> Calculates forces and stress
  !! Added interface routine to original Aradi's interface
  !! for calculating force and stress separately from the energy.
  subroutine dftd3_pbc_gdisp(this, coords, izp, latvecs, &
                            force_dftd3, stress_dftd3)

    type(dftd3_calc), intent(in) :: this
    real(wp), intent(in) :: coords(:,:)
    integer, intent(in) :: izp(:)
    real(wp), intent(in) :: latvecs(:,:)
    real(wp), intent(out) :: force_dftd3(:,:), stress_dftd3(:,:)
    integer :: natom
    real(wp) :: s6, s18, rs6, rs8, rs10, alp6, alp8, alp10
    real(wp) :: e6, e8, e10, e12, e6abc, gnorm, disp2
    real(wp) :: rtmp3(3)
    integer :: rep_cn(3), rep_vdw(3)

    natom = size(coords, dim=2)
    s6 = this%s6
    s18 = this%s18
    rs6 = this%rs6
    rs8 = this%rs18
    rs10 = this%rs18
    alp6 = this%alp
    alp8 = alp6 + 2.0_wp
    alp10 = alp8 + 2.0_wp

    call set_criteria(this%rthr, latvecs, rtmp3)
    rep_vdw(:) = int(rtmp3) + 1
    call set_criteria(this%cn_thr, latvecs, rtmp3)
    rep_cn(:) = int(rtmp3) + 1
   
    force_dftd3(:,:) = 0.0_wp
    call pbcgdisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
        & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%noabc, this%numgrad, this%version, force_dftd3, disp2, gnorm, &
        & stress_dftd3, latvecs, rep_vdw, rep_cn, this%rthr, .false., this%cn_thr)
    ! Note, the stress variable in pbcgdisp contains the *lattice derivatives*
    ! on return, so it needs to be converted to obtain the stress tensor.
    stress_dftd3(:,:) = -matmul(stress_dftd3, transpose(latvecs))&
        & / abs(determinant(latvecs))  

  end subroutine dftd3_pbc_gdisp

  subroutine dftd3_pbc_gdisp_new(this, coords, izp, latvecs, &
                            force_dftd3, rep_cn_, rep_vdw_)

    type(dftd3_calc), intent(in) :: this
    real(wp), intent(in) :: coords(:,:)
    integer, intent(in) :: izp(:)
    real(wp), intent(in) :: latvecs(:,:)
    real(wp), intent(out) :: force_dftd3(:,:)
    integer, optional, intent(out) :: rep_cn_(3), rep_vdw_(3)
    integer :: natom
    real(wp) :: s6, s18, rs6, rs8, rs10, alp6, alp8, alp10
    real(wp) :: gnorm, disp2
    real(wp) :: rtmp3(3)
    integer :: rep_cn(3), rep_vdw(3)
    real(wp), allocatable :: force_supercell_dftd3(:,:,:,:,:)

    natom = size(coords, dim=2)
    s6 = this%s6
    s18 = this%s18
    rs6 = this%rs6
    rs8 = this%rs18
    rs10 = this%rs18
    alp6 = this%alp
    alp8 = alp6 + 2.0_wp
    alp10 = alp8 + 2.0_wp

    call set_criteria(this%rthr, latvecs, rtmp3)
    rep_vdw(:) = int(rtmp3) + 1
    call set_criteria(this%cn_thr, latvecs, rtmp3)
    rep_cn(:) = int(rtmp3) + 1
    if(present(rep_cn_)) rep_cn_(:) = rep_cn(:)
    if(present(rep_vdw_)) rep_vdw_(:) = rep_vdw(:)
   
    force_dftd3(:,:) = 0.0_wp

    allocate( force_supercell_dftd3(-rep_vdw(3):rep_vdw(3),-rep_vdw(2):rep_vdw(2),-rep_vdw(1):rep_vdw(1),3,natom))
    call pbcgdisp_new(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
        & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%noabc, this%numgrad, this%version, force_dftd3, disp2, gnorm, &
        & latvecs, rep_vdw, rep_cn, this%rthr, .true., this%cn_thr, &
        & 0.0000000010d0, 2, 1, 1, force_supercell_dftd3)
    deallocate( force_supercell_dftd3 )

  end subroutine dftd3_pbc_gdisp_new

  ! Calculates the dispersion contribution to Hessian
  subroutine dftd3_pbc_hdisp(this, stdout, step, coords, izp, latvecs, rep_cn, rep_vdw, hess_dftd3_, q_gamma )
  implicit none
    type(dftd3_calc), intent(in) :: this
    integer, intent(in) :: stdout 
    real(wp), intent(in) :: step             ! Bohr
    real(wp), intent(inout) :: coords(:,:)   ! Bohr
    ! in practice this is a (in). (inout) is only to allow differentiation without further allocations
    integer, intent(in) :: izp(:)
    real(wp), intent(in) :: latvecs(:,:)     ! Bohr
    integer, intent(in)  :: rep_cn(3), rep_vdw(3)
    logical, intent(in) :: q_gamma
    real(wp), intent(out), target, contiguous :: hess_dftd3_(:,:,:,:,:,:,:)
    real(wp), pointer :: hess_dftd3(:,:,:,:,:,:,:)
    integer, allocatable :: ns(:)

    real(wp) :: s6, s18, rs6, rs8, rs10, alp6, alp8, alp10
    real(wp) :: gnorm, disp2, hnorm
    integer :: natom
    integer :: iat, ixyz, istep 
    integer :: irep, jrep, krep 
    real(wp) :: rtmp3(3), stress_dftd3(3,3)
    real(wp), allocatable :: force_dftd3(:,:), force_supercell_dftd3(:,:,:,:,:)

    write(stdout, '(/,5x,A)' ) 'Three body terms in DFT-D3 disabled for Hessian calculation'

    natom = size(coords, dim=2)

    ns = shape(hess_dftd3_)
    hess_dftd3( -ns(1)/2:ns(1)/2,-ns(2)/2:ns(2)/2,-ns(3)/2:ns(3)/2, 1:ns(4),1:ns(5),1:ns(6),1:ns(7) ) => hess_dftd3_

    if(size(hess_dftd3, dim=5) .ne. natom ) Call errore('dftd3_pbc_hdisp', "Wrong Hessian dimensions", 1)
    if(size(hess_dftd3, dim=7) .ne. natom ) Call errore('dftd3_pbc_hdisp', "Wrong Hessian dimensions", 2)

    s6 = this%s6
    s18 = this%s18
    rs6 = this%rs6
    rs8 = this%rs18
    rs10 = this%rs18
    alp6 = this%alp
    alp8 = alp6 + 2.0_wp
    alp10 = alp8 + 2.0_wp

    allocate( force_dftd3(3,natom) )
    if(.not.q_gamma) allocate( & 
      force_supercell_dftd3(-rep_vdw(3):rep_vdw(3),-rep_vdw(2):rep_vdw(2),-rep_vdw(1):rep_vdw(1),3,natom))

    hess_dftd3(:,:,:,:,:,:,:) = 0.0_wp
    !
    ! Hessian update:
    !   h = - 1/2 * ( f+ - f- ) / step
    !   f+ = -2 * force_from_pbcgdisp(+)
    !   f- = -2 * force_from_pbcgdisp(-)
    !   h = ( force_supercell_dftd3(+) - force_supercell_dftd3(-) ) / step
    !
    do iat = 1, natom
      do ixyz = 1, 3  
          !
          if(q_gamma) then ! all periodic images of iat, ixyz will be displaced inside pbcgdisp 
            !
            write(stdout, '(5x,A,2I4,A)' ) 'Displacement step: ', iat, ixyz, '   (+) ' 
            coords(ixyz,iat)=coords(ixyz,iat)+step
            force_dftd3(1:3,1:natom) = 0.0_wp
            call pbcgdisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
                & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
                & .true., .false., this%version, force_dftd3, disp2, gnorm, &
                & stress_dftd3, latvecs, rep_vdw, rep_cn, this%rthr, .true., this%cn_thr) 
            hess_dftd3(0,0,0, ixyz,iat,1:3,1:natom) = hess_dftd3(0,0,0, ixyz,iat,1:3,1:natom) + force_dftd3(1:3,1:natom)
            !
            write(stdout, '(5x,A,2I4,A)' ) 'Displacement step: ', iat, ixyz, '   (-) ' 
            coords(ixyz,iat)=coords(ixyz,iat)-2*step
            force_dftd3(1:3,1:natom) = 0.0_wp
            call pbcgdisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
                & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
                & .true., .false., this%version, force_dftd3, disp2, gnorm, &
                & stress_dftd3, latvecs, rep_vdw, rep_cn, this%rthr, .true., this%cn_thr) 
            hess_dftd3(0,0,0, ixyz,iat,1:3,1:natom) = hess_dftd3(0,0,0, ixyz,iat,1:3,1:natom) - force_dftd3(1:3,1:natom)
            !
            coords(ixyz,iat)=coords(ixyz,iat)+step
            !
          else ! only iat, ixyz in the unit cell will be displaced inside pbcgdisp  
            !
            do istep = -1, 1, 2
              write(stdout, '(5x,A,3I4)' ) 'Displacement step: ', iat, ixyz, istep
              force_dftd3(:,:) = 0.0_wp ! this is not initialized in pbcgdisp 
              !force_supercell_dftd3(:,:,:,:,:) = 0.0_wp ! this is initialized in pbcgdisp
              call pbcgdisp_new(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
                  & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
                  & .true., .false., this%version, force_dftd3, disp2, gnorm, &
                  & latvecs, rep_vdw, rep_cn, this%rthr, .true., this%cn_thr, &
                  & step, iat, ixyz, istep, force_supercell_dftd3)
              do krep = -rep_vdw(3), rep_vdw(3) 
                do jrep = -rep_vdw(2), rep_vdw(2) 
                  do irep = -rep_vdw(1), rep_vdw(1) 
                    hess_dftd3(krep,jrep,irep, ixyz,iat,1:3,1:natom) = hess_dftd3(krep,jrep,irep, ixyz,iat,1:3,1:natom) & 
                                                   + dble(istep) * force_supercell_dftd3(krep,jrep,irep,1:3,1:natom)  
                  end do 
                end do 
              end do 
              !
            end do ! istep
            !
          end if ! q_gamma
          !
          hnorm=sum(abs(hess_dftd3(0,0,0, ixyz,iat, 1:3,1:natom)))
          write(*,*)'|H(unit cell)| =',hnorm
          !
      end do ! ixyz
    end do ! iat
    hess_dftd3(:,:,:,:,:,:,:) = hess_dftd3(:,:,:,:,:,:,:) / step 

    deallocate( force_dftd3 )
    if(.not.q_gamma) deallocate( force_supercell_dftd3 )

    return

  end subroutine dftd3_pbc_hdisp


  subroutine dftd3_printout(this, input_dftd3, stdout, nsp, atm, nat, ityp,&
                  tau, at, alat )
    !
    ! Added subroutine to original interface, which computes things for printout. 
    ! Could also be done in subroutine dftd3_init.
    !
    type(dftd3_calc), intent(inout) :: this
    type(dftd3_input), intent(in) :: input_dftd3
    integer, intent(in) :: stdout
    integer, intent(in) :: nsp
    integer, intent(in) :: nat
    character(len=*), intent(in) :: atm(nsp)
    integer, intent(in) :: ityp(nat)
    real*8, intent(in) :: tau(3,nat)
    real*8, intent(in) :: at(3,3)
    real*8, intent(in) :: alat

    integer :: i,j,ata,z,iz(nat)
    real*8  :: cn(nat),rtmp3(3),c6,c8,dum,x
    real*8 :: tau_(3,nat)
    real*8 :: at_(3,3)
    !
    !
    write(stdout,'( /, 5X, "--------------------------------------------" )' )
    if ( input_dftd3%threebody ) then
       write(stdout,'(    5X, "DFT-D3 Dispersion Correction (3-body terms):")')
    else
       write(stdout,'(    5X, "DFT-D3 Dispersion Correction (no 3-body):")')
    end if
    write(stdout,'(    5X, "--------------------------------------------" , &
                  & /, 5X, "  Reference C6 values for interpolation: ",/, &
                  & /, 5X, "    atom   Coordination number   C6" )')
    !
    do i=1,94
        do ata=1,nsp
          if(i.ne.get_atomic_number(atm(ata))) cycle
            do j=1,maxc
              if(this%c6ab(i, i, j, j, 1).gt.0) then
              write(stdout,'( 9X, A3 , 7X , F6.3, 9X, F8.2)' ) &
                    atm(ata), this%c6ab(i,i,j,j,2), this%c6ab(i,i,j,j,1)*2.d0
              endif
            end do
        end do
    end do

    write(stdout,'( /, 7X, "Values used:",/, &
              & /, 7X, "  atom   Coordination number  R0_AB[au]  C6      C8" )')
    !
    do ata=1,nat
        iz(ata)=get_atomic_number(trim(atm(ityp(ata))))
    end do
    !
    tau_(:,:) = tau(:,:) * alat
    at_(:,:) = at(:,:) * alat
    CALL set_criteria(this%rthr, at_, rtmp3)
    this%rep_vdw(:) = int(rtmp3) + 1
    CALL set_criteria(this%cn_thr, at_, rtmp3)
    this%rep_cn(:) = int(rtmp3) + 1
    CALL pbcncoord(nat, rcov, iz, tau_, cn, at_, this%rep_cn, this%cn_thr)
    !    
    x = 0.d0
    do ata = 1, nat
        z = get_atomic_number(trim(atm(ityp(ata))))
        CALL getc6(maxc,max_elem,this%c6ab,this%mxc, &
                              &  iz(ata),iz(ata),cn(ata),cn(ata),c6)
        c8=r2r4(iz(ata))**2*3.0d0*c6
        do j = 1, nat
          CALL getc6(maxc,max_elem,this%c6ab,this%mxc, &
                                & iz(ata),iz(j),cn(ata),cn(j),dum)
          x = x + dum
        enddo
        write(stdout,'( 9X, A3 , 7X, F6.3, 10X, F7.3, F10.2, F10.2)') &
                atm(ityp(ata)), cn(ata), 0.5*this%r0ab(z,z),c6*2.d0,c8*2.d0
    end do
    write(stdout,'(/, 9X, "Molecular C6 ( Ry / a.u.^6 ) = ",F12.2,/)') x*2.d0

  end subroutine dftd3_printout


END MODULE dftd3_qe
