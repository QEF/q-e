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
  SAVE
  !
  type(dftd3_calc) :: dftd3
  type(dftd3_input):: dftd3_in
  real(wp) :: energy_dftd3
   
  CONTAINS

    !> Clean memory after a dftd3 calculator is run.
    !!
    SUBROUTINE dftd3_clean( )

      IF (ALLOCATED(dftd3%c6ab)) DEALLOCATE(dftd3%c6ab)
      IF (ALLOCATED(dftd3%mxc) ) DEALLOCATE(dftd3%mxc)
      IF (ALLOCATED(dftd3%r0ab)) DEALLOCATE(dftd3%r0ab)

    END SUBROUTINE dftd3_clean
    
    !> Convert XC labels from QE to those used by DFT-D3
    FUNCTION dftd3_xc ( dft )
      IMPLICIT NONE
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


  subroutine dftd3_printout(this, input_dftd3)
    !
    ! Added subroutine to original interface, which computes things for printout. 
    ! Could also be done in subroutine dftd3_init.
    !
    USE io_global, ONLY : stdout
    USE ions_base, ONLY : nsp,atm,nat,ityp,tau
    USE cell_base, ONLY : at,alat

    type(dftd3_calc), intent(inout) :: this
    type(dftd3_input), intent(in) :: input_dftd3
    integer :: i,j,ata,z,iz(nat)
    real*8  :: cn(nat),rtmp3(3),c6,c8,dum,x
    !
    !
    write(stdout,'( /, 5X, "--------------------------------------------" , &
                  & /, 5X, "Parameters for DFT-D3 Dispersion Correction:" , &
                  & /, 5X, "--------------------------------------------" , &
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
    tau(:,:) = tau(:,:) * alat
    at(:,:) = at(:,:) * alat
    CALL set_criteria(this%rthr, at, rtmp3)
    this%rep_vdw(:) = int(rtmp3) + 1
    CALL set_criteria(this%cn_thr, at, rtmp3)
    this%rep_cn(:) = int(rtmp3) + 1
    CALL pbcncoord(nat, rcov, iz, tau, cn, at, this%rep_cn, this%cn_thr)
    tau(:,:) = tau(:,:) / alat
    at(:,:) = at(:,:) / alat
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
