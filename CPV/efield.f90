!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE efield_module

  USE parameters, ONLY: natx

  IMPLICIT NONE
  SAVE

  logical      :: tefield  = .FALSE.
  integer      :: epol     = 3 !direzione campo elettrico
  real(kind=8) :: efield   = 0.d0 !intensita' campo elettrico
  real(kind=8)    evalue!strenght of electric field
  integer         ipolp  !direction of electric field

  real(kind=8) :: pberryel = 0.0d0, pberryion = 0.0d0

!***
!***  Berry phase
!***
      integer, allocatable:: ctable(:,:,:)!correspondence tables for diff. polarization
      integer, allocatable:: ctabin(:,:,:)!inverse correspondence table
      complex(kind=8), allocatable:: qmat(:,:)!inverse of matrix Q, for Barry's phase
      complex(kind=8), allocatable:: gqq(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr
      complex(kind=8), allocatable:: gqqm(:,:,:,:)! the same with exp(-iGr)
      complex(kind=8), allocatable:: gqq0(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr, at Gamma
      complex(kind=8), allocatable:: gqqm0(:,:,:,:)! the same with exp(-iGr), at Gamma
      complex(kind=8), allocatable:: df(:)
      complex(kind=8) detq
      real(kind=8) cdzp(3),cdzm(3), cdz0(3)!centers of ionic charges
      real(kind=8) taunew(3, natx)!the electric field must be always along the column direction


CONTAINS


  SUBROUTINE efield_init( epol_ , efield_ )
    USE kinds, ONLY: dbl
    REAL(dbl), INTENT(IN) :: efield_
    INTEGER, INTENT(IN)    :: epol_
    epol   = epol_
    efield = efield_
    RETURN
  END SUBROUTINE efield_init

  SUBROUTINE efield_info( )
    USE io_global, ONLY: stdout
    write (stdout,401) epol, efield
    ipolp  = epol
    evalue = efield
         
401   format (/4x,'====================================='                       &
     &        /4x,'|  CAMPO ELETTRICO                   '                       &
     &        /4x,'====================================='                       &
     &        /4x,'| direzione    =',i10,'            '                         &
     &        /4x,'| intensita    =',f10.5,' a.u.     '                         &
     &        /4x,'=====================================')

    RETURN
  END SUBROUTINE efield_info


  SUBROUTINE efield_berry_setup( eigr, tau0 )
    IMPLICIT NONE
    COMPLEX(kind=8), INTENT(IN)  :: eigr(:,:)
    REAL(kind=8), INTENT(IN)  :: tau0(:,:)
    write(6,'(''before gtable'')')
    call gtable(ipolp,ctable(1,1,ipolp))
    write(6,'(''out of gtable'')')
    call gtablein(ipolp,ctabin(1,1,ipolp))
    write(6,'(''out of gtablein'')')
    call qqberry2(gqq0,gqqm0,ipolp)!for Vanderbilt pps
    write(6,'(''out of qqberry2'')')
    call qqupdate(eigr,gqqm0,gqq,gqqm,ipolp)
    write(6,'(''out of qqupdate'')')
    call cofcharge(tau0,cdz0)
    RETURN
  END SUBROUTINE efield_berry_setup


  SUBROUTINE efield_update( eigr )
    IMPLICIT NONE
    COMPLEX(kind=8), INTENT(IN)  :: eigr(:,:)
    call qqupdate(eigr,gqqm0,gqq,gqqm,ipolp)
    RETURN
  END SUBROUTINE efield_update


  SUBROUTINE allocate_efield( ngw, nx, nhx, nas, nsp )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx, nhx, nas, nsp
      allocate( ctable(ngw,2,3))
      allocate( ctabin(ngw,2,3))
      allocate( qmat(nx,nx))
      allocate( gqq(nhx,nhx,nas,nsp))
      allocate( gqqm(nhx,nhx,nas,nsp))
      allocate( df(ngw))
      allocate( gqq0(nhx,nhx,nas,nsp))
      allocate( gqqm0(nhx,nhx,nas,nsp))
    RETURN
  END SUBROUTINE allocate_efield


  SUBROUTINE deallocate_efield( )
    IMPLICIT NONE
    IF( allocated( ctable ) )  deallocate( ctable )
    IF( allocated( ctabin ) ) deallocate( ctabin )
    IF( allocated( qmat ) ) deallocate( qmat )
    IF( allocated( gqq ) ) deallocate( gqq )
    IF( allocated( gqqm ) ) deallocate( gqqm )
    IF( allocated( df ) ) deallocate( df )
    IF( allocated( gqq0 ) ) deallocate( gqq0 )
    IF( allocated( gqqm0 ) )  deallocate( gqqm0 )
    RETURN
  END SUBROUTINE deallocate_efield


  SUBROUTINE berry_energy( enb, enbi, bec, cm, fion )
    USE uspp, ONLY: betae => vkb
    USE ions_positions, ONLY: tau0
    USE control_flags, ONLY: tfor
    IMPLICIT NONE
    real(kind=8), intent(out) :: enb, enbi
    real(kind=8) :: bec(:,:)
    real(kind=8) :: fion(:,:)
    complex(kind=8) :: cm(:,:)
    call qmatrixd(cm,bec,ctable(1,1,ipolp),gqq,qmat,detq)
    write(6,'(''out of qmatrixd'')')
    call enberry( detq, ipolp,enb)
    call berryion(tau0,fion,tfor,ipolp,evalue,enbi)
    pberryel=enb
    pberryion=enbi
    write(6,*) 'Polarizzazione',pberryel,evalue
    enb=enb*evalue
    enbi=enbi*evalue
  END SUBROUTINE berry_energy


  SUBROUTINE dforce_efield (bec,i,cm,c2,c3,rhos)
    USE uspp, ONLY: betae => vkb, deeq
    USE gvecw, ONLY: ngw
    IMPLICIT NONE
    complex(kind=8), intent(out) :: c2(:), c3(:)
    complex(kind=8), intent(in) :: cm(:)
    REAL(kind=8) :: rhos(:,:)
    real(kind=8) :: bec(:,:)
    integer :: i
    integer :: ig
    call dforceb (cm, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
    do ig=1,ngw
      c2(ig)=c2(ig)+evalue*df(ig)
    enddo
    call dforceb (cm, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
    do ig=1,ngw
      c3(ig)=c3(ig)+evalue*df(ig)
    enddo
  END SUBROUTINE dforce_efield

END MODULE efield_module
