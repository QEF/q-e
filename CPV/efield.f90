!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE efield_module

  IMPLICIT NONE
  SAVE

  logical      :: tefield  = .FALSE.
  logical      :: tefield2 = .FALSE.
  integer      :: epol     = 3 !direction electric field
  real(8) :: efield   = 0.d0 !intensity electric field
  real(8)  :: efield2 =0.d0 
  real(8)    evalue!strenght of electric field
  real(8)  evalue2
  integer epol2,ipolp2
  integer         ipolp  !direction of electric field

  real(8) :: pberryel = 0.0d0, pberryion = 0.0d0
  real(8) :: pberryel2 = 0.0d0, pberryion2 = 0.0d0

!***
!***  Berry phase
!***
      integer, allocatable:: ctable(:,:,:)!correspondence tables for diff. polarization
      integer, allocatable:: ctabin(:,:,:)!inverse correspondence table
      complex(8), allocatable:: qmat(:,:)!inverse of matrix Q, for Barry's phase
      complex(8), allocatable:: gqq(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr
      complex(8), allocatable:: gqqm(:,:,:,:)! the same with exp(-iGr)
      complex(8), allocatable:: gqq0(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr, at Gamma
      complex(8), allocatable:: gqqm0(:,:,:,:)! the same with exp(-iGr), at Gamma
      complex(8), allocatable:: df(:)
      integer, allocatable:: ctable2(:,:,:)!correspondence tables for diff. polarization
      integer, allocatable:: ctabin2(:,:,:)!inverse correspondence table
      complex(8), allocatable:: qmat2(:,:)!inverse of matrix Q, for Barry's phase
      complex(8), allocatable:: gqq2(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr
      complex(8), allocatable:: gqqm2(:,:,:,:)! the same with exp(-iGr)
      complex(8), allocatable:: gqq02(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr, at Gamma
      complex(8), allocatable:: gqqm02(:,:,:,:)! the same with exp(-iGr), at Gamma
      complex(8) detq
      complex(8) detq2
      real(8) cdzp(3),cdzm(3), cdz0(3)!centers of ionic charges


CONTAINS


  SUBROUTINE efield_init( epol_ , efield_ )
    USE kinds, ONLY: DP
    REAL(DP), INTENT(IN) :: efield_
    INTEGER, INTENT(IN)    :: epol_
    epol   = epol_
    efield = efield_
    RETURN
  END SUBROUTINE efield_init

  SUBROUTINE efield_info( )
    USE io_global, ONLY: ionode,stdout
    if(ionode) write (stdout,401) epol, efield
         
401   format (/4x,'====================================='                       &
     &        /4x,'|  BERRY PHASE ELECTRIC FIELD 1        '                       &
     &        /4x,'====================================='                       &
     &        /4x,'| direction    =',i10,'            '                         &
     &        /4x,'| intensity    =',f10.5,' a.u.     '                         &
     &        /4x,'=====================================')

    RETURN
  END SUBROUTINE efield_info


  SUBROUTINE efield_berry_setup( eigr, tau0 )
    USE io_global, ONLY: ionode,stdout
    IMPLICIT NONE
    COMPLEX(8), INTENT(IN)  :: eigr(:,:)
    REAL(8), INTENT(IN)  :: tau0(:,:)
    if(ionode) write(stdout,'(''Initialize Berry phase electric field'')')
    ipolp = epol
    evalue = efield 
    call gtable(ipolp,ctable(1,1,ipolp))
    call gtablein(ipolp,ctabin(1,1,ipolp))
    call qqberry2(gqq0,gqqm0,ipolp)!for Vanderbilt pps
    call qqupdate(eigr,gqqm0,gqq,gqqm,ipolp)
    !the following line was to keep the center of charge fixed
    !when performing molecular dynamics in the presence of an electric
    !field
    !call cofcharge(tau0,cdz0)
    RETURN
  END SUBROUTINE efield_berry_setup


  SUBROUTINE efield_update( eigr )
    IMPLICIT NONE
    COMPLEX(8), INTENT(IN)  :: eigr(:,:)
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
    USE control_flags, ONLY: tfor, tprnfor
    IMPLICIT NONE
    real(8), intent(out) :: enb, enbi
    real(8) :: bec(:,:)
    real(8) :: fion(:,:)
    complex(8) :: cm(:,:)
    call qmatrixd(cm,bec,ctable(1,1,ipolp),gqq,qmat,detq)
    call enberry( detq, ipolp,enb)
    call berryion(tau0,fion,tfor.or.tprnfor,ipolp,evalue,enbi)
    pberryel=enb
    pberryion=enbi
    enb=enb*evalue
    enbi=enbi*evalue
  END SUBROUTINE berry_energy


  SUBROUTINE dforce_efield (bec,i,cm,c2,c3,rhos)
    USE uspp, ONLY: betae => vkb, deeq
    USE gvecw, ONLY: ngw
    IMPLICIT NONE
    complex(8), intent(out) :: c2(:), c3(:)
    complex(8), intent(in) :: cm(:,:)
    REAL(8) :: rhos(:,:)
    real(8) :: bec(:,:)
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

 SUBROUTINE efield_init2( epol_ , efield_ )
    USE kinds, ONLY: DP
    REAL(DP), INTENT(IN) :: efield_
    INTEGER, INTENT(IN)    :: epol_
    epol2   = epol_
    efield2 = efield_
    RETURN
  END SUBROUTINE efield_init2

  SUBROUTINE efield_info2( )
    USE io_global, ONLY: ionode,stdout
    if(ionode) write (stdout,402) epol2, efield2
         
402   format (/4x,'====================================='                       &
     &        /4x,'|  BERRY PHASE ELECTRIC FIELD 2        '                       &
     &        /4x,'====================================='                       &
     &        /4x,'| direction    =',i10,'            '                         &
     &        /4x,'| intensity    =',f10.5,' a.u.     '                         &
     &        /4x,'=====================================')

    RETURN
  END SUBROUTINE efield_info2


  SUBROUTINE efield_berry_setup2( eigr, tau0 )
    USE io_global, ONLY: ionode,stdout
    IMPLICIT NONE
    COMPLEX(8), INTENT(IN)  :: eigr(:,:)
    REAL(8), INTENT(IN)  :: tau0(:,:)
    if(ionode) write(stdout,'(''Initialize Berry phase electric field'')')
    ipolp2 = epol2
    evalue2 = efield2 
    call gtable(ipolp2,ctable2(1,1,ipolp2))
    call gtablein(ipolp2,ctabin2(1,1,ipolp2))
    call qqberry2(gqq02,gqqm02,ipolp2)!for Vanderbilt pps
    call qqupdate(eigr,gqqm02,gqq2,gqqm2,ipolp2)
    !the following line was to keep the center of charge fixed
    !when performing molecular dynamics in the presence of an electric
    !field
    !call cofcharge(tau0,cdz0)
    RETURN
  END SUBROUTINE efield_berry_setup2


  SUBROUTINE efield_update2( eigr )
    IMPLICIT NONE
    COMPLEX(8), INTENT(IN)  :: eigr(:,:)
    call qqupdate(eigr,gqqm02,gqq2,gqqm2,ipolp2)
    RETURN
  END SUBROUTINE efield_update2


  SUBROUTINE allocate_efield2( ngw, nx, nhx, nas, nsp )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx, nhx, nas, nsp
      allocate( ctable2(ngw,2,3))
      allocate( ctabin2(ngw,2,3))
      allocate( qmat2(nx,nx))
      allocate( gqq2(nhx,nhx,nas,nsp))
      allocate( gqqm2(nhx,nhx,nas,nsp))
      allocate( gqq02(nhx,nhx,nas,nsp))
      allocate( gqqm02(nhx,nhx,nas,nsp))
    RETURN
  END SUBROUTINE allocate_efield2


  SUBROUTINE deallocate_efield2( )
    IMPLICIT NONE
    IF( allocated( ctable2 ) )  deallocate( ctable2 )
    IF( allocated( ctabin2 ) ) deallocate( ctabin2 )
    IF( allocated( qmat2 ) ) deallocate( qmat2 )
    IF( allocated( gqq2 ) ) deallocate( gqq2 )
    IF( allocated( gqqm2 ) ) deallocate( gqqm2 )
    IF( allocated( gqq02 ) ) deallocate( gqq02 )
    IF( allocated( gqqm02 ) )  deallocate( gqqm02 )
    RETURN
  END SUBROUTINE deallocate_efield2


  SUBROUTINE berry_energy2( enb, enbi, bec, cm, fion )
    USE uspp, ONLY: betae => vkb
    USE ions_positions, ONLY: tau0
    USE control_flags, ONLY: tfor, tprnfor
    IMPLICIT NONE
    real(8), intent(out) :: enb, enbi
    real(8) :: bec(:,:)
    real(8) :: fion(:,:)
    complex(8) :: cm(:,:)
    call qmatrixd(cm,bec,ctable2(1,1,ipolp2),gqq2,qmat2,detq2)
    call enberry( detq2, ipolp2,enb)
    call berryion(tau0,fion,tfor.or.tprnfor,ipolp2,evalue2,enbi)
    pberryel2=enb
    pberryion2=enbi
    enb=enb*evalue2
    enbi=enbi*evalue2
  END SUBROUTINE berry_energy2


  SUBROUTINE dforce_efield2 (bec,i,cm,c2,c3,rhos)
    USE uspp, ONLY: betae => vkb, deeq
    USE gvecw, ONLY: ngw
    IMPLICIT NONE
    complex(8), intent(out) :: c2(:), c3(:)
    complex(8), intent(in) :: cm(:,:)
    REAL(8) :: rhos(:,:)
    real(8) :: bec(:,:)
    integer :: i
    integer :: ig
    call dforceb (cm, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
    do ig=1,ngw
      c2(ig)=c2(ig)+evalue2*df(ig)
    enddo
    call dforceb (cm, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
    do ig=1,ngw
      c3(ig)=c3(ig)+evalue2*df(ig)
    enddo
  END SUBROUTINE dforce_efield2

END MODULE efield_module
