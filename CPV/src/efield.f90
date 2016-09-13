!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE efield_module

  USE kinds, ONLY : DP

  IMPLICIT NONE
  SAVE

  logical      :: tefield  = .FALSE.
  logical      :: tefield2 = .FALSE.
  integer      :: epol     = 3 !direction electric field
  real(kind=DP) :: efield   = 0.d0 !intensity electric field
  real(kind=DP)  :: efield2 =0.d0 
  real(kind=DP)    evalue!strenght of electric field
  real(kind=DP)  evalue2
  integer epol2,ipolp2
  integer         ipolp  !direction of electric field

  real(kind=DP) :: pberryel = 0.0d0, pberryion = 0.0d0
  real(kind=DP) :: pberryel2 = 0.0d0, pberryion2 = 0.0d0

!***
!***  Berry phase
!***
      integer, allocatable:: ctable(:,:,:)!correspondence tables for diff. polarization
      integer, allocatable:: ctabin(:,:,:)!inverse correspondence table
      complex(DP), allocatable:: qmat(:,:)!inverse of matrix Q, for Barry's phase
      complex(DP), allocatable:: gqq(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr
      complex(DP), allocatable:: gqqm(:,:,:,:)! the same with exp(-iGr)
      complex(DP), allocatable:: gqq0(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr, at Gamma
      complex(DP), allocatable:: gqqm0(:,:,:,:)! the same with exp(-iGr), at Gamma
      complex(DP), allocatable:: df(:)
      integer, allocatable:: ctable2(:,:,:)!correspondence tables for diff. polarization
      integer, allocatable:: ctabin2(:,:,:)!inverse correspondence table
      complex(DP), allocatable:: qmat2(:,:)!inverse of matrix Q, for Barry's phase
      complex(DP), allocatable:: gqq2(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr
      complex(DP), allocatable:: gqqm2(:,:,:,:)! the same with exp(-iGr)
      complex(DP), allocatable:: gqq02(:,:,:,:)!factors int beta_Ri^*beta_Rj exp(iGr)dr, at Gamma
      complex(DP), allocatable:: gqqm02(:,:,:,:)! the same with exp(-iGr), at Gamma
      complex(DP) detq
      complex(DP) detq2
      real(DP) cdzp(3),cdzm(3), cdz0(3)!centers of ionic charges

!for parallelization for direcions 1 and 2
      integer :: n_g_missing_p(2)!number of g vector with correspondence G-->G+1 is missing
      integer :: n_g_missing_m(2)!number of g vector with correspondence G-->G-1 is missing
      integer, allocatable :: whose_is_g(:) !correspondence G(plane waves, global) ---> processor
      integer, allocatable :: ctable_missing_1(:,:,:)!correspondence G(plane waves local)--> array for mpi_alltoall
                                              !n_g_missing*nproc
      integer, allocatable :: ctable_missing_rev_1(:,:,:)!missing_g --> G (plane waves local)
      integer, allocatable :: ctable_missing_2(:,:,:)!correspondence G(plane waves local)--> array for mpi_alltoall
                                              !n_g_missing*nproc
      integer, allocatable :: ctable_missing_rev_2(:,:,:)!missing_g --> G (plane waves local)
      integer, allocatable :: ctabin_missing_1(:,:,:)!correspondence G(plane waves local)--> array for mpi_alltoall
                                              !n_g_missing*nproc
      integer, allocatable :: ctabin_missing_rev_1(:,:,:)!missing_g --> G (plane waves local)
      integer, allocatable :: ctabin_missing_2(:,:,:)!correspondence G(plane waves local)--> array for mpi_alltoall
                                              !n_g_missing*nproc
      integer, allocatable :: ctabin_missing_rev_2(:,:,:)!missing_g --> G (plane waves local)

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
    COMPLEX(DP), INTENT(IN)  :: eigr(:,:)
    REAL(DP), INTENT(IN)  :: tau0(:,:)
    if(ionode) write(stdout,'("Initialize Berry phase electric field")')
    ipolp = epol
    evalue = efield 
!set up for parallel calculations

#if defined(__MPI)
    call find_whose_is_g
    call gtable_missing
    call gtable_missing_inv
#endif

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
    COMPLEX(DP), INTENT(IN)  :: eigr(:,:)
    call qqupdate(eigr,gqqm0,gqq,gqqm,ipolp)
    RETURN
  END SUBROUTINE efield_update


  SUBROUTINE allocate_efield( ngw, ngw_g, nx, nhx, nax, nsp )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, ngw_g, nx, nhx, nax, nsp
      allocate( ctable(ngw,2,3))
      allocate( ctabin(ngw,2,3))
      allocate( qmat(nx,nx))
      allocate( gqq(nhx,nhx,nax,nsp))
      allocate( gqqm(nhx,nhx,nax,nsp))
      allocate( df(ngw))
      allocate( gqq0(nhx,nhx,nax,nsp))
      allocate( gqqm0(nhx,nhx,nax,nsp))
      allocate( whose_is_g(ngw_g))
   
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
    IF( allocated( whose_is_g) ) deallocate(whose_is_g)
    IF( allocated( ctable_missing_1) ) deallocate( ctable_missing_1)
    IF( allocated( ctable_missing_2) ) deallocate( ctable_missing_2)
    IF( allocated( ctable_missing_rev_1) ) deallocate( ctable_missing_rev_1)
    IF( allocated( ctable_missing_rev_1) ) deallocate( ctable_missing_rev_2)
    IF( allocated( ctabin_missing_1) ) deallocate( ctabin_missing_1)
    IF( allocated( ctabin_missing_2) ) deallocate( ctabin_missing_2)
    IF( allocated( ctabin_missing_rev_1) ) deallocate( ctabin_missing_rev_1)
    IF( allocated( ctabin_missing_rev_1) ) deallocate( ctabin_missing_rev_2)
    RETURN
  END SUBROUTINE deallocate_efield


  SUBROUTINE berry_energy( enb, enbi, bec, cm, fion )
    USE ions_positions, ONLY: tau0
    USE control_flags, ONLY: tfor, tprnfor
    IMPLICIT NONE
    real(DP), intent(out) :: enb, enbi
    real(DP) :: bec(:,:)
    real(DP) :: fion(:,:)
    complex(DP) :: cm(:,:)
    real(dp), external :: enberry

    call qmatrixd(cm,bec,ctable(1,1,ipolp),gqq,qmat,detq,ipolp)
    enb =  enberry( detq, ipolp )
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
    complex(DP), intent(out) :: c2(:), c3(:)
    complex(DP), intent(in) :: cm(:,:)
    REAL(DP) :: rhos(:,:)
    real(DP) :: bec(:,:)
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
    COMPLEX(DP), INTENT(IN)  :: eigr(:,:)
    REAL(DP), INTENT(IN)  :: tau0(:,:)
    if(ionode) write(stdout,'("Initialize Berry phase electric field")')
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
    COMPLEX(DP), INTENT(IN)  :: eigr(:,:)
    call qqupdate(eigr,gqqm02,gqq2,gqqm2,ipolp2)
    RETURN
  END SUBROUTINE efield_update2


  SUBROUTINE allocate_efield2( ngw, nx, nhx, nax, nsp )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx, nhx, nax, nsp
      allocate( ctable2(ngw,2,3))
      allocate( ctabin2(ngw,2,3))
      allocate( qmat2(nx,nx))
      allocate( gqq2(nhx,nhx,nax,nsp))
      allocate( gqqm2(nhx,nhx,nax,nsp))
      allocate( gqq02(nhx,nhx,nax,nsp))
      allocate( gqqm02(nhx,nhx,nax,nsp))
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
    USE ions_positions, ONLY: tau0
    USE control_flags, ONLY: tfor, tprnfor
    IMPLICIT NONE
    real(DP), intent(out) :: enb, enbi
    real(DP) :: bec(:,:)
    real(DP) :: fion(:,:)
    complex(DP) :: cm(:,:)
    real(dp), external :: enberry

    call qmatrixd(cm,bec,ctable2(1,1,ipolp2),gqq2,qmat2,detq2,ipolp2)
    enb =  enberry( detq2, ipolp2 )
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
    complex(DP), intent(out) :: c2(:), c3(:)
    complex(DP), intent(in) :: cm(:,:)
    REAL(DP) :: rhos(:,:)
    real(DP) :: bec(:,:)
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
