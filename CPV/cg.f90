!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE cg_module

  IMPLICIT NONE
  SAVE

      logical      :: tcg        = .false.   ! se vero fa gradiente coniugato
      integer      :: maxiter    = 20      !numero massimo interazioni c.g.
      real(kind=8) :: etresh    = 1.d-5    !soglia convergenza c.g.
      real(kind=8) :: passop    =0.3d0    !passetto per gradiente coniugato

!***
!***  Conjugate Gradient
!***
      real(kind=8)  esse,essenew !fattori cg
      COMPLEX(kind=8), ALLOCATABLE :: gi(:,:)!coniugati
      COMPLEX(kind=8), ALLOCATABLE :: hi(:,:)!gradienti di ricerca
      COMPLEX(kind=8), ALLOCATABLE :: c0old(:,:)!vecchie funzioni d'onda, per estrapolazione
      COMPLEX(kind=8), ALLOCATABLE :: hpsi(:,:) !termini H|Psi_i>
      real(kind=8)  ene0,ene1,dene0,enever,enesti !energie per minimizazzione lineare lungo hi
      real(kind=8)  passof,passov !passo effettivo stimato durante minimizzazione
      integer itercg !numero iterazione
      logical ltresh!flag per convergenza su energia
      real(kind=8) passo!passo per arrivare a  minimo
      real(kind=8) etotnew,etotold!per vedere convergenza
      real(kind=8) spasso!segno passetto
      logical tcutoff!convergenza energia per togliere cutoff
      real(kind=8), ALLOCATABLE :: emme(:,:)!matrice per convergenza spinta
      logical restartcg!se vero ricomincia gradiente coniugato da steepest descend
      integer numok!numero volte differenza energia sotto treshhold
      real(kind=8) pcnum,pcden!per calcolare preconditioning
      integer iter3!per tentativi ciclo3
      real(kind=8) ebanda!energia banda per preconditioning
      logical ene_ok!se vero l'energia ha passato il test, non deve ricalcolare
      integer ninner_ef!per conteggio ciclo interno, usato per stati eccitati


CONTAINS


  SUBROUTINE cg_init( tcg_ , maxiter_ , etresh_ , passop_ )
    USE kinds, ONLY: dbl
    LOGICAL, INTENT(IN) :: tcg_
    INTEGER, INTENT(IN) :: maxiter_
    REAL(dbl), INTENT(IN) :: etresh_ , passop_
    RETURN
  END SUBROUTINE cg_init

  SUBROUTINE cg_info()
    USE io_global, ONLY: stdout 
    if(tcg) then
       write (stdout,400) maxiter,etresh,passop                         
    endif
400 format (/4x,'====================================='                          &
   &        /4x,'|  GRADIENTE CONIUGATO              |'                          &
   &        /4x,'====================================='                          &
   &        /4x,'| iterazioni   =',i10,'          |'                             &
   &        /4x,'| etresh       =',f10.5,' a.u.     |'                           &
   &        /4x,'| passop       =',f10.5,' a.u.     |'                           &
   &        /4x,'=====================================')
    RETURN
  END SUBROUTINE cg_info


  SUBROUTINE allocate_cg( ngw, nx )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx
    allocate(hi(ngw,nx))!poi bisogna fare che uno semplicemnte punti su cm
    allocate(gi(ngw,nx))
    allocate(c0old(ngw,nx))
    allocate( emme(nx,nx))
    allocate( hpsi(ngw,nx))
    RETURN
  END SUBROUTINE allocate_cg

  SUBROUTINE deallocate_cg( )
    IMPLICIT NONE
    IF( ALLOCATED( hi ) ) deallocate(hi )
    IF( ALLOCATED( gi ) ) deallocate(gi )
    IF( ALLOCATED( c0old ) ) deallocate(c0old )
    IF( ALLOCATED( emme ) ) deallocate( emme )
    IF( ALLOCATED( hpsi ) ) deallocate( hpsi )
    RETURN
  END SUBROUTINE deallocate_cg

  SUBROUTINE cg_update( tfirst, nfi, c0 )
    use gvecw, only: ngw
    use electrons_base, only: n => nbsp
    IMPLICIT NONE
    COMPLEX(kind=8) :: c0( :, :, :, : )
    INTEGER :: nfi
    LOGICAL :: tfirst
    INTEGER :: i, ig
    if(.not. tfirst.and.(mod(nfi,10).ne.1)) then
      call DSWAP(2*ngw*n,c0,1,c0old,1)
      do i=1,n
        do ig=1,ngw
          c0(ig,i,1,1)=-c0(ig,i,1,1)+2.d0*c0old(ig,i)
        enddo
      enddo
    else
      do i=1,n
        do ig=1,ngw
          c0old(ig,i)=c0(ig,i,1,1)
        enddo
      enddo
    endif
    RETURN
  END SUBROUTINE cg_update


END MODULE
