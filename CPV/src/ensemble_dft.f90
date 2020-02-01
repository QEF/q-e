!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ensemble_dft

  USE kinds, ONLY: DP

  IMPLICIT NONE
  SAVE

      logical      :: tens       = .false. ! whether to do ensemble calculations.
      logical      :: tgrand     = .false. ! whether to do grand canonical
                                           ! ensemble calculations.
      integer      :: ninner     = 0       ! number of inner loops per CP step.
      integer      :: ismear     = 2       ! type of smearing:
                                           !  1 => gaussian
                                           !  2 => fermi-dirac
                                           !  3 => hermite-delta_function
                                           !  4 => gaussian splines
                                           !  5 => cold smearing i
                                           !  6 => cold smearing ii
      real(DP) :: etemp      = 0       ! smearing temperature.
      real(DP) :: ef         = 0       ! Fermi energy (relevant if tgrand=.true.).
      
      integer :: niter_cold_restart !frequency for accuarate cold smearing (in iterations)
      real(DP) :: lambda_cold !step for cold smearing for not accurate iterations     
      INTEGER :: nrlx    ! first dimension of z0t, fmat0
      INTEGER :: nrcx    ! first two dimensions of psihpsi
!***ensemble-DFT
      real(DP), allocatable::                 z0t(:,:,:)   ! transpose of z0
      complex(DP), allocatable::             c0diag(:,:)
      real(DP), allocatable::               becdiag(:,:)
      real(DP), allocatable::                      e0(:)
      real(DP), allocatable::               fmat0(:,:,:)
      real(DP) :: gibbsfe
! variables for cold-smearing
      real(DP), allocatable ::               psihpsi(:,:,:)!it contains the matrix <Psi|H|Psi>

CONTAINS


  SUBROUTINE compute_entropy( entropy, f, nspin )
    implicit none
    real(DP), intent(out) :: entropy
    real(DP), intent(in) :: f
    integer, intent(in) :: nspin
    real(DP) :: f2
    entropy=0.0d0
    if ((f.gt.1.0d-20).and.(f.lt.(2.0/DBLE(nspin)-1.0d-20))) then
       f2=DBLE(nspin)*f/2.0d0
       entropy=-f2*log(f2)-(1.d0-f2)*log(1.d0-f2)
    end if
    entropy=-etemp*2.0d0*entropy/DBLE(nspin)
  END SUBROUTINE compute_entropy


  SUBROUTINE compute_entropy2( entropy, f, n, nspin )
    implicit none
    real(DP), intent(out) :: entropy
    real(DP), intent(in) :: f(:)
    integer, intent(in) :: n, nspin
    real(DP) :: f2
    integer :: i
    entropy=0.0d0
    do i=1,n
      if ((f(i).gt.1.0d-20).and.(f(i).lt.(2.0/DBLE(nspin)-1.0d-20))) then
        f2=DBLE(nspin)*f(i)/2.0d0
        entropy=entropy-f2*log(f2)-(1.d0-f2)*log(1.d0-f2)
      end if
    end do
    entropy=-etemp*2.0d0*entropy/DBLE(nspin)
    return
  END SUBROUTINE compute_entropy2


  SUBROUTINE compute_entropy_der( ex, fx, n, nspin )
    implicit none
    real(DP), intent(out) :: ex(:)
    real(DP), intent(in) :: fx(:)
    integer, intent(in) :: n, nspin
    real(DP) :: f2,xx
    integer :: i
    !     calculation of the entropy derivative at x
    do i=1,n
    if ((fx(i).gt.1.0d-200).and.(fx(i).lt.(2.0/DBLE(nspin)-1.0d-200))) then
      ex(i)=(log((2.0d0/DBLE(nspin)-fx(i))/fx(i)))
    else if (fx(i).le.1.0d-200) then
      xx=1.0d-200
      ex(i)=log(2.0d0/DBLE(nspin)/xx-1)
    else
      !      the calculation of ex_i is done using ex_i=-log(mf/(1-f_i)-1)
      !                                      instead of ex_i=log(mf/f_i-1)
      !      to avoid numerical errors
      xx=1.0d-200
      ex(i)=-log(2.0d0/DBLE(nspin)/xx-1)
    end if
    end do

    return
  END SUBROUTINE compute_entropy_der



  SUBROUTINE id_matrix_init( descla, nspin )
    ! initialization of the matrix identity
    USE descriptors
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nspin
    TYPE(la_descriptor), INTENT(IN) :: descla( nspin )
    INTEGER :: is, i, ii
    INTEGER :: np, me
    z0t(:,:,:)=0.0d0
    do is = 1, nspin
       np = descla( is )%npc * descla( is )%npr
       me = descla( is )%mype
       IF( descla( is )%active_node > 0 ) THEN
          ii = me + 1
          DO i = 1, descla( is )%nrl
             z0t( i, ii , is ) = 1.d0
             ii = ii + np
          END DO
       END IF
    end do
    RETURN
  END SUBROUTINE id_matrix_init


  SUBROUTINE h_matrix_init( descla, nspin )
    ! initialization of the psihpsi matrix 
    USE descriptors
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nspin
    TYPE(la_descriptor), INTENT(IN) :: descla( nspin )
    INTEGER :: is, i, nr
      psihpsi(:,:,:)=0.0d0
      do is = 1, nspin
         IF( descla( is )%active_node > 0 ) THEN
            !
            nr = descla( is )%nr
            !
!            IF( descla( is )%la_myr == descla( is )%la_myc ) THEN
               DO i = 1, nr
                  psihpsi(i,i,is) = 1.0d0
               END DO
!            END IF
         END IF
      end do
    RETURN
  END SUBROUTINE h_matrix_init



  SUBROUTINE ensemble_initval &
    ( occupations_ , n_inner_ , fermi_energy_ ,niter_cold_restart_, lambda_cold_, rotmass_ , occmass_ , rotation_damping_ , &
      occupation_damping_ , occupation_dynamics_ , rotation_dynamics_ ,  degauss_ , smearing_) 
            
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: occupations_
    CHARACTER(LEN=*), INTENT(IN) :: rotation_dynamics_
    CHARACTER(LEN=*), INTENT(IN) :: occupation_dynamics_
    CHARACTER(LEN=*), INTENT(IN) :: smearing_
    INTEGER, INTENT(IN) :: n_inner_
    REAL(DP), INTENT(IN) :: fermi_energy_ , rotmass_ , occmass_ , rotation_damping_
    REAL(DP), INTENT(IN) :: occupation_damping_ , degauss_
    INTEGER, INTENT(in) :: niter_cold_restart_
    REAL(DP), INTENT(in) :: lambda_cold_

      SELECT CASE ( TRIM( occupations_ ) )
          !
      CASE ('bogus')
          !
      CASE ('from_input')
          !
      CASE ('fixed')
          !
      CASE ('grand-canonical','g-c','gc')
          !
          tens    =.true.
          tgrand  =.true.
          CALL errore(' ensemble_initval ','grand-canonical not yet implemented ', 1 )
          !
      CASE ('ensemble','ensemble-dft','edft')
          !
          tens    =.true.
          ninner  = n_inner_
          etemp   = degauss_
          ef      = fermi_energy_
          niter_cold_restart = niter_cold_restart_
          lambda_cold = lambda_cold_

          SELECT CASE ( TRIM( smearing_ ) )
            CASE ( 'gaussian','g' )
              ismear = 1
            CASE ( 'fermi-dirac','f-d', 'fd' )
              ismear = 2
            CASE ( 'hermite-delta','h-d','hd' )
              ismear = 3
            CASE ( 'gaussian-splines','g-s','gs' )
              ismear = 4
            CASE ( 'cold-smearing','c-s','cs','cs1' )
              ismear = 5
            CASE ( 'marzari-vanderbilt','m-v','mv','cs2' )
              ismear = 6
            CASE ( '0')
              ismear = 0
            CASE ( '-1')
              ismear = -1
            CASE DEFAULT
              CALL errore(' ensemble_initval ',' smearing not implemented', 1 )
          END SELECT
          !
      CASE DEFAULT
          !
          CALL errore(' ensemble_initval ',' occupation method not implemented', 1 )
          !
      END SELECT
     
      IF(tens) CALL ensemble_dft_info()

    RETURN
  END SUBROUTINE ensemble_initval


  SUBROUTINE ensemble_dft_info()
    USE io_global, ONLY: stdout
      write(stdout,250) tens
      write(stdout,252) tgrand
250   format (4x,'  ensemble-DFT calculation     =',l5)
252   format (4x,'  grand-canonical calculation  =',l5)

      if(tens) then
         write (stdout,251) ninner,etemp,ismear,ef                         
      endif
251   format (/4x,'====================================='                          &
     &        /4x,'|      ensemble-DFT parameters      |'                          &
     &        /4x,'====================================='                          &
     &        /4x,'| ninner       =',i10,'          |'                             &
     &        /4x,'| etemp        =',f10.5,' a.u.     |'                           &
     &        /4x,'| ismear       =',i10,'          |'                             &
     &        /4x,'| fermi energy =',f10.5,' a.u.     |'                           &
     &        /4x,'=====================================')

      if(tens.and. ismear /= 2) then
         write(stdout,*) 'Full inner-cycle every: ', niter_cold_restart, ' Iterations'
         write(stdout, *) 'With step :', lambda_cold
      endif

    RETURN
  END SUBROUTINE ensemble_dft_info


  SUBROUTINE allocate_ensemble_dft( nhsa, n, ngw, nudx, nspin, nx, nrxxs, nat,&
                                    descla )
    USE descriptors
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nhsa, n, ngw, nudx, nspin, nx, nrxxs, nat
    TYPE(la_descriptor), INTENT(IN) :: descla( nspin )
    
    nrcx = MAXVAL (descla(:)%nrcx )
    nrlx = MAXVAL (descla(:)%nrlx )
    
    allocate(c0diag(ngw,nx))
    allocate(z0t(nrlx,nudx,nspin))
    allocate(becdiag(nhsa,n))
    allocate(e0(nx))
    allocate(fmat0(nrlx,nudx,nspin))
    allocate(psihpsi(nrcx,nrcx,nspin))
    RETURN
  END SUBROUTINE allocate_ensemble_dft



  SUBROUTINE deallocate_ensemble_dft( )
    IMPLICIT NONE
    IF( ALLOCATED( c0diag ) )  deallocate(c0diag )
    IF( ALLOCATED( z0t ) )  deallocate(z0t )
    IF( ALLOCATED( becdiag ) )  deallocate(becdiag )
    IF( ALLOCATED( e0 ) )  deallocate(e0 )
    IF( ALLOCATED( fmat0 ) )  deallocate(fmat0 )
    IF( ALLOCATED( psihpsi ) ) deallocate(psihpsi)
   RETURN
  END SUBROUTINE deallocate_ensemble_dft
  

END MODULE ensemble_dft
