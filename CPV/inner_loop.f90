!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!====================================================================
   SUBROUTINE inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter, vpot )
!====================================================================
      !
      ! minimizes the total free energy with respect to the
      ! occupation matrix f_ij at fixed Kohn-Sham orbitals
      ! Cf. Marzari, Vanderbilt, Payne PRL 79, 1337 (1997)
      !

      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, taurdr, &
                                tprnfor, ndr, ndw, nbeg, nomore, &
                                tsde, tortho, tnosee, tnosep, trane, &
                                tranp, tsdp, tcp, tcap, ampre, &
                                amprp, tnoseh
      USE atom,           ONLY: nlcc
      USE core,           ONLY: nlcc_any
      USE energies,       ONLY: eht, epseu, exc, etot, eself, enl, &
                                ekin, atot, entropy, egrand
      USE electrons_base, ONLY: f, nspin, nel, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE ensemble_dft,   ONLY: tens, tgrand, ninner, ismear, etemp, &
                                ef, z0, c0diag, becdiag, &
                                fmat0,  &
                                e0,   &
                                 compute_entropy2, &
                                compute_entropy_der, compute_entropy
      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: ng0 => gstart
      USE cvan,           ONLY: nvb, ish
      USE ions_base,      ONLY: na, nat, pmass, nax, nsp, rcmax
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE cell_base,      ONLY: ainv, a1, a2, a3
      USE cell_base,      ONLY: omega, alat
      USE cell_base,      ONLY: h, hold, deth, wmass, tpiba2
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx, nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, &
                          ONLY: nnrb => nnrbx, nr1b, nr2b, nr3b
      USE local_pseudo,   ONLY: vps, rhops
      USE io_global,      ONLY: io_global_start, stdout, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm
      USE dener
      USE derho
      USE cdvan
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE cg_module,      ONLY:  ene_ok, &
                                enever
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast
      use cp_interfaces,  only: rhoofr, dforce

      !
      IMPLICIT NONE

      ! declares local variables and counters
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(DP)            :: eigr( ngw, nat )
      COMPLEX(DP)            :: c0( ngw, n )
      REAL(DP)               :: bec( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (DP)           :: eigrb( ngb, nat )
      REAL(DP)               :: rhor( nnr, nspin )
      REAL(DP)               :: vpot( nnr, nspin )
      COMPLEX(DP)            :: rhog( ngm, nspin )
      REAL(DP)               :: rhos( nnrsx, nspin )
      REAL(DP)               :: rhoc( nnr )
      COMPLEX(DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(DP)            :: sfac( ngs, nsp )
      INTEGER                :: i
      INTEGER                :: j
      INTEGER                :: ig
      INTEGER                :: k
      INTEGER                :: is
      INTEGER                :: ia
      INTEGER                :: iv
      INTEGER                :: jv
      INTEGER                :: inl
      INTEGER                :: jnl
      REAL(DP)               :: enb
      REAL(DP)               :: enbi
      COMPLEX(DP)            :: c2( ngw )
      COMPLEX(DP)            :: c3( ngw )
      REAL(DP)               :: gamma
      REAL(DP)               :: entmp
      REAL(DP)               :: sta
      INTEGER                :: npt
      REAL(DP)               :: deltax
      REAL(DP)               :: deltaxmin
      REAL(DP)               :: xinit
      INTEGER                :: ii,jj,kk,nrot
      REAL(kind=DP), ALLOCATABLE :: ztmp(:,:), dval(:), ex(:), dx(:)
      REAL(kind=DP), ALLOCATABLE :: fion2(:,:),c0hc0(:,:,:), z1(:,:,:), zx(:,:,:), zxt(:,:,:), zaux(:,:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: h0c0(:,:)
      REAL(kind=DP), ALLOCATABLE :: f0(:),f1(:),fx(:), faux(:), fmat1(:,:,:),fmatx(:,:,:)
      REAL(kind=DP), ALLOCATABLE :: dfmat(:,:,:), epsi0(:,:,:)
      REAL(kind=DP) :: atot0,atot1,atotmin,etot0,etot1, ef1, enocc
      REAL(kind=DP) :: dadx1,dedx1,dentdx1,eqa,eqb,eqc, etot2, entropy2
      REAL(kind=DP) :: f2,x,xmin
      INTEGER ::  niter,nss,istart,il
   
      CALL start_clock( 'inner_loop' )
      ! initializes variables
      allocate(fion2(3,nat))
      allocate(c0hc0(nudx,nudx,nspin))
      allocate(h0c0(ngw,nx))
      allocate(z1(nudx,nudx,nspin),zx(nudx,nudx,nspin),zxt(nudx,nudx,nspin),zaux(nudx,nudx,nspin))
      allocate(dval(nx),ex(nx),dx(nx),f0(nx),f1(nx),fx(nx),faux(nx))
      allocate(fmat1(nudx,nudx,nspin),fmatx(nudx,nudx,nspin),dfmat(nudx,nudx,nspin))
      allocate(epsi0(nudx,nudx,nspin))

      fion2( :, : )= 0.D0
      npt=10
      deltaxmin=1.D-8
      ! calculates the initial free energy if necessary
      IF( .not. ene_ok ) THEN

        ! calculates the overlaps bec between the wavefunctions c0
        ! and the beta functions
        CALL calbec( 1, nsp, eigr, c0, bec )
 
        ! rotates the wavefunctions c0 and the overlaps bec
        ! (the occupation matrix f_ij becomes diagonal f_i)      
        CALL rotate( z0, c0(:,:), bec, c0diag, becdiag, firstiter)
  
        ! calculates the electronic charge density
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 )
        IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
        vpot = rhor

        ! calculates the SCF potential, the total energy
        ! and the ionic forces
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
        
        ! calculates the entropy
        CALL compute_entropy2( entropy, f, n, nspin )

      END IF
      ene_ok=.FALSE.
      atot0=etot+entropy
      etot0=etot

      ! calculates the occupation matrix 
      ! fmat_ij = \sum_k z_ki * f_k * z_kj
      CALL calcmt( f, z0, fmat0,firstiter )
   
      ! calculateas the energy contribution associated with 
      ! the augmentation charges and the 
      ! corresponding contribution to the ionic force
      CALL newd( vpot, irb, eigrb, rhovan, fion2 )
      CALL prefor( eigr, betae ) ! ATTENZIONE
  
      ! iterates on niter
      INNERLOOP : DO niter= 1, ninner

        ! operates the Hamiltonian on the wavefunction c0
        h0c0( :, : )= 0.D0
          DO i= 1, n, 2                      
            CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
          END DO
        
        ! calculates the Hamiltonian matrix in the basis {c0}           
        c0hc0(:,:,:)=0.d0
        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              DO k= 1, nss
                DO ig= 1, ngw
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    - 2.0d0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
                    * h0c0( ig, i+istart-1 ) )
                END DO
                IF( ng0 .eq. 2 ) THEN
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    + DBLE( CONJG( c0( 1, k+istart-1 ) ) &
                    * h0c0( 1, i+istart-1 ) )
                END IF
              END DO
            END DO
            CALL mp_sum( c0hc0( 1:nss, 1:nss, is ), intra_image_comm )
        END DO
 
        DO is= 1, nspin
          nss= nupdwn( is )
          epsi0( 1:nss, 1:nss, is )= c0hc0( 1:nss, 1:nss, is ) ! ATTENZIONE
        END DO
          
        ! diagonalizes the Hamiltonian matrix
        !e0( : )= 0.D0 This is not set to 0 anymore ecause of blockocc
        DO is= 1, nspin
          istart= iupdwn( is )
          nss= nupdwn( is ) 
            IF( ionode ) THEN
              CALL ddiag( nss, nss, epsi0(1,1,is), dval(1), &
                          z1(1,1,is), 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( z1(:,:,is), ionode_id, intra_image_comm )
            DO i= 1, nss
              e0( i+istart-1 )= dval( i )
            END DO
        END DO
  

        ! calculates the occupations and the fermi energy at the 
        ! end of the search direction
        CALL efermi( nelt, n, etemp, 1, f1, ef1, e0, enocc, ismear, & 
                     nspin )

        ! fmat1_ij = \sum_k z_ik * f_k * z_jk
        CALL calcm( f1, z1, fmat1, firstiter)
             
        ! calculates of dfmat_ij
        ! ( dfmat defines the search direction in occupation space)
        DO is= 1, nspin
          nss= nupdwn( is )
            dfmat( 1:nss, 1:nss, is )= - fmat0( 1:nss, 1:nss, is ) &
                                       + fmat1( 1:nss, 1:nss, is )
        END DO
            
        ! 
        f0( 1:n )= f( 1:n )
       
        ! calculates fmatx= fmat0 + x* dfmat
        ! (here x=1)      
        x=1.D0
        DO is= 1, nspin
          nss= nupdwn( is )
            fmatx( 1:nss, 1:nss, is)= fmat0( 1:nss, 1:nss, is) &
                                    + x * dfmat( 1:nss, 1:nss, is )
        END DO
                      
        ! diagonalizes fmatx 
        fx( : ) = 0.0d0
        DO is=1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            IF( ionode ) THEN
              CALL ddiag( nss, nss, fmatx(1,1,is), dval(1), &
                          zaux(1,1,is), 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( zaux(:,:,is), ionode_id, intra_image_comm )
            DO i= 1, nss
              faux( i+istart-1 )= dval( i )
            END DO
            DO i= 1, nss
              fx( i+istart-1 )= faux( nss-i+istart )
              DO j=1, nss
                zx( i, j, is )= zaux( i, nss-j+1, is )
              END DO
            END DO
            DO i= 1, nss
              DO k= 1, nss
                zxt( k, i, is )= zx( i, k, is )
              END DO
            END DO
        END DO
           
 
        ! updates f
        f( 1:n )= fx( 1:n )

        ! re-calculates fmatx
        CALL calcmt( f, zxt, fmatx, firstiter)
              
        ! calculates the entropy and its derivatives with respect
        ! to each occupation at x
        CALL compute_entropy2( entropy, fx, n, nspin )
        CALL compute_entropy_der( ex, fx, n, nspin )
  
        ! calculates the free energy at x    
        CALL rotate( zxt, c0(:,:), bec, c0diag, becdiag, firstiter)
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 ) 
        IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
        vpot = rhor
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( vpot, irb, eigrb, rhovan, fion2 )
        CALL prefor( eigr, betae )
        atot1=etot+entropy
        etot1=etot
  
        ! calculates the Hamiltonian matrix
        h0c0( :, : )= 0.D0
          DO i= 1, n, 2
            CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
          END DO
        
        c0hc0(:,:,:)=0.d0

        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              DO k= 1, nss
                DO ig= 1, ngw
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    - 2.0d0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
                    * h0c0( ig, i+istart-1 ) )
                END DO
                IF( ng0 .eq. 2 ) THEN
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    + DBLE( CONJG( c0( 1, k+istart-1 ) ) &
                    * h0c0( 1, i+istart-1 ) )
                END IF
              END DO
            END DO
            CALL mp_sum( c0hc0( 1:nss, 1:nss, is ), intra_image_comm )
      ENDDO


        DO is= 1, nspin
          nss= nupdwn( is )
          epsi0( 1:nss, 1:nss, is )= c0hc0( 1:nss, 1:nss, is )
        END DO
                   
        !     calculates 
        !     (1) the energy derivative 
        !         dE / dx (1) = dE / df_ji * (f1_ji - f0_ji) 
        !     (2) the entropy derivative
        !         d Tr S(f) / df_ij = S'(f)_ji
        !         ( d( -T S )/dx is calculated as 
        !           (ex)_j [(zt)_ji (dfmat)_ik (zt)_jk] 
        !           instead of as
        !           [(zt)_jk (ex)_j (zt)_ji] (dfmat)_ik )
        !     (3) the free energy derivative
        dedx1= 0.D0
        dentdx1= 0.D0
        DO is= 1,nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              dx( i+istart-1 )= 0.D0
              DO k= 1, nss
                DO j= 1, nss
                  dx( i+istart-1 )= dx( i+istart-1 ) &
                                  - zxt(i,k,is) * fmat0(k,j,is) &
                                  * zxt(i,j,is)
                END DO
              END DO
              dx( i+istart-1 )= dx( i+istart-1 ) + fx( i+istart-1 )
            END DO
        END DO
        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              dentdx1= dentdx1 - etemp * dx( i+istart-1 ) &
                       * ex(i+istart-1)
              DO k= 1, nss
                dedx1= dedx1 + dfmat( i, k, is ) * epsi0( k, i, is )
              END DO
            END DO
        END DO
        dadx1 = dedx1 + dentdx1
                   
        ! performs the line minimization
        ! (a) the energy contribution is approximated 
        !     by a second order polynomial
        ! (b) the entropic term is approcimated by a function 
        !     of the form \sum_i s(a_i*x**2+b_i*x+c_i)
        !    (where s(f)=-f*ln(f)-(1-f)*ln(1-f) ).
        !    The coefficients a_i, b_i and c_i are calculated
        !    by first-order perturbation
        eqc= etot0
        eqa= dedx1 - etot1 + etot0
        eqb= etot1 - etot0 - eqa
        atotmin= atot0
        xmin= 0.D0
        xinit=0.D0
        deltax= 1.D0 / DBLE( npt )
        DO WHILE ( deltax .gt. deltaxmin )
          SAMPLING : DO il= 0, npt
            x= xinit + deltax * DBLE( il )
            IF( x .gt. 1.D0 ) EXIT SAMPLING
            entropy2= 0.D0
            DO is=1,nspin
              nss= nupdwn( is )
              istart= iupdwn( is )
              DO i= 1, nss
                f2= fx( i+istart-1 ) + ( x-1 ) * dx( i+istart-1 ) &
                  + ( - fx( i+istart-1 ) + dx( i+istart-1 ) + &
                  f0( i+istart-1 ) ) * ( x-1 )**2
                CALL compute_entropy( entmp, f2, nspin )
                entropy2 = entropy2 + entmp
              END DO
            END DO
            etot2= eqa * x ** 2 + eqb * x + eqc
            IF ( ( etot2 + entropy2 ) .lt. atotmin ) THEN
              xmin= x
              atotmin= etot2+ entropy2
            END IF
          END DO SAMPLING
          xinit= MAX( 0.D0, xmin - deltax )
          deltax= 2.D0 * deltax / DBLE( npt )
        END DO

        IF( ionode ) THEN
          WRITE(37,'(a5,3f15.10)') &
                             'XMIN', xmin, atotmin, atotmin-atot0
          IF ( atotmin-atot0 .gt. 0.D0 ) &
          WRITE(37,*) "INNER LOOP, WARNING : increasing free energy"
        END IF
       
        !
        IF ( xmin .eq. 1.0D0 ) GOTO 300
  
        ! calculates the occupation matrix at xmin
        DO is= 1, nspin
          nss= nupdwn( is )
            DO i= 1, nss
              DO j= 1, nss
                fmatx( i, j, is )= fmat0( i, j, is ) &
                                   + xmin * dfmat( i, j, is )
              END DO
            END DO
        END DO
               
        ! diagonalizes the occupation matrix at xmin
        fx( : )= 0.D0
        DO is= 1, nspin
          nss= nupdwn( is ) 
          istart= iupdwn( is )
            IF(ionode) CALL ddiag( nss, nss, fmatx(1,1,is), &
                                 dval(1), zaux(1,1,is), 1 )  
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( zaux(:,:,is), ionode_id, intra_image_comm )
            DO i= 1, n
              faux( i+istart-1 )= dval( i )
            END DO
            DO i= 1, nss
              fx( i+istart-1 )= faux( nss-i+istart )
              DO j= 1, nss
                zx( i, j, is )= zaux( i, nss-j+1, is )
              END DO
            END DO
        END DO
    
        ! updates f
        f( 1:n )= fx( 1:n )
  
 300    CONTINUE
         
        ! updates z0
        DO is= 1, nspin
          nss= nupdwn( is )
          DO i= 1, nss
            DO k= 1, nss
              z0( k, i, is )= zx( i, k, is )
            END DO
          END DO
        END DO
        
        ! calculates the total free energy
        CALL calcmt( f, z0, fmat0, firstiter)
        CALL rotate( z0, c0(:,:), bec, c0diag, becdiag, firstiter )
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, &
                     rhovan, rhor, rhog, rhos, enl, denl, ekin, dekin6 ) 
        IF(nlcc_any) CALL set_cc(irb,eigrb,rhoc)
        vpot = rhor
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( vpot, irb, eigrb, rhovan, fion2 )
        CALL compute_entropy2( entropy, f, n, nspin )
        CALL prefor( eigr, betae )
        ene_ok= .TRUE.
        atotmin = etot + entropy      
        IF (ionode) write(37,'(a3,i2,2f15.10)') 'CI',niter,atot0,atotmin
        atot0=atotmin
        atot=atotmin
        etot0=etot
        enever=etot
        if(xmin==0.d0) exit 
     END DO INNERLOOP

     deallocate(fion2,c0hc0,h0c0,z1)
     deallocate(zx,zxt,zaux,dval,ex,dx)
     deallocate(f0,f1,fx,faux,fmat1,fmatx)
     deallocate(dfmat,epsi0)
     CALL stop_clock( 'inner_loop' )
     return
!====================================================================      
   END SUBROUTINE inner_loop
!====================================================================
