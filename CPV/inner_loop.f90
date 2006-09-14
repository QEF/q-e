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
                          sfac, c0, bec, firstiter)
!====================================================================
      !
      ! minimizes the total free energy with respect to the
      ! occupation matrix f_ij at fixed Kohn-Sham orbitals
      ! Cf. Marzari, Vanderbilt, Payne PRL 79, 1337 (1997)
      !

      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, tvlocw, taurdr, &
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
                                fmat0, fion2, atot0, etot0, h0c0, &
                                c0hc0, epsi0, e0, dval, z1, f1, &
                                dfmat, fmat1, ef1, enocc, f0, fmatx, &
                                fx, zaux, zx, ex, zxt, atot1, etot1, &
                                dedx1, dentdx1, dx, dadx1, faux, eqa, &
                                eqb, eqc, atotmin, xmin, entropy2, &
                                f2, etot2, compute_entropy2, &
                                compute_entropy_der, compute_entropy, &
                                l_blockocc, n_blockocc
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
      USE stre
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE cg_module,      ONLY: ltresh, itercg, etotnew, etotold, &
                                tcutoff, restartcg, passof, passov, &
                                passop, ene_ok, numok, maxiter, &
                                enever, conv_thr, ene0, &
                                esse, essenew, dene0, spasso, ene1, &
                                passo, iter3, enesti, ninner_ef
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast
      use cp_interfaces,  only: rhoofr

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
      INTEGER                :: il
      INTEGER                :: inl
      INTEGER                :: jnl
      INTEGER                :: niter
      INTEGER                :: istart
      INTEGER                :: nss
      REAL(DP)               :: enb
      REAL(DP)               :: enbi
      REAL(DP)               :: x
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
      REAL(kind=DP), ALLOCATABLE :: ztmp(:,:)
 
      CALL start_clock( 'inner_loop' )
      ! initializes variables
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
  
        ! calculates the SCF potential, the total energy
        ! and the ionic forces
        CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
        
        ! calculates the entropy
        CALL compute_entropy2( entropy, f, n, nspin )

      END IF
      ene_ok=.FALSE.
      atot0=etot+entropy
      write(37,*) 'energy entropy 0', etot, entropy!ATTENZIONE
      etot0=etot

      ! calculates the occupation matrix 
      ! fmat_ij = \sum_k z_ki * f_k * z_kj
      CALL calcmt( f, z0, fmat0,firstiter )
   
      ! calculateas the energy contribution associated with 
      ! the augmentation charges and the 
      ! corresponding contribution to the ionic force
      CALL newd( rhor, irb, eigrb, rhovan, fion2 )
      CALL prefor( eigr, betae ) ! ATTENZIONE
  
      ! iterates on niter
      INNERLOOP : DO niter= 1, ninner

        ! operates the Hamiltonian on the wavefunction c0
        h0c0( :, : )= 0.D0
        if(firstiter .or. .not. l_blockocc) then
          DO i= 1, n, 2                      
            CALL dforce( bec, betae, i, c0(1,i), &
                         c0(1,i+1), h0c0(1,i), h0c0(1,i+1), &
                         rhos, ispin, f, n, nspin )
          END DO
        else
          do is=1,nspin
             nss= nupdwn( is )
            istart= iupdwn( is )
            DO i=istart+n_blockocc(is),istart+nss-1,2
                CALL dforce( bec, betae, i, c0(1,i), &
                         c0(1,i+1), h0c0(1,i), h0c0(1,i+1), &
                         rhos, ispin, f, n, nspin )
            END DO
          enddo
        endif
        
        ! calculates the Hamiltonian matrix in the basis {c0}           
        c0hc0(:,:,:)=0.d0
        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
          if(firstiter .or. .not. l_blockocc) then
            DO i= 1, nss
              DO k= 1, nss
                DO ig= 1, ngw
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    - 2.0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
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
          else
            nrot=nss-n_blockocc(is)
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              DO k= n_blockocc(is)+1,nss
                kk=k-n_blockocc(is)
                DO ig= 1, ngw
                  c0hc0( kk, ii, is )= c0hc0( kk, ii, is ) &
                    - 2.0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
                    * h0c0( ig, i+istart-1 ) )
                END DO
                IF( ng0 .eq. 2 ) THEN
                  c0hc0( kk, ii, is )= c0hc0( kk, ii, is ) &
                    + DBLE( CONJG( c0( 1, k+istart-1 ) ) &
                    * h0c0( 1, i+istart-1 ) )
                END IF
              END DO
            END DO
            CALL mp_sum( c0hc0( 1:nrot, 1:nrot, is ), intra_image_comm )
         endif
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
          if(firstiter .or. .not.l_blockocc) then
            IF( ionode ) THEN
              CALL ddiag( nss, nss, epsi0(1,1,is), dval(1), &
                          z1(1,1,is), 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( z1(:,:,is), ionode_id, intra_image_comm )
            DO i= 1, nss
              e0( i+istart-1 )= dval( i )
            END DO
          else
           nrot=nss-n_blockocc(is)
           allocate(ztmp(nss,nrot))
           IF( ionode ) THEN
              CALL ddiag( nss, nrot, epsi0(1,1,is), dval(1), &
                          ztmp, 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( ztmp(:,:), ionode_id, intra_image_comm )
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              e0( i+istart-1 )= dval( ii )
            END DO
            z1(:,:,is)=0.d0
            do i=1,n_blockocc(is)
              z1(i,i,is)=1.d0
            enddo
            do i=n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              do j=n_blockocc(is)+1,nss
                 jj=j-n_blockocc(is)
                 z1(i,j,is)=ztmp(ii,jj)
              enddo
            enddo
            deallocate(ztmp)
          endif
        END DO
  

        ! calculates the occupations and the fermi energy at the 
        ! end of the search direction
        CALL efermi( nelt, n, etemp, 1, f1, ef1, e0, enocc, ismear, & 
                     nspin )
        if(l_blockocc) then
!check occupancy of blocked states
           do is=1,nspin
             istart= iupdwn( is )
             nss= nupdwn( is )
             do i=1,n_blockocc(is)
               write(37,*) 'Occupancy: ',i,is,f1(i+istart-1)!ATTENZIONE
             enddo
           enddo
        endif


        ! fmat1_ij = \sum_k z_ik * f_k * z_jk
        CALL calcm( f1, z1, fmat1, firstiter)
             
        ! calculates of dfmat_ij
        ! ( dfmat defines the search direction in occupation space)
        DO is= 1, nspin
          nss= nupdwn( is )
          if(firstiter .or. .not. l_blockocc) then  
            dfmat( 1:nss, 1:nss, is )= - fmat0( 1:nss, 1:nss, is ) &
                                       + fmat1( 1:nss, 1:nss, is )
          else
            nrot=nss-n_blockocc(is)
            dfmat( 1:nrot, 1:nrot, is )= - fmat0( 1:nrot, 1:nrot, is ) &
                                       + fmat1( 1:nrot, 1:nrot, is )
          endif
        END DO
            
        ! 
        f0( 1:n )= f( 1:n )
       
        ! calculates fmatx= fmat0 + x* dfmat
        ! (here x=1)      
        x=1.D0
        DO is= 1, nspin
          nss= nupdwn( is )
          if(firstiter .or. .not. l_blockocc) then
            fmatx( 1:nss, 1:nss, is)= fmat0( 1:nss, 1:nss, is) &
                                    + x * dfmat( 1:nss, 1:nss, is )
          else
            nrot=nss-n_blockocc(is)
            fmatx( 1:nrot, 1:nrot, is)= fmat0( 1:nrot, 1:nrot, is) &
                                    + x * dfmat( 1:nrot, 1:nrot, is )
          endif
        END DO
                      
        ! diagonalizes fmatx 
        fx( : ) = 0.0d0
        DO is=1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
          if(firstiter .or. .not. l_blockocc) then
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
          else
            nrot=nss-n_blockocc(is)
            allocate(ztmp(nss,nrot))
            IF( ionode ) THEN
              CALL ddiag( nss, nrot, fmatx(1,1,is), dval(1), &
                          ztmp, 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( ztmp(:,:), ionode_id, intra_image_comm )
!now the occupancy are inverted from 0. to 1. ...
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              faux( i+istart-1 )= dval(ii)
            END DO
            zx(:,:,is)=0.d0
            DO i= n_blockocc(is)+1,nss
              ii=i- n_blockocc(is)
              fx( i+istart-1 )= faux( nss-ii+istart )
              DO j=n_blockocc(is)+1,nss
                jj=j- n_blockocc(is)
                zx( i, j, is )= ztmp( ii, nrot-jj+1)
              END DO
            END DO
            do i=1,n_blockocc(is)
               zx(i,i,is)=1.d0
               fx(i+istart-1)=f(i+istart-1)
            enddo
            DO i= 1, nss
              DO k= 1, nss
                zxt( k, i, is )= zx( i, k, is )
              END DO
            END DO
            deallocate(ztmp)
          endif
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
        CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( rhor, irb, eigrb, rhovan, fion2 )
        CALL prefor( eigr, betae )
        atot1=etot+entropy
        write(37,*) 'energy entropy 1', etot, entropy!ATTENZIONE
        etot1=etot
  
        ! calculates the Hamiltonian matrix
        h0c0( :, : )= 0.D0
        if(firstiter .or. .not. l_blockocc) then
          DO i= 1, n, 2
            CALL dforce( bec, betae, i, c0(1,i), &
                         c0(1,i+1), h0c0(1,i), h0c0(1,i+1), &
                         rhos, ispin, f, n, nspin )
          END DO
        else
          do is=1,nspin
            nss= nupdwn( is )
            istart= iupdwn( is )
            DO i=istart+n_blockocc(is),istart+nss-1,2
                CALL dforce( bec, betae, i, c0(1,i), &
                         c0(1,i+1), h0c0(1,i), h0c0(1,i+1), &
                         rhos, ispin, f, n, nspin )
            END DO
          enddo
        endif


        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
          if(firstiter .or. .not. l_blockocc) then
            DO i= 1, nss
              DO k= 1, nss
                DO ig= 1, ngw
                  c0hc0( k, i, is )= c0hc0( k, i, is ) &
                    - 2.0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
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
          else
            nrot=nss-n_blockocc(is)
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              DO k= n_blockocc(is)+1,nss
                kk=k-n_blockocc(is)
                DO ig= 1, ngw
                  c0hc0( kk, ii, is )= c0hc0( kk, ii, is ) &
                    - 2.0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
                    * h0c0( ig, i+istart-1 ) )
                END DO
                IF( ng0 .eq. 2 ) THEN
                  c0hc0( kk, ii, is )= c0hc0( kk, ii, is ) &
                    + DBLE( CONJG( c0( 1, k+istart-1 ) ) &
                    * h0c0( 1, i+istart-1 ) )
                END IF
              END DO
            END DO
            CALL mp_sum( c0hc0( 1:nrot, 1:nrot, is ), intra_image_comm )
          endif
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
          if(firstiter .or. .not.l_blockocc) then
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
          else
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              dx( i+istart-1 )= 0.D0
              DO k= n_blockocc(is)+1,nss
                kk=k-n_blockocc(is)
                DO j= n_blockocc(is)+1,nss
                  jj=j- n_blockocc(is)
                  dx( i+istart-1 )= dx( i+istart-1 ) &
                                  - zxt(i,k,is) * fmat0(kk,jj,is) &
                                  * zxt(i,j,is)
                END DO
              END DO
              dx( i+istart-1 )= dx( i+istart-1 ) + fx( i+istart-1 )
            END DO
         endif
        END DO
        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
          if(firstiter .or. .not.l_blockocc) then
            DO i= 1, nss
              dentdx1= dentdx1 - etemp * dx( i+istart-1 ) &
                       * ex(i+istart-1)
              DO k= 1, nss
                dedx1= dedx1 + dfmat( i, k, is ) * epsi0( k, i, is )
              END DO
            END DO
          else
            DO i= n_blockocc(is)+1,nss
              ii=i- n_blockocc(is)
              dentdx1= dentdx1 - etemp * dx( i+istart-1 ) &
                       * ex(i+istart-1)
              DO k= n_blockocc(is)+1,nss
                kk=k-n_blockocc(is)
                dedx1= dedx1 + dfmat( ii, kk, is ) * epsi0( kk, ii, is )
              END DO
            END DO
          endif
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
          if(firstiter .or. .not.l_blockocc) then
            DO i= 1, nss
              DO j= 1, nss
                fmatx( i, j, is )= fmat0( i, j, is ) &
                                   + xmin * dfmat( i, j, is )
              END DO
            END DO
          else
            nrot=nss-n_blockocc(is)
            DO i= 1, nrot
              DO j= 1, nrot
                fmatx( i, j, is )= fmat0( i, j, is ) &
                                   + xmin * dfmat( i, j, is )
              END DO
            END DO
          endif
        END DO
               
        ! diagonalizes the occupation matrix at xmin
        fx( : )= 0.D0
        DO is= 1, nspin
          nss= nupdwn( is ) 
          istart= iupdwn( is )
          if(firstiter .or. .not. l_blockocc) then
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
          else
            nrot=nss-n_blockocc(is)
            allocate(ztmp(nss,nrot))
            IF(ionode) CALL ddiag( nss, nrot, fmatx(1,1,is), &
                                 dval(1), ztmp, 1 )
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( ztmp(:,:), ionode_id, intra_image_comm )


!now the occupancy are inverted from 0. to 1. ...
            DO i= n_blockocc(is)+1,nss
              ii=i-n_blockocc(is)
              faux( i+istart-1 )= dval(ii)
            END DO
            zx(:,:,is)=0.d0
            DO i= n_blockocc(is)+1,nss
              ii=i- n_blockocc(is)
              fx( i+istart-1 )= faux( nss-ii+istart )
              DO j=n_blockocc(is)+1,nss
                jj=j- n_blockocc(is)
                zx( i, j, is )= ztmp( ii, nrot-jj+1 )
              END DO
            END DO
            do i=1,n_blockocc(is)
               zx(i,i,is)=1.d0
               fx(i+istart-1)=f(i+istart-1)
            enddo
            deallocate(ztmp)
         endif
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
        CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( rhor, irb, eigrb, rhovan, fion2 )
        CALL compute_entropy2( entropy, f, n, nspin )
        CALL prefor( eigr, betae )
        ene_ok= .TRUE.
        atotmin = etot + entropy      
        write(37,*) 'energy entropy 2',etot,entropy!ATTENZIONE
        IF (ionode) write(37,'(a3,i2,2f15.10)') 'CI',niter,atot0,atotmin
        atot0=atotmin
        atot=atotmin
        etot0=etot
        enever=etot
        if(xmin==0.d0) exit 
     END DO INNERLOOP
!if we use the blockocc and firstiteration bring everything into the diagonal representation

     if(firstiter .and. l_blockocc) then
        c0(:,:)=c0diag(:,:)
        bec(:,:)=becdiag(:,:)
        z0(:,:,:)=0.d0
        do is=1,nspin
           nss=nupdwn(is)
           do i=1,nss
              z0(i,i,is)=1.d0
           enddo
        enddo
      endif



     CALL stop_clock( 'inner_loop' )
     return
!====================================================================      
   END SUBROUTINE inner_loop
!====================================================================
