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
   SUBROUTINE inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter)
!====================================================================
      !
      ! minimizes the total free energy
      ! using cold smearing,
      !
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
                                fmat0,  &
                                e0, psihpsi, compute_entropy2, &
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
      use cp_interfaces,  only: rhoofr, dforce

      !
      IMPLICIT NONE

!input variables
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(kind=DP)            :: eigr( ngw, nat )
      COMPLEX(kind=DP)            :: c0( ngw, n )
      REAL(kind=DP)               :: bec( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (kind=DP)           :: eigrb( ngb, nat )
      REAL(kind=DP)               :: rhor( nnr, nspin )
      COMPLEX(kind=DP)            :: rhog( ngm, nspin )
      REAL(kind=DP)               :: rhos( nnrsx, nspin )
      REAL(kind=DP)               :: rhoc( nnr )
      COMPLEX(kind=DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(kind=DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(kind=DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(kind=DP)            :: sfac( ngs, nsp )
  

!local variables
      REAL(kind=DP) :: atot0, atot1, atotl, atotmin
      REAL(kind=DP), ALLOCATABLE :: fion2(:,:), c0hc0(:,:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: h0c0(:,:)
      INTEGER :: niter
      INTEGER :: i,k, is, nss, istart, ig
      REAL(kind=DP) :: lambda, lambdap
      REAL(kind=DP), ALLOCATABLE :: epsi0(:,:,:), dval(:), zaux(:,:,:)

      CALL start_clock( 'inner_loop')

      allocate(fion2(3,nat))
      allocate(c0hc0(nudx,nudx,nspin))
      allocate(h0c0(ngw,nx))
      allocate(epsi0(nudx,nudx,nspin))
      allocate(dval(nx))
      allocate(zaux(nudx,nudx,nspin))


      lambdap=0.3d0!small step for free-energy calculation
      
 


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
       !entropy value already  been calculated

      END IF
     
      atot0=etot+entropy
      
    
!starts inner loop
      do niter=1,ninner
!calculates c0hc0, which defines the search line (1-\labda)* psihpsi+\labda*c0hc0

      ! calculateas the energy contribution associated with 
      ! the augmentation charges and the 
      ! corresponding contribution to the ionic force
       
         CALL newd( rhor, irb, eigrb, rhovan, fion2 )

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

       

!calculates free energy at lamda=1.
         CALL inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter,psihpsi,c0hc0,1.d0,atot1)
!calculates free energy at lamda=lambdap
       
         CALL inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter,psihpsi,c0hc0,lambdap,atotl)
!find minimum point lambda
        
         CALL three_point_min(atot0,atotl,atot1,lambdap,lambda,atotmin)
        


!calculates free energy and rho at lambda
        

!calculates the new matrix psihpsi
         
         DO is= 1, nspin
            nss= nupdwn( is )
            psihpsi(:,:,is)=(1.d0-lambda)*psihpsi(:,:,is)+lambda*c0hc0(:,:,is)
         END DO
!diagonalize and calculates energies
         
         DO is= 1, nspin
            nss= nupdwn( is )
            epsi0( 1:nss, 1:nss, is )= psihpsi( 1:nss, 1:nss, is )
         END DO
          
     
         
        ! diagonalizes the Hamiltonian matrix
         e0( : )= 0.D0 
         DO is= 1, nspin
            istart= iupdwn( is )
            nss= nupdwn( is )
!only one processor diagonalizes TO BE PARALLELIZED!!! 
            IF( ionode ) THEN
               CALL ddiag( nss, nss, epsi0(:,:,is), dval, &
                          zaux(:,:,is), 1 )
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL mp_bcast( zaux(:,:,is), ionode_id, intra_image_comm )
            DO i= 1, nss
               e0( i+istart-1 )= dval( i )
            END DO
         END DO
  
!NB zaux is transposed

!calculates fro e0 the new occupation and the entropy        

       

         CALL efermi( nelt, n, etemp, 1, f, ef, e0, entropy, ismear, & 
                     nspin )

!alculates z0
         do is=1,nspin
            nss=nupdwn(is)
            do i=1,nss
               do k=1,nss
                  z0(k,i,is)=zaux(i,k,is)
               enddo
            enddo
         enddo


!calculates new charge and new energy
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
       !entropy value already  been calculated

         write(37,*) niter
         write(37,*) atot0,atotl,atot1
         write(37,*) lambda,atotmin,etot+entropy
       
         atotl=atot0
         atot0=etot+entropy

   
         if(lambda==0.d0) exit 


         
      enddo


      atot=atot0

!ATTENZIONE
!the following is of capital importance
      ene_ok= .FALSE.
!but it would be better to avoid it


      DEALLOCATE(fion2)
      DEALLOCATE(c0hc0)
      DEALLOCATE(h0c0)
      DEALLOCATE(epsi0)
      deallocate(dval)
      deallocate(zaux)


      CALL stop_clock( 'inner_loop' )
      return
!====================================================================      
    END SUBROUTINE inner_loop_cold
!====================================================================





   SUBROUTINE inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter,c0hc0,c1hc1,lambda,free_energy)
    
!this subroutine for the energy matrix (1-lambda)c0hc0+labda*c1hc1
!calculates the corresponding free energy


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
                                 c0diag, becdiag
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
      use cp_interfaces,  only: rhoofr, dforce

      !
      IMPLICIT NONE

!input variables
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(kind=DP)            :: eigr( ngw, nat )
      COMPLEX(kind=DP)            :: c0( ngw, n )
      REAL(kind=DP)               :: bec( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (kind=DP)           :: eigrb( ngb, nat )
      REAL(kind=DP)               :: rhor( nnr, nspin )
      COMPLEX(kind=DP)            :: rhog( ngm, nspin )
      REAL(kind=DP)               :: rhos( nnrsx, nspin )
      REAL(kind=DP)               :: rhoc( nnr )
      COMPLEX(kind=DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(kind=DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(kind=DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(kind=DP)            :: sfac( ngs, nsp )
  
      REAL(kind=DP), INTENT(in)   :: c0hc0(nudx,nudx,nspin),c1hc1(nudx,nudx,nspin)
      REAL(kind=DP), INTENT(in)   :: lambda
      REAL(kind=DP), INTENT(out)  :: free_energy


!local variables
      REAL(kind=DP), ALLOCATABLE  :: clhcl(:,:,:), fion2(:,:)
      INTEGER :: i,k, is, nss, istart, ig
      REAL(kind=DP), ALLOCATABLE :: eaux(:),faux(:),zaux(:,:,:), dval(:),zauxt(:,:,:)
      REAL(kind=DP) :: entropy_aux, ef_aux

      
      allocate(clhcl(nudx,nudx,nspin))
      allocate(eaux(nx))
      allocate(faux(nx))
      allocate(zaux(nudx,nudx,nspin))
      allocate(zauxt(nudx,nudx,nspin))
      allocate(dval(nx))
      allocate(fion2(3,nat))


!calculates the matrix clhcl
         
      DO is= 1, nspin
         nss= nupdwn( is )
         clhcl(:,:,is)=(1.d0-lambda)*c0hc0(:,:,is)+lambda*c1hc1(:,:,is)
      END DO

      
! diagonalizes the Hamiltonian matrix
      DO is= 1, nspin
         istart= iupdwn( is )
         nss= nupdwn( is )
!only one processor diagonalizes TO BE PARALLELIZED!!! 
         IF( ionode ) THEN
            CALL ddiag( nss, nss, clhcl(:,:,is), dval, &
                 zaux(:,:,is), 1 )
         END IF
         CALL mp_bcast( dval, ionode_id, intra_image_comm )
         CALL mp_bcast( zaux(:,:,is), ionode_id, intra_image_comm )
         DO i= 1, nss
            eaux( i+istart-1 )= dval( i )
         END DO
      END DO
  
!NB zaux is transposed

!calculates fro e0 the new occupation and the entropy        

!rhoofr uses array f
!  we save it

      faux(:)=f(:)

      CALL efermi( nelt, n, etemp, 1, f, ef_aux, eaux, entropy_aux, ismear, & 
           nspin )

!calculates transposed matrix
         do is=1,nspin
            nss=nupdwn(is)
            do i=1,nss
               do k=1,nss
                  zauxt(k,i,is)=zaux(i,k,is)
               enddo
            enddo
         enddo


!calculates new charge and new energy
! rotates the wavefunctions c0 and the overlaps bec
! (the occupation matrix f_ij becomes diagonal f_i)      

         CALL rotate( zauxt, c0(:,:), bec, c0diag, becdiag, firstiter)
  


        ! calculates the electronic charge density
         CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
        ! calculates the SCF potential, the total energy
        ! and the ionic forces
         CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
       !entropy value already  been calculated

         free_energy=etot+entropy_aux
         
         f(:)=faux(:)


      deallocate(clhcl)
      deallocate(dval)
      deallocate(zaux)
      deallocate(faux)
      deallocate(eaux)
      deallocate(zauxt)
      deallocate(fion2)

      return


    END SUBROUTINE inner_loop_lambda


    SUBROUTINE three_point_min(y0,yl,y1,l,lambda,amin)
!calculates the estimate for the minimum 

      USE kinds,  ONLY : DP


      implicit none

      REAL(kind=DP), INTENT(in) :: y0,yl,y1, l
      REAL(kind=DP), INTENT(out) :: lambda, amin


      REAL(kind=DP) :: a,b,c, x_min, y_min

! factors for f(x)=ax**2+b*x+c
      c=y0
      b=(yl-y0-y1*l**2.d0+y0*l**2.d0)/(l-l**2.d0)
      a=y1-y0-b
   
      
      x_min=-b/(2.d0*a)
      if( x_min <= 1.d0 .and. x_min >= 0.d0) then
         y_min=a*x_min**2.d0+b*x_min+c
         if(y_min <= y0 .and. y_min <= y1) then
            lambda=x_min
            amin=y_min
         else
            if(y0 < y1) then
                lambda=0.d0
                amin=y0
             else
                lambda=1.d0
                amin=y1
             endif
          endif
       else
          if(y0 < y1) then
             lambda=0.d0
             amin=y0
          else
             lambda=1.d0
             amin=y1
          endif
       endif


       return

     END SUBROUTINE three_point_min
