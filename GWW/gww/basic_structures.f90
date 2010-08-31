! P.Umari GW code
! Modified by G. Stenuit
!
MODULE basic_structures
!this module describes the basis structures
!which are obtained from a DFT code
  USE kinds, ONLY : DP

  TYPE wannier_u
!this structure describes the  transformation
!from KS to Wannier states together with the KS eigenenergies
    INTEGER :: nums!number of states
    INTEGER :: nums_occ!number of occupied states
    REAL(kind=DP), DIMENSION(:), POINTER :: ene!KS energies
    REAL(kind=DP), DIMENSION(:), POINTER :: ene_xc!LDA exchange and correlation terms
    REAL(kind=DP), DIMENSION(:), POINTER :: ene_u!LDA Hubbard terms
    REAL(kind=DP), DIMENSION(:), POINTER :: ene_lda_h!LDA exchange and correlation terms
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: umat!INVERSE transformation matrix to wannier: Psi_i=U_{i,j}w_j
  END TYPE wannier_u


  TYPE wannier_P
!this structure described the localized and normalized products of wanniers  w_P
    INTEGER :: numij!number of (unordered)couples of wannier w_i w_j which overlap both with w_P
    INTEGER, DIMENSION(:,:), POINTER :: ij!array for i,j (ij,:)
    REAL(kind=DP),DIMENSION(:), POINTER :: o!overlap <w_P|w_i*w_j>
  END TYPE wannier_P

  TYPE v_pot
!this structure describes the coulomb potential on the base of (orthonormalized)
!products of wanniers
    INTEGER :: numpw!number of states(products)
    REAL(kind=DP),DIMENSION(:,:), POINTER :: vmat!potentail 1/|r-r'|
  END TYPE v_pot

 TYPE q_mat
!this structures describes the set of overlap of othonormalized products
!of wanniers with products of wannier
    INTEGER :: numpw!number of states(orthonormalized products)
!parameters used for  parallelization
    LOGICAL :: is_parallel!if true is a parallel part of the global q_mat matrix
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor
    TYPE(wannier_P), DIMENSION(:), POINTER :: wp!arrays of wannier products descriptors
 END TYPE q_mat


 TYPE ortho_polaw
!this structure describe the orthonormalization matrix
!w^P_i=A_{i,j}\tilde{w^P}_j
!it is put here because it's read from a PW file
    INTEGER :: numpw!number of states (products of wanniers)
    LOGICAL :: inverse!if true, the inverse transform is stored
    REAL(kind=DP), DIMENSION(:,:), POINTER :: on_mat!the transformation
 END TYPE ortho_polaw

 TYPE wp_psi
!this structure describe the product of KS wavefunctions with unorthonormalized
!products of wannier \int dr w^P_i(r) w^P_j(r) Psi_v(r) Psi_v(r)
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    REAL(kind=DP), DIMENSION(:,:,:), POINTER :: wwp!terms
 END TYPE wp_psi

 TYPE wannier_u_prim
!this structure describes the  transformation
!from KS to Wannier states in the manifold C'
    INTEGER :: nums!total number of states
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_prim!number of states in manifold C'
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: umat!INVERSE transformation matrix to wannier: Psi_c'=U_{c',j}w_j
 END TYPE wannier_u_prim

 TYPE v_pot_prim
!this structure describes the terms <\tilde{w}^P'_i|V|\tilde{w}^P_j'>
    INTEGER ::  numpw!number of states (products of wanniers)
    INTEGER ::  numpw_prim!number of states in manifold C' (products of wanniers)
    INTEGER, DIMENSION(:,:), POINTER :: ij!array of dimesion(2,numpw_prim) defining the product of w^C'*w^C
    REAL(kind=DP), DIMENSION(:,:), POINTER :: vmat!coulumbian matrix
    LOGICAL :: is_parallel!if true is a parallel part of the global cprim_prod matrix on polarization basis
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor
 END TYPE v_pot_prim

 TYPE wp_psi_cutoff_index
!this structure contains the indices for the description of terms \int Psi_i(r)\tilde{w}^P_i\tilde{w}^P_jPsi_i(r)dr
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    INTEGER :: numpwpw!number of products of wannier products
    INTEGER, DIMENSION(:,:), POINTER :: index! of dimension (2,numpwpw) indices to wannier products
 END TYPE wp_psi_cutoff_index

 TYPE wp_psi_cutoff_data
!this structure contains the data for the description of terms \int Psi_i(r)\tilde{w}^P_i\tilde{w}^P_jPsi_i(r)dr
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    INTEGER :: numpwpw!number of products of wannier products
    REAL(kind=DP), DIMENSION(:,:), POINTER :: wwp!  of dimension (numpwpw,nums_psi)
 END TYPE wp_psi_cutoff_data

 TYPE  head_epsilon
!this structure contains the data for the descrpition of the head of the dielectric matrix
!calculated accurately with k_points sampling
!it also contains the data for the treatment of the wings
    INTEGER :: n!number of frequency steps
    REAL(kind=DP) :: omega!frequency range
    REAL(kind=DP), DIMENSION(:), POINTER :: freqs!frequency steps 2n+1
    REAL(kind=DP), DIMENSION(:), POINTER :: head!elements G=0,G=0 of the dielectric matrix
    INTEGER :: numpw!number of products of wanniers
    REAL(kind=DP), DIMENSION(:), POINTER :: gzero!G=0 elements of non orthogonal products of wanniers \tilde{w^P}
    REAL(kind=DP), DIMENSION(:,:), POINTER :: wing!contains the terms \Sum_G <G|\tilde{w^P_i}>\epsilon(G,G'=0; iw)
    REAL(kind=DP), DIMENSION(:,:), POINTER :: wing_c!contains the terms \Sum_G <G|\tilde{w^P_i}>\epsilon(G=0,G'; iw)
 END TYPE head_epsilon

 TYPE  cprim_prod
!this structure contains the terms \int Psi_c'(r) Psi_c(r) v(r,r') \tilde{w^P_i}dr dr'
!it can contain also the terms \int Psi_i(r) Psi_v,c(r) v(r,r') \tilde{w^P_i}dr dr'
!it can contain also the terms \int Psi_v(r) Psi_c(r) \tilde{w^P_i} dr
    INTEGER :: cprim!conduction band considered
    INTEGER :: nums!total number of states
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_cond!total number of conduction states
    INTEGER :: numpw!number of products of wanniers
    REAL(kind=DP), DIMENSION(:,:), POINTER :: cpmat!product terms
    INTEGER  :: lda!leading dimension of cpmat important for parallel execution
    LOGICAL :: is_parallel!if true is a parallel part of the global cprim_prod matrix on polarization basis
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor

 END TYPE cprim_prod


 TYPE upper_states
!this structure contains the data for the reduced upper states
    INTEGER :: nums!total number of REGULAR states
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_reduced!number of reduced states
    INTEGER :: nums_tot!number of TOTAL states
    REAL(kind=DP), DIMENSION(:), POINTER :: ene!KS energies of reduced states
 END TYPE upper_states


!these routines deallocate the allocates structures
      INTERFACE free_memory

         MODULE PROCEDURE free_wannier_u,free_wannier_P,free_v_pot, free_q_mat,  free_memory_ortho_polaw, &
              &free_memory_wp_psi, free_memory_wannier_u_prim, free_memory_v_pot_prim, free_memory_wp_psi_cutoff_index,&
              &free_memory_wp_psi_cutoff_data, free_memory_head_epsilon, free_cprim_prod, free_memory_upper_states

      END INTERFACE


      INTERFACE initialize_memory
         MODULE PROCEDURE initialize_memory_cprim_prod, initialize_memory_upper_states
      END INTERFACE



    CONTAINS

      subroutine free_wannier_u( r)
        implicit none
        type(wannier_u) :: r
        if(associated(r%ene))  deallocate(r%ene)
        nullify(r%ene)
        if(associated(r%umat)) deallocate(r%umat)
        nullify(r%umat)
        if(associated(r%ene_xc))  deallocate(r%ene_xc)
        nullify(r%ene_xc)
        if(associated(r%ene_u))  deallocate(r%ene_u)
        nullify(r%ene_u)
        if(associated(r%ene_lda_h))  deallocate(r%ene_lda_h)
        nullify(r%ene_lda_h)
        return
      end subroutine free_wannier_u


      subroutine free_wannier_P( w_P)
        implicit none
        type(wannier_P) :: w_P
        if(associated(w_P%ij)) deallocate(w_P%ij)
        nullify(w_P%ij)
        if(associated(w_P%o))  deallocate(w_P%o)
        nullify(w_P%o)
        return
      end subroutine free_wannier_P

      subroutine free_v_pot(vp)
        implicit none
        type(v_pot) :: vp
        if(associated(vp%vmat)) deallocate(vp%vmat)
        nullify(vp%vmat)
        return
      end subroutine free_v_pot


      subroutine free_q_mat( qm)
        implicit none
        type(q_mat) :: qm
        integer :: iw
        if(associated(qm%wp)) then
           do iw=1,qm%numpw_para
              call free_wannier_P(qm%wp(iw))
           enddo
           deallocate(qm%wp)
           nullify(qm%wp)
        endif
        return
      end subroutine free_q_mat

      SUBROUTINE free_memory_ortho_polaw(op)
        !this subroutine deallocates the green descriptor
        implicit none
        TYPE(ortho_polaw) op
        if(associated(op%on_mat)) deallocate(op%on_mat)
        nullify(op%on_mat)
        return
      END SUBROUTINE free_memory_ortho_polaw

      SUBROUTINE free_memory_wp_psi( wp)
        implicit none
        TYPE(wp_psi) :: wp
        if(associated(wp%wwp))  deallocate(wp%wwp)
        nullify(wp%wwp)
        return
      END SUBROUTINE free_memory_wp_psi

      SUBROUTINE  free_memory_wannier_u_prim( ww)
        implicit none
        TYPE(wannier_u_prim) :: ww
        if(associated(ww%umat)) deallocate(ww%umat)
        nullify(ww%umat)
        return
      END SUBROUTINE free_memory_wannier_u_prim

      SUBROUTINE free_memory_v_pot_prim(vp)
        implicit none
        TYPE(v_pot_prim) :: vp

        if(associated(vp%ij)) deallocate(vp%ij)
        nullify(vp%ij)
        if(associated(vp%vmat)) deallocate(vp%vmat)
        nullify(vp%vmat)
      END SUBROUTINE free_memory_v_pot_prim

      SUBROUTINE free_memory_wp_psi_cutoff_index(wpi)
        implicit none
        TYPE(wp_psi_cutoff_index) :: wpi
        if(associated(wpi%index)) deallocate(wpi%index)
        nullify(wpi%index)
        return
      END SUBROUTINE free_memory_wp_psi_cutoff_index

      SUBROUTINE free_memory_wp_psi_cutoff_data(wp)
        implicit none
        TYPE(wp_psi_cutoff_data) :: wp
        if(associated(wp%wwp)) deallocate(wp%wwp)
        nullify(wp%wwp)
        return
      END SUBROUTINE free_memory_wp_psi_cutoff_data

      SUBROUTINE free_memory_head_epsilon(he)
        implicit none
        TYPE(head_epsilon) :: he
        if(associated(he%freqs)) deallocate(he%freqs)
        nullify(he%freqs)
        if(associated(he%head)) deallocate(he%head)
        nullify(he%head)
        if(associated(he%gzero)) deallocate(he%gzero)
        nullify(he%gzero)
        if(associated(he%wing)) deallocate(he%wing)
        nullify(he%wing)
        if (associated(he%wing_c)) deallocate(he%wing_c)
        nullify(he%wing_c)


      END SUBROUTINE free_memory_head_epsilon


      SUBROUTINE free_cprim_prod(cpp)
        implicit none
        TYPE(cprim_prod) :: cpp

        if(associated(cpp%cpmat)) deallocate( cpp%cpmat)
        nullify(cpp%cpmat)

        return
      END SUBROUTINE free_cprim_prod

      SUBROUTINE initialize_memory_cprim_prod(cpp)
         implicit none
        TYPE(cprim_prod) :: cpp

        nullify(cpp%cpmat)

        return
      END SUBROUTINE initialize_memory_cprim_prod

      SUBROUTINE free_memory_upper_states(us)
        implicit none
        TYPE(upper_states) :: us
        if(associated(us%ene)) deallocate(us%ene)
        nullify(us%ene)
      END SUBROUTINE free_memory_upper_states

       SUBROUTINE initialize_memory_upper_states(us)
        implicit none
        TYPE(upper_states) :: us
        nullify(us%ene)
      END SUBROUTINE initialize_memory_upper_states

 END MODULE basic_structures



