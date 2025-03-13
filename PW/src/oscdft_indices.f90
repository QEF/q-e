MODULE oscdft_indices
#if defined (__OSCDFT)
   USE kinds, ONLY : DP

   PRIVATE
   PUBLIC oscdft_constr_indices_type,&
          oscdft_orbital_indices_type,&
          oscdft_indices_type

   TYPE oscdft_constr_indices_type
      INTEGER              :: norbs
      INTEGER, ALLOCATABLE :: icorb2iorb(:),& ! constr%norbitals
                              iorb2icorb(:),& ! inp%norbs
                              icorb_start(:),& ! nconstr
                              icorb_end(:),& ! nconstr
                              ins2iorb(:,:),& ! (max_ns_dim,nconstr) rows/cols of occup matrix to constr orbital index (icorb)
                              ins2ioff(:,:) ! (max_ns_dim,nconstr) rows/cols of occup matrix to quantum number m
   END TYPE oscdft_constr_indices_type

   TYPE oscdft_orbital_indices_type
      INTEGER              :: norbs
      INTEGER, ALLOCATABLE :: iat(:),& ! norbitals atom index
                              n(:),&! norbitals
                              l(:),&! norbitals
                              iat_sym(:,:),& ! (nsym,norbs)
                              iorb_start(:),& ! noscdft
                              iorb_end(:),& !noscdft
                              iorb2ioscdft(:),& ! norbs
                              ins2iorb(:,:),& ! (max_ns_dim,noscdft) rows/cols of occup matrix to orbital index (iorb)
                              ins2ioff(:,:) ! (max_ns_dim,noscdft) rows/cols of occup matrix to quantum number m
   END TYPE oscdft_orbital_indices_type

   TYPE oscdft_indices_type
      INTEGER              :: max_ns_dim, nconstr
      INTEGER, ALLOCATABLE :: iconstr2ioscdft(:),&
                              ioscdft2iconstr(:),&
                              ns_dim(:),&
                              nchi(:,:)
      TYPE(oscdft_constr_indices_type)  :: constr
      TYPE(oscdft_orbital_indices_type) :: orbs

      REAL(DP),    ALLOCATABLE :: overlap_gam(:,:,:,:,:) ! max_ns_dim,max_ns_dim,nsym,noscdft,nkstot
      COMPLEX(DP), ALLOCATABLE :: overlap_k  (:,:,:,:,:) ! max_ns_dim,max_ns_dim,nsym,noscdft,nkstot
      COMPLEX(DP), ALLOCATABLE :: coeffs(:,:,:,:,:) ! max_ns_dim,max_ns_dim,nsym,noscdft,nks
   END TYPE oscdft_indices_type
#endif
END MODULE oscdft_indices
