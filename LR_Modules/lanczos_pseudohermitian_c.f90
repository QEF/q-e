subroutine lanczos_pseudohermitian_c(j, npwx_npol, nbnd_occ, nksq, qj_r, Aqj_r, qjold_r, n_ipol, u, alpha, beta, gamma, zeta)
   !
   !! Bi-Orthogonal Lanczos algorithm for magnons 
   !! 
   !! $$ g(w) =\sum_j (u,q_j){q_j,(w-A)^(-1)q}
   !! Algorithm 1 in "Computer Physics Communications 185 (2014) 2080-2089"
   !!
   !! this subroutine generates alpha, beta and gamma coefficients (the
   !! tridiagonal matrix elements), z = (u,q_j) elements. And update 
   !! Lanczos vectors
   !
   USE kinds,                    ONLY : dp
   !
   INTEGER, INTENT(IN)               :: j
   !! iteration index
   INTEGER, INTENT(IN)               :: npwx_npol
   !! firts dimension of qj, Aqj, SAqj, qjold, u in qe npwx*npol
   INTEGER, INTENT(IN)               :: nbnd_occ
   !! second dimension of qj, Aqj, SAqj, qjold, n_ipol, u in qe nbnd
   INTEGER, INTENT(IN)               :: nksq
   !! third dimension of qj, Aqj, SAqj, qjold, n_ipol, u in qe nksq
   INTEGER, INTENT(IN)               :: n_ipol
   !! polarization, forth dimension of u and dimension of zeta
   COMPLEX(kind=dp), INTENT(INOUT)   :: Aqj_r(npwx_npol,nbnd_occ,nksq, 2)
   !! operator applied to right qj vector 
   COMPLEX(kind=dp), INTENT(IN)      :: u(npwx_npol,nbnd_occ,nksq, 2, n_ipol)
   !! second lanczos vector, in qe O_psi for magnons
   COMPLEX(kind=dp), INTENT(INOUT)   :: qj_r(npwx_npol,nbnd_occ,nksq, 2)
   !! right qj vector, become right qj+1 vector
   COMPLEX(kind=dp), INTENT(INOUT)   :: qjold_r(npwx_npol,nbnd_occ,nksq, 2)
   !! right qj-1 vector, become right qj vector
   COMPLEX(kind=dp), INTENT(OUT)     :: alpha
   !! diagonal cofficient of the tridiagonal matrix
   REAL(kind=dp),    INTENT(OUT)     :: beta
   !! lower coefficient of the tridiagonal matrix
   COMPLEX(kind=dp), INTENT(OUT)     :: gamma
   !! upper coefficient of the tridiagonal matrix
   COMPLEX(kind=dp), INTENT(OUT)     :: zeta(n_ipol)
   !! (u,q_j) products
   ! dummy variable
   COMPLEX(kind=dp)   :: Aqj_l(npwx_npol,nbnd_occ,nksq, 2)
   !   
   COMPLEX(kind=dp),EXTERNAL :: lr_dot_magnons
   !
   INTEGER                :: size_evc, ip
   !
   size_evc = npwx_npol*nbnd_occ*nksq*2
   !
   ! Orthogonality requirement: <v|\bar{L}|v> = 1
   ! computes gamma and beta
   !
   qj_r(:,:,:,2) = - qj_r(:,:,:,2)
   gamma = lr_dot_magnons(qj_r, Aqj_r)
   qj_r(:,:,:,2) = - qj_r(:,:,:,2)
   !
   beta = sqrt(abs(gamma))
   !
   gamma = gamma / beta
   !
   ! Renormalize q(i) and Lq(i), also p(i) and Lp(i) in the non-Hermitian case
   !
   CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),qj_r(1,1,1,1),1)
   CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),Aqj_r(1,1,1,1),1)
   !
   !
   ! computes alpha 
   !
   alpha = (0.0d0, 0.0d0)
   Aqj_l = Aqj_r
   Aqj_l(:,:,:,2) = -Aqj_l(:,:,:,2)
   !
   alpha = lr_dot_magnons(Aqj_l, Aqj_r)
   !
   ! Calculation of zeta coefficients.
   ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
   !
   DO ip = 1, n_ipol
      !
      ! Optics: In the ultrasoft case, the S operator was already
      ! applied to d0psi, so we have <S*d0psi|evc1>.
      !
      zeta(ip) = (0.0d0,0.0d0)
      zeta(ip) = lr_dot_magnons(qj_r, u(:,:,:,:,ip))
      !
   ENDDO
   !
   ! Aqj_r = Aqj_r -alpha*qj_r -gamma*qjold_r
   ! Renormalization will be done in the begining of the next iteration.
   ! In the non-Hermitian case, similar operation needs to be done also for p(i).
   !
   CALL zaxpy(size_evc,-alpha,qj_r(1,1,1,1),1,Aqj_r(1,1,1,1),1)
   CALL zaxpy(size_evc,-gamma,qjold_r(1,1,1,1),1,Aqj_r(1,1,1,1),1)
   !
   ! X. Ge: Throw away q(i-1), and make q(i+1) to be the current vector,
   ! be ready for the next iteration. Aqj will be free again after this
   ! step
   !
   CALL zcopy(size_evc,qj_r(1,1,1,1),1,qjold_r(1,1,1,1),1) ! qjold_r = qj_r
   CALL zcopy(size_evc,Aqj_r(1,1,1,1),1,qj_r(1,1,1,1),1)    ! qj_r = Aqj_r
   !
end subroutine lanczos_pseudohermitian_c
