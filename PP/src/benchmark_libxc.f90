!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
PROGRAM benchmark_libxc
  !
  !**--------------------------------------------------------------------------------**
  !** REMEMBER to comment the libxc blocks in 'Modules/correlation_lda_lsda.f90' and **
  !** 'Modules/exchange_lda_lsda.f90' in order to run consistent tests.              ** 
  !**--------------------------------------------------------------------------------**
  !
#if defined(__LIBXC)
  !
  USE xc_f90_types_m
  USE xc_f90_lib_m
  !
  USE xc_lda_lsda
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER     :: DP = selected_real_kind(14,200)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info1, xc_info2
  REAL(DP), ALLOCATABLE  :: ex_lda(:) , vx_lda(:), ec_lda(:) , vc_lda(:)
  !-----------------------------
  REAL(DP), ALLOCATABLE  :: rho(:), ex(:), vx(:), ec(:), vc(:)
  CHARACTER(LEN=120) :: name1, name2
  INTEGER :: ii,nnr
  INTEGER :: iexch_qe, icorr_qe, iexch_xc, icorr_xc
  !
  !
  nnr=10
  ALLOCATE( rho(1:nnr) )
  ALLOCATE( vx(1:nnr), ex(1:nnr) )
  ALLOCATE( vc(1:nnr), ec(1:nnr) )
  ALLOCATE( vx_lda(1:nnr), ex_lda(1:nnr) )
  ALLOCATE( vc_lda(1:nnr), ec_lda(1:nnr) )
  !
  DO ii = 1, nnr
     rho(ii)=DBLE(ii)/DBLE(nnr)
  ENDDO
  !
  ex_lda = 0._dp  ;  ec_lda = 0._dp
  vx_lda = 0._dp  ;  vc_lda = 0._dp
  ex = 0._dp      ;  ec = 0._dp
  vx = 0._dp      ;  vc = 0._dp
  !
  !**--------------------------------------------------------------------------------**
  !** libxc funct. indexes: http://bigdft.org/Wiki/index.php?title=XC_codes          **
  !** qe      "       "   : see comments in Modules/funct.f90                        **
  !**--------------------------------------------------------------------------------**
  !**                                                                                **
  !**  ... some examples:                                                            **
  !**                                                                                **
  !**                 |   q-e     |     libxc    |                                   **
  !**                 |___________|______________|                                   **
  !**     exchange:   |    1      |       1      |                                   **
  !**     pz          |    1      |       9      |                                   **
  !**     wigner      |    5      |       2      |                                   **
  !**     vwn         |    2      |       7      |                                   **
  !**     pw          |    4      |      12      |                                   **
  !**                                                                                **
  !**--------------------------------------------------------------------------------**
  !
  !
  PRINT *, "Insert functional indexes"
  PRINT *, " "
  PRINT *, "iexch_xc, icorr_xc: "
  READ(*,*) iexch_xc, icorr_xc
  PRINT *, "iexch_qe, icorr_qe: "
  READ(*,*) iexch_qe, icorr_qe
  !
  !------ LIBXC ------
  !
  CALL xc_f90_func_init( xc_func, xc_info1, iexch_xc, XC_UNPOLARIZED )
   CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho(1), ex_lda(1), vx_lda(1) )
  CALL xc_f90_func_end( xc_func )
  !
  CALL xc_f90_func_init( xc_func, xc_info2, icorr_xc, XC_UNPOLARIZED )
   CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho(1), ec_lda(1), vc_lda(1) )
  CALL xc_f90_func_end( xc_func )
  !
  !----- QE ----------
  !
  CALL select_functionals( iexch_qe, icorr_qe )
  CALL xc( nnr, rho, ex, ec, vx, vc )
  !CALL xc_spin( nnr, rho, zeta, ex, ec, vx, vc )
  !
  !------------------
  !
  CALL xc_f90_info_name(xc_info1,name1)
  CALL xc_f90_info_name(xc_info2,name2)
  PRINT *, " "
  PRINT *, "libxc:  ", trim(name1), " + ", trim(name2)
  PRINT *, " "
  DO ii = 1, nnr, nnr
     PRINT *, "rho:",rho(ii)
     PRINT *, "--- ex, vx ---"
     PRINT *, "libxc:",ex_lda(ii),vx_lda(ii)
     PRINT *, "libqe:",ex(ii),vx(ii)
     PRINT *, "diffs:",ex_lda(ii)-ex(ii),vx_lda(ii)-vx(ii)
     PRINT *, "--- ec, vc ---"
     PRINT *, "libxc:",ec_lda(ii),vc_lda(ii)
     PRINT *, "libqe:",ec(ii),vc(ii)
     PRINT *, "diffs:",ec_lda(ii)-ec(ii),vc_lda(ii)-vc(ii)
     PRINT *, "--- ---"
  ENDDO
  !
  !
  DEALLOCATE(rho)
  DEALLOCATE(ex)
  DEALLOCATE(ec)
  DEALLOCATE(vx)
  DEALLOCATE(vc)
  DEALLOCATE(ex_lda)
  DEALLOCATE(vx_lda)
  DEALLOCATE(ec_lda)
  DEALLOCATE(vc_lda)
  !
#else
  !
  PRINT *, "ERROR: library libxc not found. Check Makefile."
  !
#endif
  !
END PROGRAM benchmark_libxc
