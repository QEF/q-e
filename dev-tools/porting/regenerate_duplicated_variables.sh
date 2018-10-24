#!/bin/bash

python gen_intrinsic.py gvect "gg:REAL(DP):1" "g:REAL(DP):2" "mill:INTEGER:2" "eigts1:COMPLEX(DP):2" "eigts2:COMPLEX(DP):2" "eigts3:COMPLEX(DP):2" 
mv gvect_gpu.f90 ../../Modules/recvec_gpu.f90

python gen_intrinsic.py wavefunctions "evc:COMPLEX(DP):2" "psic:COMPLEX(DP):1" "psic_nc:COMPLEX(DP):2"
mv wavefunctions_gpu.f90 ../../Modules/wavefunctions_gpu.f90

python gen_intrinsic.py wvfct "g2kin:REAL(DP):1" "et:REAL(DP):2" "wg:REAL(DP):2"
mv wvfct_gpu.f90 ../../PW/src/pwcom_gpu.f90

python gen_intrinsic.py us "qrad:REAL(DP):4" "tab:REAL(DP):3" "tab_at:REAL(DP):3" "tab_d2y:REAL(DP):3"
cat us_gpu.f90 >> ../../PW/src/pwcom_gpu.f90
rm us_gpu.f90

python gen_intrinsic.py spin_orb "fcoef:COMPLEX(DP):5"
cat spin_orb_gpu.f90 >> ../../PW/src/pwcom_gpu.f90
rm spin_orb_gpu.f90

python gen_intrinsic.py g_psi_mod "h_diag:REAL(DP):2" "s_diag:REAL(DP):2"
mv g_psi_mod_gpu.f90 ../../PW/src/g_psi_mod_gpu.f90

python gen_intrinsic.py scf "vrs:REAL(DP):2"
mv scf_gpu.f90 ../../PW/src/scf_mod_gpu.f90

python gen_intrinsic.py uspp "indv:INTEGER:2" "nhtol:INTEGER:2" "nhtolm:INTEGER:2" \
                             "ijtoh:INTEGER:3" "indv_ijkb0:INTEGER:1" \
                             "vkb:COMPLEX(DP):2" "becsum:REAL(DP):3" "ebecsum:REAL(DP):3" \
                             "dvan:REAL(DP):3" "deeq:REAL(DP):4" \
                             "qq_nt:REAL(DP):3" "qq_at:REAL(DP):3" \
                             "nhtoj:REAL(DP):2" \
                             "qq_so:COMPLEX(DP):4" "dvan_so:COMPLEX(DP):4" \
                             "deeq_nc:COMPLEX(DP):4"
mv uspp_gpu.f90 ../../Modules/uspp_gpu.f90

python gen_derived.py becmod
mv becmod_gpu.f90 ../../Modules/becmod_gpu.f90


