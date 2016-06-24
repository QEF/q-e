#############################
# Test 01 : B-doped diamond #
#############################
The following feature of the code are tested
- Correct unfolding from IBZ to full BZ
- Correct Wannier interpolation
- Phonon & electron self-energy
- Eliashberg a2F
- Homogeneous fine k and q-grid integration
- Test nesting function
- Test spectral function
- Test parallel_k (epw1.in)
- Test parallel_q (epw2.in)
- Test restart feature epwread = .true. (epw2.in)
- Test band_plot (epw3.in)
- Test iverbosity = 1 (epw4.in)

#################
# Test 02 : SiC #
#################
The following feature of the code is tested
- Time-reversal symmetry when inversion sym. is not part of the small group of q.
- Test memory optimization etf_mem = .false. 

################
# Test 03 : Pb #
################
The following feature of the code is tested
- Test metals
- Wannier BS plotting to check Wannier interpolation

########################
# Test 04 : Pb with SO #
########################
The following feature of the code is tested
- Spin-orbit coupling
- Wannier BS plotting to check Wannier interpolation

##################
# Test 05 : MgB2 #
##################
The following feature of the code is tested
- Superconducting properties within the Migdal-Eliashberg framework
- Isotropic superconducting properties with contour integration
- Isotropic superconducting properties along the real axis
- Anisotropic superconducting properties with Pade continuation


#################
# Test 06 : GaN #
#################
The following feature of the code is tested
- Polar correction
- Test reading fine q-grid from file





