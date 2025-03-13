# Common setup
* scf.in, ph.in, q2r.in, matdyn.in, scf.in

# Testset 1: 4 bands (valence only)
* nscf1.in, ahc1.in, epw11.in, epw12.in
* lifc = .false.
* vme = 'dipole'
* epw11.in: Use `lopt_w2b = .true.`
* epw12.in: use elecselfen_type = 'adiabatic'

# Testset 2: 8 bands (valence + conduction)
* nscf2.in, ahc2.in, epw21.in, epw22.in, epw23.in
* Set ahc_win_max to a finite value, use bands_skipped, use disentanglement
* lifc = .true.
* vme = 'wannier'
* epw23.in: Compute self-energy and spectral function. Use `lopt_w2b = .true.`