###########################################################
# atom_info.x
###########################################################
set(sources SRC/atom_info.f90)
qe_add_executable(qe_qha_atominfo_exe ${sources})
set_target_properties(qe_qha_atominfo_exe PROPERTIES OUTPUT_NAME atom_info.x)

###########################################################
# Atom_projected_properties.x
###########################################################
set(sources SRC/Atom_projected_properties.f90)
qe_add_executable(qe_qha_atomprojectedproperties_exe ${sources})
set_target_properties(qe_qha_atomprojectedproperties_exe PROPERTIES OUTPUT_NAME Atom_projected_properties.x)

###########################################################
# F_QHA.x
###########################################################
set(sources SRC/F_QHA.f90)
qe_add_executable(qe_qha_fqha_exe ${sources})
set_target_properties(qe_qha_fqha_exe PROPERTIES OUTPUT_NAME F_QHA.x)

###########################################################
# Ghost_DOS.x
###########################################################
set(sources SRC/Ghost_DOS.f90)
qe_add_executable(qe_qha_ghostdos_exe ${sources})
set_target_properties(qe_qha_ghostdos_exe PROPERTIES OUTPUT_NAME Ghost_DOS.x)

###########################################################
# Mean_square_displacement.x
###########################################################
set(sources SRC/Mean_square_displacement.f90)
qe_add_executable(qe_qha_meansquaredisplacement_exe ${sources})
set_target_properties(qe_qha_meansquaredisplacement_exe PROPERTIES OUTPUT_NAME Mean_square_displacement.x)

###########################################################
# Partial_phonon_DOS.x
###########################################################
set(sources SRC/Partial_phonon_DOS.f90)
qe_add_executable(qe_qha_partialphonondos_exe ${sources})
set_target_properties(qe_qha_partialphonondos_exe PROPERTIES OUTPUT_NAME Partial_phonon_DOS.x)

###########################################################
# Debye.x
###########################################################
set(sources
    Debye/cheval.f
    Debye/d1mach.f
    Debye/debye3.f
    Debye/Debye.f90
    Debye/Debye_T.f)
qe_add_executable(qe_qha_debye_exe ${sources})
set_target_properties(qe_qha_debye_exe PROPERTIES OUTPUT_NAME Debye.x)

###########################################################
# phonon_dos.x
###########################################################
set(sources
    Phonon_DOS/phonon_dos.f
    Phonon_DOS/k_brillouin.f
    Phonon_DOS/generate_tetra.f 
    Phonon_DOS/det3.f
    Phonon_DOS/det4.f
    Phonon_DOS/Tetrahedra.f
    Phonon_DOS/Integration.f)
qe_add_executable(qe_qha_phonondos_exe ${sources})
set_target_properties(qe_qha_phonondos_exe PROPERTIES OUTPUT_NAME phonon_dos.x)

###########################################################
# tetra.x
###########################################################
set(sources
    Phonon_DOS/tetra.f
    Phonon_DOS/k_brillouin.f
    Phonon_DOS/generate_tetra.f 
    Phonon_DOS/det3.f
    Phonon_DOS/det4.f)
qe_add_executable(qe_qha_tetra_exe ${sources})
set_target_properties(qe_qha_tetra_exe PROPERTIES OUTPUT_NAME tetra.x)

###########################################################

qe_install_targets(
    qe_qha_atominfo_exe
    qe_qha_atomprojectedproperties_exe
    qe_qha_fqha_exe
    qe_qha_ghostdos_exe
    qe_qha_meansquaredisplacement_exe
    qe_qha_partialphonondos_exe
    qe_qha_debye_exe
    qe_qha_phonondos_exe
    qe_qha_tetra_exe)