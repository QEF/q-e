###########################################################
# atom_info.x
###########################################################
set(src_atom_info_x SRC/atom_info.f90)
qe_add_executable(qe_qha_atominfo_exe ${src_atom_info_x})
set_target_properties(qe_qha_atominfo_exe 
    PROPERTIES 
        OUTPUT_NAME atom_info.x)

###########################################################
# Atom_projected_properties.x
###########################################################
set(src_Atom_projected_properties_X SRC/Atom_projected_properties.f90)
qe_add_executable(qe_qha_atomprojectedproperties_exe ${src_Atom_projected_properties_X})
set_target_properties(qe_qha_atomprojectedproperties_exe 
    PROPERTIES 
        OUTPUT_NAME Atom_projected_properties.x)

###########################################################
# F_QHA.x
###########################################################
set(src_F_QHA_x SRC/F_QHA.f90)
qe_add_executable(qe_qha_fqha_exe ${src_F_QHA_x})
set_target_properties(qe_qha_fqha_exe 
    PROPERTIES 
        OUTPUT_NAME F_QHA.x)

###########################################################
# Ghost_DOS.x
###########################################################
set(src_Ghost_DOS_x SRC/Ghost_DOS.f90)
qe_add_executable(qe_qha_ghostdos_exe ${src_Ghost_DOS_x})
set_target_properties(qe_qha_ghostdos_exe 
    PROPERTIES 
        OUTPUT_NAME Ghost_DOS.x)

###########################################################
# Mean_square_displacement.x
###########################################################
set(src_Mean_square_displacement_x SRC/Mean_square_displacement.f90)
qe_add_executable(qe_qha_meansquaredisplacement_exe ${src_Mean_square_displacement_x})
set_target_properties(qe_qha_meansquaredisplacement_exe 
    PROPERTIES 
        OUTPUT_NAME Mean_square_displacement.x)

###########################################################
# Partial_phonon_DOS.x
###########################################################
set(src_Partial_phonon_DOS_x SRC/Partial_phonon_DOS.f90)
qe_add_executable(qe_qha_partialphonondos_exe ${src_Partial_phonon_DOS_x})
set_target_properties(qe_qha_partialphonondos_exe 
    PROPERTIES 
        OUTPUT_NAME Partial_phonon_DOS.x)

###########################################################
# Debye.x
###########################################################
set(src_Debye_x
    Debye/cheval.f
    Debye/d1mach.f
    Debye/debye3.f
    Debye/Debye.f90
    Debye/Debye_T.f)
qe_add_executable(qe_qha_debye_exe ${src_Debye_x})
set_target_properties(qe_qha_debye_exe 
    PROPERTIES 
        OUTPUT_NAME Debye.x)

###########################################################
# phonon_dos.x
###########################################################
set(src_phonon_dos_x
    Phonon_DOS/phonon_dos.f
    Phonon_DOS/k_brillouin.f
    Phonon_DOS/generate_tetra.f 
    Phonon_DOS/det3.f
    Phonon_DOS/det4.f
    Phonon_DOS/Tetrahedra.f
    Phonon_DOS/Integration.f)
qe_add_executable(qe_qha_phonondos_exe ${src_phonon_dos_x})
set_target_properties(qe_qha_phonondos_exe 
    PROPERTIES 
        OUTPUT_NAME phonon_dos.x)

###########################################################
# tetra.x
###########################################################
set(src_tetra_x
    Phonon_DOS/tetra.f
    Phonon_DOS/k_brillouin.f
    Phonon_DOS/generate_tetra.f 
    Phonon_DOS/det3.f
    Phonon_DOS/det4.f)
qe_add_executable(qe_qha_tetra_exe ${src_tetra_x})
set_target_properties(qe_qha_tetra_exe 
    PROPERTIES 
        OUTPUT_NAME tetra.x)

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
