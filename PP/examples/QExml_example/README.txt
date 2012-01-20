
===========================================
QEXML example (courtesy of Andrea Ferretti)
===========================================

list of actions to perform (~espresso = root directory of QE):

- compile pw.x
- go to ~espresso/PP/, compile qexml.x  (type "make qexml.x")
- come back to this directory
- run scf.in using pw.x to produce a .save dir (verify that pseudo_dir 
  points to the directory containing pseudopotentials)

- run qexml.x: : ~espresso/PP/qexml.x < example.in

- convert to textual some of the *.dat files to double check
  that everything is consistent; in particular
    silicon.save/gvectors.dat
    silicon.save/K00002/gkvectors.dat
    silicon.save/K00002/evc.dat

  to do this run the iotk script
    ~espresso/bin/iotk convert file.dat  file.xml

- have a look at ~espresso/PP/qexml_example.f90 to see the use of
  qexml.f90 routines
  
