#!/bin/bash

cd espresso/02-Wfn
mkdir -p lif.save
cd ../../

cd espresso/03-Wfnq
mkdir -p lif.save
cd ../../

cd espresso/04-Wfn_co
mkdir -p lif.save
cd ../../

cd espresso/05-Wfn_fi
mkdir -p lif.save
cd ../../

cd espresso/06-Wfnq_fi
mkdir -p lif.save
cd ../../

ln -nfs ../../01-Density/lif.save/charge-density.hdf5 espresso/02-Wfn/lif.save/charge-density.hdf5
ln -nfs ../../01-Density/lif.save/spin-polarization.hdf5 espresso/02-Wfn/lif.save/spin-polarization.hdf5

ln -nfs ../../01-Density/lif.save/charge-density.hdf5 espresso/03-Wfnq/lif.save/charge-density.hdf5
ln -nfs ../../01-Density/lif.save/spin-polarization.hdf5 espresso/03-Wfnq/lif.save/spin-polarization.hdf5

ln -nfs ../../01-Density/lif.save/charge-density.hdf5 espresso/04-Wfn_co/lif.save/charge-density.hdf5
ln -nfs ../../01-Density/lif.save/spin-polarization.hdf5 espresso/04-Wfn_co/lif.save/spin-polarization.hdf5

ln -nfs ../../01-Density/lif.save/charge-density.hdf5 espresso/05-Wfn_fi/lif.save/charge-density.hdf5
ln -nfs ../../01-Density/lif.save/spin-polarization.hdf5 espresso/05-Wfn_fi/lif.save/spin-polarization.hdf5

ln -nfs ../../01-Density/lif.save/charge-density.hdf5 espresso/06-Wfnq_fi/lif.save/charge-density.hdf5
ln -nfs ../../01-Density/lif.save/spin-polarization.hdf5 espresso/06-Wfnq_fi/lif.save/spin-polarization.hdf5


ln -nfs ../espresso/02-Wfn/wfn.cplx 11-epsilon/WFN
ln -nfs ../espresso/03-Wfnq/wfn.cplx 11-epsilon/WFNq
ln -nfs ../espresso/04-Wfn_co/wfn.cplx 12-sigma/WFN_inner
ln -nfs ../espresso/04-Wfn_co/rho.real 12-sigma/RHO
ln -nfs ../espresso/04-Wfn_co/vxc.dat 12-sigma/vxc.dat


ln -nfs ../11-epsilon/eps0mat.h5 12-sigma/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 12-sigma/epsmat.h5
ln -nfs ../espresso/04-Wfn_co/wfn.cplx 13-kernel/WFN_co
ln -nfs ../11-epsilon/eps0mat.h5 13-kernel/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 13-kernel/epsmat.h5
ln -nfs ../espresso/04-Wfn_co/wfn.cplx 14-absorption/WFN_co
ln -nfs ../espresso/05-Wfn_fi/wfn.cplx 14-absorption/WFN_fi
ln -nfs ../espresso/06-Wfnq_fi/wfn.cplx 14-absorption/WFNq_fi
ln -nfs ../12-sigma/sigma_hp.log 14-absorption/sigma_hp.log
ln -nfs ../12-sigma/eqp1.dat 14-absorption/eqp_co.dat
ln -nfs ../11-epsilon/eps0mat.h5 14-absorption/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 14-absorption/epsmat.h5
ln -nfs ../13-kernel/bsemat.h5 14-absorption/bsemat.h5
