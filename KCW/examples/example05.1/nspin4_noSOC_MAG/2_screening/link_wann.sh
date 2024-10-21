#!/bin/bash

rm -fr wann_*
ln -s ../1_wannier/occ/wann_u.mat wann_u.mat
ln -s ../1_wannier/occ/wann_centres.xyz wann_centres.xyz
ln -s ../1_wannier/emp/wann_u.mat wann_emp_u.mat
ln -s ../1_wannier/emp/wann_u_dis.mat wann_emp_u_dis.mat
ln -s ../1_wannier/emp/wann_centres.xyz wann_emp_centres.xyz

