#
# SP: This test is a bit specific and involves a Python script with spglib dependence
#

The dependendencies have been checked with a clean conda environment
conda create --name test python=3.9
conda activate test

pip install spglib

multipole.py can be run with (suggsted) or without (default) ase

In case you are on a cluster without remote access, you may have to install it differently.
In which case, you also need to comment out the pip install line in ../run-ph.sh and look for case number 12.

