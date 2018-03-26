# FFTXlib

Implements real space grid parallelization of FFT and task groups. 

## Testing and Benchmarking

This library also provides a testing and timing code to asses the performance of your FFT, estimate the
scalability and the optimal parameters for your simulation.

To compile the test program, once you have properly configure QE within a parallel environment,
go inside the directory FFTXlib and type:

    make TEST

Then you can run your FFT tests using command like:

    mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4

Command line arguments:

    -ecutwfc  Plane wave energy cut off
    -alat     Lattice parameter (for hard coded lattice structure)
    -nbnd     Number of bands (fft cycles)
    -ntg      Number of task groups
    -av1  x y z    First lattice vector, in atomic units. N.B.: when using -av1, -alat is ignored!
    -av2  x y z    Second lattice vector, in atomic units. N.B.: when using -av2, -alat is ignored!
    -av3  x y z    Third lattice vector, in atomic units. N.B.: when using -av3, -alat is ignored!
    -kmax kx ky kz    Reciprocal lattice vector inside the BZ with maximum norm. Used to calculate max(|G+K|). (2pi/a)^2 units.

A python script to extract the parameters from an output file of pw.x is also available. Example usage:

    $ python gen_test_params.py a_pw_output
    To analize performances run with:
    mpirun -np X ./fft_test.x -ntg Y -ecutwfc 36.7500 -ecutrho 147.0000 -av1 36.6048 0.0 0.0 -av2 -18.3024 31.70067192 0.0 -av3 0.0 0.0 18.3024 -nbnd 400 -gamma .true.

Replace `X` and `Y` with appropriate values for your simualtion.
    
