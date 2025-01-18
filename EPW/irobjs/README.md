# How to Install IR Object Files

To use `sparse-ir` sampling for anisotropic Migdal-Eliashberg calculations, set the input value `gridsamp` to 2 and specify the file containing the IR functions in `filirobj`. 

By running the following command in this directory, you can download several files that contain precomputed IR functions:
```
make irobjs
```
The files are organized by the parameters Λ and ε, so make sure to set the corresponding file in `filirobj` according to the parameters you intend to use.

If you want to modify the parameters and calculate the IR functions yourself, you will need to install the `sparse-ir` Python package. After installation, you can run the following command:
```
python dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
```
For more details, please refer to the [sparse-ir-fortran GitHub repository](https://github.com/SpM-lab/sparse-ir-fortran).
