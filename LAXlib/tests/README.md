# LAXlib testing suite

In order to run the tests first run `./configure` in QE topdir and generate a
valid `make.inc`.
You may also download large eigenvalue problems (in this directory) with:

    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1EAB3xkoD-i9p4nW6NJDED3WaEK8ZCcf4' -O SiGeK1.bin
    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=13lFkDbv99V8fqiXER1N2IzoJ_EhuGtt9' -O SiGeK2.bin

Finally `make` and run all `.x` executable files.
