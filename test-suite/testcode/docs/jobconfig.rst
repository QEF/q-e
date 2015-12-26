.. _jobconfig:

jobconfig
=========

The jobconfig file defines the tests to run.  If a section named 'categories'
exists, then it gives labels to sets of tests.  All other sections are assumed
to individually define a test.

Tests
-----

A test is assumed to reside in the directory given by the name of the test
section.  For example::

    [carbon_dioxide_ccsd]
    inputs_args = ('co2.inp','')

would define a test in the ``carbon_dioxide_ccsd`` subdirectory relative to the
``jobconfig`` configuration file, with the input file as ``co2.inp`` (in the
``carbon_dioxide_ccsd`` subdirectory) with no additional arguments to be passed
to the test program.  All input and output files related to the test are
assumed to be contained within the test subdirectory.

The following options are permitted:

inputs_args [inputs and arguments format (see :ref:`below <inputs>`)]
    Input filename and associated arguments to be passed to the test program.
    No default.
min_nprocs [integer]
    Minimum number of processors to run test on.  Cannot be overridden by the
    '--processors' command-line option.  Default: 0.
max_nprocs [integer]
    Maximum number of processors to run test on.  Cannot be overridden by the
    '--processors' command-line option.  Default: 2^31-1 or 2^63-1.
nprocs [integer]
    Number of processors to run the test on.  Zero indicates to run the test
    purely in serial, without using an external program such as mpirun to
    launch the test program.  Default: 0.
output [string]
    Filename to which the output is written if the output is not written to
    standard output.  The output file is moved to the specific testcode test
    filename at the end of the calculation before the test output is validated
    against the benchmark output.  Wildcards are allowed so long as the pattern
    only matches a single file at the end of the calculation.  Default:
    inherits from setting in :ref:`userconfig`.
path [string]
    Set path (relative to the directory containing the ``jobconfig``
    configuration file) of the test.  The test is run in this directory and so
    input filenames need to be relative to it.  If the given path contains
    wildcards, then this is expanded and an individual test is created for each
    path that maches the pattern.  Note that Python's configparser restricts
    the use of special characters in section names and hence some patterns can
    only be accomplished by explicitly using the path option.  Default: test
    name (i.e.  the name of the section defining the test).
run_concurrent [boolean]
    If true then subtests defined by the inputs_args option are allowed to run
    concurrently rather than consecutively, assuming enough processors are
    available.  Default: false.
submit_template [string]
    Path to a template of a submit script used to submit jobs to a queueing
    system.  testcode will replace the string given in submit_pattern with the
    command(s) to run the test.  The submit script must do all other actions (e.g.
    setting environment variables, loading modules, copying files from the test
    directory to a local disk and copying files back afterwards).  No default.
program [string]
    Program name (appropriate section heading in :ref:`userconfig`) to use to
    run the test.  Default: specified in the [user] section of
    :ref:`userconfig`.
tolerance [tolerance format (see :ref:`tolerance`)]
    Tolerances for comparing test output to the benchmark output.  Default:
    inherits from the settings in :ref:`userconfig`.

If a test is defined via a category/path containing wildcards and explicitly,
then the explicit category will inherit any settings from the wildcard
definition.  For example, given the subdirectories ``t1`` and ``t2``, each
containing tests, the definition::

    [t*]
    inputs_args = ('test.in', '')
    [t1]
    nprocs = 2

is entirely equivalent to::

    [t1]
    nprocs = 2
    inputs_args = ('test.in', '')
    [t2]
    inputs_args = ('test.in', '')

.. note::

    Explicitly defining a test multiple times, e.g.::

        [t1]
        inputs_args = ('inp1', '')
        [t1]
        inputs_args = ('inp2', '')

    is not permitted and the resultant settings are not uniquely defined.

Test categories
---------------

For the purposes of selecting a subset of the tests in :ref:`testcode.py`, each
test is automatically placed in two separate categories, one labelled by the
test's name and the other by the test's path.  A test can hence be referred to
by either its path or by its name (which are identical by default).  

Additional categories can be specified in the [categories] section.  This makes
it very easy to select subsets of the tests to run.  For example::

    [categories]
    cat1 = t1 t2
    cat2 = t3 t4
    cat3 = cat1 t3

defines three categories (`cat`, `cat2` and `cat3`), each containing a subset
of the overall tests.  A category may contain another category so long as
circular dependencies are avoided.  There are two special categories, `_all_`
and `_default_`.  The `_all_` category contains, by default, all tests and
should not be changed under any circumstances.  The `_default_` category can
be set; if it is not specified then it is set to be the `_all_` category.

.. _inputs:

Program inputs and arguments
----------------------------

The inputs and arguments must be given in a specific format.  As with the
:ref:`tolerance format <tolerance>`,  the inputs and arguments are specified
using a comma-separated list of python tuples.  Each tuple (basically
a comma-separated list enclosed in parantheses) contains two elements: the name
of an input file and the associated arguments, in that order, represents
a subtest belonging to the given test.  Both elements must be quoted.  If the
input filename contains wildcard, then those wildcards are expanded to find all
files in the test subdirectory which match that pattern; the expanded list is
sorted in alphanumerical order.  A separate subtest (with the same arguments
string) is then created for each file matching the pattern.  used to construct
the command to run.  A null string (``''``) should be used to represent the
absence of an input file or arguments.  By default subtests run in the order
they are specified.  For example::

    inputs_args = ('test.inp', '')

defines a single subtest, with input filename ``test.inp`` and no arguments,

::

    inputs_args = ('test.inp', ''), ('test2.inp', '--verbose')

defines two subtests, with an additional argument for the second subtest, and

::

    inputs_args = ('test*.inp', '')

defines a subtest for each file matching the pattern ``test*inp`` in the
subdirectory of the test.
