.. _userconfig:

userconfig
==========

The userconfig file must contain at least two sections.  One section must be
entitled 'user' and contains various user settings.  Any other section is
assumed to define a program to be tested, where the program is referred to
internally by its section name.  This makes it possible for a set of tests to
cover multiple, heavily intertwined, programs.  It is, however, far better to
have a distinct set of tests for each program where possible.

[user] section
--------------

The following options are allowed in the [user] section:

benchmark [string]
    Specify the ID of the benchmark to compare to.  This should be set running

    .. code-block bash

        $ testcode.py make-benchmarks

    The format of the benchmark files is'benchmark.out.ID.inp=INPUT_FILE.arg=ARGS'.  
    The 'inp' and/or 'arg' section is not included if it is empty.

    Multiple benchmarks can be used by providing a space-separated list of IDs.  The first
    ID in the list which corresponds to an existing benchmark filename is used to
    validate the test.
date_fmt [string]
    Format of the date string used to uniquely label test outputs.  This must
    be a valid date format string (see `Python documenation
    <http://docs.python.org/library/time.html>`_).  Default: %d%m%Y.
default_program [string]
    Default program used to run each test.  Only needs to be set if
    multiple program sections are specified.  No default.
diff [string]
    Program used to diff test and benchmark outputs.  Default: diff.
tolerance [tolerance format (see :ref:`below <tolerance>`.)]
    Default tolerance(s) used to compare all tests to their respective
    benchmarks.  Default: absolute tolerance 10^-10; no relative tolerance set.

[program_name] section(s)
-------------------------

The following options are allowed to specify a program (called 'program_name')
to be tested:

data_tag [string]
    Data tag to be used to extract data from test and benchmark output.  See
    :ref:`verification` for more details.  No default.
ignore_fields [space-separated list of strings]
    Specify the fields (e.g. column headings in the output from the extraction
    program) to ignore.  This can be used to include, say, timing information
    in the test output for performance comparison without causing failure of
    tests.  No default.
exe [string]
    Path to the program executable.  No default.
extract_args [string]
    Arguments to supply to the extraction program.  Default: null string. 
extract_cmd_template [string]
    Template of command used to extract data from output(s) with the following
    substitutions made:

        tc.extract
            replaced with the extraction program.
        tc.args
            replaced with extract_args.
        tc.file
            replaced with (as required) the filename of the test output or the
            filename of the benchmark output.
        tc.bench
            replaced with the filename of the benchmark output.
        tc.test
            replaced with the filename of the test output.

    Default: tc.extract tc.args tc.file if verify is False and
    tc.extract tc.args tc.test tc.bench if verify is True.
extract_program [string]
    Path to program to use to extract data from test and benchmark output.
    See :ref:`verification` for more details.  No default.
extract_fmt [string]
    Format of the data returned by extraction program. See :ref:`verification`
    for more details.  Can only take values table or yaml.  Default: table.
launch_parallel [string]
    Command template used to run the test program in parallel.  tc.nprocs is
    replaced with the number of processors a test uses (see run_cmd_template).
    If tc.nprocs does not appear, then testcode has no control over the number
    of processors a test is run on.  Default: mpirun -np tc.nprocs.
run_cmd_template [string]
    Template of command used to run the program on the test with the following
    substitutions made:

        tc.program
            replaced with the program to be tested.
        tc.args
            replaced with the arguments of the test.
        tc.input
            replaced with the input filename of the test.
        tc.output
            replaced with the filename for the standard output.  The filename
            is selected at runtime.
        tc.error
            replaced with the filename for the error output.  The filename is
            selected at runtime.
        tc.nprocs
            replaced with the number of processors the test is run on.

    Default: 'tc.program tc.args tc.input > tc.output 2> tc.error' in serial
    and 'launch_command tc.program tc.args tc.input > tc.output 2> tc.error' in
    parallel, where launch_command is specified above.  The parallel version is
    only used if the number of processors to run a test on is greater than
    zero.
skip_args [string]
    Arguments to supply to the program to test whether to skip the comparison
    of the test and benchmark.  Default: null string.
skip_cmd_template [string]
    Template of command used to test whether test was successfully run or
    whether the comparison of the benchmark and test output should be skipped.
    See :ref:`below <skip>` for more details.  The following strings in the
    template are replaced:

        tc.skip
            replaced with skip_program.
        tc.args
            replaced with skip_args.
        tc.test
            replaced with the filename of the test output.

    Default: tc.skip tc.args tc.test.
skip_program [string]
    Path to the program to test whether to skip the comparison of the test and
    benchmark.  If null, then this test is not performed.  Default: null string.
submit_pattern [string]
    String in the submit template to be replaced by the run command.  Default:
    testcode.run_cmd.
tolerance [tolerance format (see :ref:`below <tolerance>`.)]
    Default tolerance for tests of this type.  Default: inherits from
    [user].
verify [boolean]
    True if the extraction program compares the benchmark and test
    outputs directly.  See :ref:`verification` for more details.  Default:
    False.
vcs [string]
    Version control system used for the source code.  This is used to
    label the benchmarks.  The program binary is assumed to be in the same
    directory tree as the source code.  Supported values are: hg, git and svn
    and None.  If vcs is set to None, then the version id of the program is
    requested interactively when benchmarks are produced.  Default: None.

Most settings are optional and need only be set if certain functionality is
required or the default is not appropriate.  Note that either data_tag or
extract_program must be supplied.

In addition, the following variables are used, if present, as default settings
for all tests of this type:

* inputs_args (no default)
* nprocs (default: 0)
* min_nprocs (default: 0)
* max_nprocs (default: 2^31-1 or 2^63-1)
* output (no default)
* run_concurrent (defailt: false)
* submit_template
 
See :ref:`jobconfig` for more details.

All other settings are assumed to be paths to other versions of the program
(e.g. a stable version).  Using one of these versions instead of the one listed
under the 'exe' variable can be selected by an option to :ref:`testcode.py`.

.. _tolerance:

Tolerance format
----------------

The format for the tolerance for the data is very specific.  Individual
tolerance elements are specified in a comma-separated list.  Each individual
tolerance element is a python tuple (essentially a comma-separated list
enclosed in parentheses) consisting of, in order, the absolute tolerance, the
relative tolerance, the label of the field to which the tolerances apply and
a boolean value specifying the strictness of the tolerance (see below).  The
labels must be quoted.  If no label is supplied (or is set to None) then the
setting is taken to be the default tolerance to be applied to all data.  If the
strictness value is not given, the tolerance is assumed to be strict.  For
example, the setting::

    (1e-8, 1.e-6), (1.e-4, 1.e-4, 'Force')

uses an absolute tolerance of 10^-8 and a relative tolerance of 10^-6 by
default and an absolte tolerance and a relative tolerance of 10^-4 for data
items labelled with 'Force' (i.e. in columns headed by 'Force' using an
external data extraction program or labelled 'Force' by the internal data
extraction program using data tags).  If a tolerance is set to None, then it is
ignored.  At least one of the tolerances must be set.

A strict tolerance requires both the test value to be within the absolute and
relative tolerance of the benchmark value in order to be considered to pass.
This is the default behaviour.  A non-strict tolerance only requires the test
value to be within the absolute or relative tolerance of the benchmark value.
For example::

    (1e-8, 1e-6, None, False), (1e-10, 1e-10, 'Energy')

sets the default absolute and relative tolerances to be 10^-8 and 10^-6
respectively and sets the default tolerance to be non-strict except for the
'Energy' values, which have a strict absolute and relative tolerances of
10^-10.  If only one of the tolerances is set, then the strict and non-strict
settings are equivalent.

Alternatively, the tolerance can be labelled by a regular expression, in which case any
data labels which match the regular expression will use that tolerance unless there is
a tolerance with that specific label (i.e. exact matches override a regular
expression match).  Note that this is the case even if the tolerance using the exact
tolerance is defined in :ref:`userconfig` and the regular expression match is
defined in :ref:`jobconfig`.

.. _skip:

Skipping tests
--------------

Sometimes a test should not be compared to the benchmark---for example, if the
version of the program does not support a given feature or can only be run in
parallel.  testcode supports this by running a command to detect whether a test
should be skipped.

If the skipped program is set, then the skipped command is ran before
extracting data from output files.  For example, if

skip_program = grep
skip_args = "is not implemented."

are set, then testcode will run:

.. code-block:: bash

    grep "is not implemented." test_file

where test_file is the test output file.  If grep returns 0 (i.e.
test_file contains the string "is not implemented") then the test is
marked as skipped and the test file is not compared to the benchmark.
