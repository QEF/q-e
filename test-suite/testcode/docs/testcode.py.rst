.. _testcode.py:

testcode.py
===========

.. only:: html

    testcode.py - a command-line interface to testcode.

Synopsis
--------

testcode.py [options] [action1 [action2...]]

Description
-----------

Run a set of actions on a set of tests.

Requires two configuration files, :ref:`jobconfig` and :ref:`userconfig`.  See
testcode documentation for further details.

testcode.py provides a command-line interface to testcode, a simple framework
for comparing output from (principally numeric) programs to previous output to
reveal regression errors or miscompilation.

Actions
-------

''run'' is th default action.

compare
    compare set of test outputs from a previous testcode run against the
    benchmark outputs.
diff
    diff set of test outputs from a previous testcode run against the benchmark
    outputs.
make-benchmarks
    create a new set of benchmarks and update the :ref:`userconfig` file with
    the new benchmark id.  Also runs the 'run' action unless the 'compare'
    action is also given.
recheck
    compare set of test outputs from a previous testcode run against
    benchmark outputs and rerun any failed tests.
run
    run a set of tests and compare against the benchmark outputs.
tidy
    Remove files from previous testcode runs from the test directories.

Options
-------

-h, --help
    show this help message and exit
-b BENCHMARK, --benchmark=BENCHMARK
    Set the file ID of the benchmark files.  If BENCHMARK is in the format
    t:ID, then the test files with the corresponding ID are used.  This
    allows two sets of tests to be compared.  Default: specified in the [user]
    section of the :ref:`userconfig` file.
-c CATEGORY, --category=CATEGORY
    Select the category/group of tests.  Can be specified multiple times.
    Wildcards or parent directories can be used to select multiple directories
    by their path.  Default: use the `_default_` category if run is an action
    unless make-benchmarks is an action.  All other cases use the `_all_`
    category by default.  The `_default_` category contains all  tests unless
    otherwise set in the :ref:`jobconfig` file.
-e EXECUTABLE, --executable=EXECUTABLE
    Set the executable(s) to be used to run the tests.  Can be  a path or name
    of an option in the :ref:`userconfig` file, in which case all test programs are
    set to use that value, or in the format program_name=value, which affects
    only the specified program.  Only relevant to the run action.  Default: exe
    variable set for each program listed in the :ref:`userconfig` file.
-i, --insert
    Insert the new benchmark into the existing list of benchmarks in userconfig
    rather than overwriting it.  Only relevant to the make-benchmarks action.
    Default: False.
--jobconfig=JOBCONFIG
    Set path to the job configuration file.  Default: jobconfig.
--job-option=JOB_OPTION
    Override/add setting to :ref:`jobconfig`.  Takes three arguments.  Format:
    section_name option_name value.  Default: none.
--older-than=OLDER_THAN
    Set the age (in days) of files to remove.  Only relevant to the tidy
    action.  Default: 14 days.
-p NPROCS, --processors=NPROCS
    Set the number of processors to run each test on.  Only relevant to the run
    action.  Default: run tests as serial jobs.
-q, --quiet
    Print only minimal output.  Default: False.
-s QUEUE_SYSTEM, --submit=QUEUE_SYSTEM
    Submit tests to a queueing system of the specified type.  Only PBS system
    is currently implemented.  Only relevant to the run action.  Default: none.
-t TEST_ID, --test-id=TEST_ID
    Set the file ID of the test outputs.  If TEST_ID is in the format b:ID, then
    the benchmark files with the corresponding ID are used.  This allows two
    sets of benchmarks to be compared.  Default: unique filename based upon
    date if running tests and most recent test_id if comparing tests.
--total-processors=TOT_NPROCS
    Set the total number of processors to use to run as many tests as possible
    at the same time.  Relevant only to the run option.  Default: run all tests
    concurrently run if --submit is used; run tests sequentially otherwise.
--userconfig=USERCONFIG
    Set path to the user configuration file.  Default: userconfig.
--user-option=USER_OPTION
    Override/add setting to :ref:`userconfig`.  Takes three arguments.  Format:
    section_name option_name value.  Default: none.
-v, --verbose
    Increase verbosity of output.  Can be specified up to two times.
    The default behaviour is to print out the test and its status.  (See the
    --quiet option to suppress even this.)  Specify -v or --verbose once to
    show (if relevant) which data values caused warnings or failures.
    Specify -v or --verbose twice to see all (external) commands run and all
    data extracted from running the tests.  Using the maximum verbosity level
    is highly recommended for debugging.

Exit status
-----------

1 if one or more tests fail (run and compare actions only) and 0 otherwise.

License
-------

Modified BSD License.  See LICENSE in the source code for more details.

Bugs
----

Contact James Spencer (j.spencer@imperial.ac.uk) regarding bug reports,
suggestions for improvements or code contributions.
