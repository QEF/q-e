#!/usr/bin/env python
'''testcode [options] [action1 [action2...]]

testcode is a simple framework for comparing output from (principally numeric)
programs to previous output to reveal regression errors or miscompilation.

Run a set of actions on a set of tests.

Available actions:
  compare               compare set of test outputs from a previous testcode
                        run against the benchmark outputs.
  diff                  diff set of test outputs from a previous testcode
                        run against the benchmark outputs.
  make-benchmarks       create a new set of benchmarks and update the userconfig
                        file with the new benchmark id.  Also runs the 'run'
                        action unless the 'compare' action is also given.
  recheck               compare a set of test outputs and rerun failed tests.
  run                   run a set of tests and compare against the benchmark
                        outputs.  Default action.
  tidy                  Remove files from previous testcode runs from the test
                        directories.

Requires two configuration files, jobconfig and userconfig.  See testcode
documentation for further details.'''

# copyright: (c) 2012 James Spencer
# license: modified BSD; see LICENSE for more details

import glob
import optparse
import os
import re
import subprocess
import sys
import threading
import time

try:
    import testcode2
except ImportError:
    # try to find testcode2 assuming it is being run directly from the source
    # layout.
    SCRIPT_DIR = os.path.abspath(os.path.dirname(sys.argv[0]))
    TESTCODE2_LIB = os.path.join(SCRIPT_DIR, '../lib/')
    sys.path.extend([TESTCODE2_LIB])
    import testcode2

import testcode2.config
import testcode2.util
import testcode2.compatibility
import testcode2.exceptions
import testcode2.validation

#--- testcode initialisation ---

def init_tests(userconfig, jobconfig, test_id, reuse_id, executables=None,
        categories=None, nprocs=-1, benchmark=None, userconfig_options=None,
        jobconfig_options=None):
    '''Initialise tests from the configuration files and command-line options.

userconfig, executables, test_id and userconfig_options are passed to
testcode2.config.userconfig.

jobconfig and jobconfig_options are passed to testcode2.config.parse_jobconfig.

categories is passed to testcode2.config.select_tests.

test_id is used to set the test identifier.  If test_id is null and reused_id
is true, then the identifier is set to that of the last tests ran by testcode
otherwise a unique identifier based upon the date is used.

nprocs is the number of processors each test is run on.  If negative, the
defaults in the configuration files are used.

benchmark is the benchmark id labelling the set of benchmarks to compare the
tests too.  If None, the default in userconfig is used.

Returns:

user_options: dictionary containing user options specified in userconfig.
test_programs: dict of the test programs defined in userconfig.
tests: list of selected tests.
'''

    config_exists = os.path.exists(userconfig) and os.path.exists(jobconfig)

    try:
        (user_options, test_programs) = testcode2.config.parse_userconfig(
                userconfig, executables, test_id, userconfig_options)
    except testcode2.exceptions.TestCodeError:
        err = str(sys.exc_info()[1])
        if not config_exists:
            err += (' Please run from a directory containing (or specify) the '
                    'userconfig file. Use ``--help`` to see available options.')
        raise testcode2.exceptions.TestCodeError(err)

    # Set benchmark if required.
    if benchmark:
        for key in test_programs:
            test_programs[key].benchmark = [benchmark]

    try:
        (tests, test_categories) = testcode2.config.parse_jobconfig(
                jobconfig, user_options, test_programs, jobconfig_options)
    except testcode2.exceptions.TestCodeError:
        err = str(sys.exc_info()[1])
        if not config_exists:
            err += (' Please run from a directory containing (or specify) the '
                    'jobconfig file. Use ``--help`` to see available options.')
        raise testcode2.exceptions.TestCodeError(err)

    # Set number of processors...
    if nprocs >= 0:
        for test in tests:
            test.nprocs = nprocs
            if test.nprocs < test.min_nprocs:
                test.nprocs = test.min_nprocs
            if test.nprocs > test.max_nprocs:
                test.nprocs = test.max_nprocs

    # parse selected job categories from command line
    # Remove those tests which weren't run most recently if comparing.
    if categories:
        tests = testcode2.config.select_tests(tests, test_categories,
                categories, os.path.abspath(os.path.dirname(userconfig)))

    # Sort by path (as that's how they appear in the user's directory).
    tests.sort(key=lambda test: test.path)

    if not test_id:
        test_id = testcode2.config.get_unique_test_id(tests, reuse_id,
                user_options['date_fmt'])
        for key in test_programs:
            test_programs[key].test_id = test_id

    return (user_options, test_programs, tests)

#--- create command line interface ---

def parse_cmdline_args(args):
    '''Parse command line arguments.

args: list of supplied arguments.

Returns:

options: object returned by optparse containing the options.
actions: list of testcode2 actions to run.
'''

    # Curse not being able to use argparse in order to support python <= 2.7!
    parser = optparse.OptionParser(usage=__doc__)

    allowed_actions = ['compare', 'run', 'diff', 'tidy', 'make-benchmarks',
                       'recheck']

    parser.add_option('-b', '--benchmark', help='Set the file ID of the '
            'benchmark files.  Default: specified in the [user] section of the '
            'userconfig file.')
    parser.add_option('-c', '--category', action='append', default=[],
            help='Select the category/group of tests.  Can be specified '
            'multiple times.  Default: use the _default_ category if run is an '
            'action unless make-benchmarks is an action.  All other cases use '
            'the _all_ category by default.  The _default_ category contains '
            'all  tests unless otherwise set in the jobconfig file.')
    parser.add_option('-e', '--executable', action='append', default=[],
            help='Set the executable(s) to be used to run the tests.  Can be'
            ' a path or name of an option in the userconfig file, in which'
            ' case all test programs are set to use that value, or in the'
            ' format program_name=value, which affects only the specified'
            ' program.')
    parser.add_option('-i', '--insert', action='store_true', default=False,
            help='Insert the new benchmark into the existing list of benchmarks'
            ' in userconfig rather than overwriting it.  Only relevant to the'
            ' make-benchmarks action.  Default: %default.')
    parser.add_option('--jobconfig', default='jobconfig', help='Set path to the'
            ' job configuration file.  Default: %default.')
    parser.add_option('--job-option', action='append', dest='job_option',
            default=[], nargs=3, help='Override/add setting to jobconfig.  '
            'Takes three arguments.  Format: section_name option_name value.  '
            'Default: none.')
    parser.add_option('--older-than', type='int', dest='older_than', default=14,
            help='Set the age (in days) of files to remove.  Only relevant to '
            'the tidy action.  Default: %default days.')
    parser.add_option('-p', '--processors', type='int', default=-1,
            dest='nprocs', help='Set the number of processors to run each test '
            'on.  Default: use settings in configuration files.')
    parser.add_option('-q', '--quiet', action='store_const', const=0, 
            dest='verbose', default=1, help='Print only minimal output.  '
            'Default: False.')
    parser.add_option('-s', '--submit', dest='queue_system', default=None,
            help='Submit tests to a queueing system of the specified type.  '
            'Only PBS system is currently implemented.  Default: %default.')
    parser.add_option('-t', '--test-id', dest='test_id', help='Set the file ID '
            'of the test outputs.  Default: unique filename based upon date '
            'if running tests and most recent test_id if comparing tests.')
    parser.add_option('--total-processors', type='int', default=-1,
            dest='tot_nprocs', help='Set the total number of processors to use '
            'to run tests concurrently.  Relevant only to the run option.  '
            'Default: run all tests concurrently run if --submit is used; run '
            'tests sequentially otherwise.')
    parser.add_option('--userconfig', default='userconfig', help='Set path to '
            'the user configuration file.  Default: %default.')
    parser.add_option('--user-option', action='append', dest='user_option',
            default=[], nargs=3, help='Override/add setting to userconfig.  '
            'Takes three arguments.  Format: section_name option_name value.  '
            'Default: none.')
    parser.add_option('-v', '--verbose', default=1, action="count", 
            dest='verbose', help='Increase verbosity of output.  Can be '
            'specified multiple times.')

    (options, args) = parser.parse_args(args)

    # Default action.
    if not args or ('make-benchmarks' in args and 'compare' not in args
            and 'run' not in args):
        # Run tests by default if no action provided.
        # Run tests before creating benchmark by default.
        args.append('run')

    # Default category.
    if not options.category:
        # We quietly filter out tests which weren't run last when diffing
        # or comparing.
        options.category = ['_all_']
        if 'run' in args and 'make-benchmarks' not in args:
            options.category = ['_default_']

    test_args = (arg not in allowed_actions for arg in args)
    if testcode2.compatibility.compat_any(test_args):
        print('At least one action is not understood: %s.' % (' '.join(args)))
        parser.print_usage()
        sys.exit(1)

    # Parse executable option to form dictionary in format expected by
    # parse_userconfig.
    exe = {}
    for item in options.executable:
        words = item.split('=')
        if len(words) == 1:
            # setting executable for all programs (unless otherwise specified)
            exe['_tc_all'] = words[0]
        else:
            # format: program_name=executable
            exe[words[0]] = words[1]
    options.executable = exe

    # Set FILESTEM if test_id refers to a benchmark file or the benchmark
    # refers to a test_id.
    filestem = testcode2.FILESTEM.copy()
    if options.benchmark and options.benchmark[:2] == 't:':
        filestem['benchmark'] = testcode2.FILESTEM['test']
        options.benchmark = options.benchmark[2:]
    if options.test_id and options.test_id[:2] == 'b:':
        filestem['test'] = testcode2.FILESTEM['benchmark']
        options.test_id = options.test_id[2:]
    if filestem['test'] != testcode2.FILESTEM['test'] and 'run' in args:
        print('Not allowed to set test filename to be a benchmark filename '
                'when running calculations.')
        sys.exit(1)
    testcode2.FILESTEM = filestem.copy()

    # Convert job-options and user-options to dict of dicsts format.
    for item in ['user_option', 'job_option']:
        uj_opt = getattr(options, item)
        opt = dict( (section, {}) for section in
                testcode2.compatibility.compat_set(opt[0] for opt in uj_opt) )
        for (section, option, value) in uj_opt:
            opt[section][option] = value
        setattr(options, item, opt)

    return (options, args)

#--- actions ---

def run_tests(tests, verbose=1, cluster_queue=None, tot_nprocs=0):
    '''Run tests.

tests: list of tests.
verbose: level of verbosity in output.
cluster_queue: name of cluster system to use.  If None, tests are run locally.
    Currently only PBS is implemented.
tot_nprocs: total number of processors available to run tests on.  As many
    tests (in a LIFO fashion from the tests list) are run at the same time as
    possible without using more processors than this value.  If less than 1 and
    cluster_queue is specified, then all tests are submitted to the cluster at
    the same time.  If less than one and cluster_queue is not set, then
    tot_nprocs is ignored and the tests are run sequentially (default).
'''
    def run_test_worker(semaphore, semaphore_lock, tests, *run_test_args):
        '''Launch a test after waiting until resources are available to run it.

semaphore: threading.Semaphore object containing the number of cores/processors
    which can be used concurrently to run tests.
semaphore.lock: threading.Lock object used to restrict acquiring the semaphore
    to one thread at a time.
tests: list of (serialized) tests to run in this thread.
run_test_args: arguments to pass to test.run_test method.
'''

        # Ensure that only one test attempts to register resources with the
        # semaphore at a time.  This restricts running the tests to a LIFO
        # fashion which is not perfect (we don't attempt to backfill with
        # smaller tests, for example) but is a reasonable and (most
        # importantly) simple first-order approach.
        for test in tests:
            semaphore_lock.acquire()
            # test.nprocs is <1 when program is run in serial.
            nprocs_used = max(1, test.nprocs)
            for i in range(nprocs_used):
                semaphore.acquire()
            semaphore_lock.release()

            test.run_test(*run_test_args)

            for i in range(nprocs_used):
                semaphore.release()

    # Check executables actually exist...
    compat = testcode2.compatibility
    executables = [test.test_program.exe for test in tests]
    executables = compat.compat_set(executables)
    for exe in executables:
        mswin = sys.platform.startswith('win') or sys.platform.startswith('cyg')
        # The test is not reliable if there's an unholy combination of windows
        # and cygwin being used to run the program.  We've already warned the
        # user (in config.set_program_name) that we struggled to find the
        # executable.
        if not os.path.exists(exe) and not mswin:
            err = 'Executable does not exist: %s.' % (exe)
            raise testcode2.exceptions.TestCodeError(err)

    if tot_nprocs <= 0 and cluster_queue:
        # Running on cluster.  Default to submitting all tests at once.
        tot_nprocs = sum(test.nprocs for test in tests)

    if tot_nprocs > 0:
        # Allow at most tot_nprocs cores to be used at once by tests.
        max_test_nprocs = max(test.nprocs for test in tests)
        if max_test_nprocs > tot_nprocs:
            err = ('Number of available cores less than the number required by '
                   'the largest test: at least %d needed, %d available.'
                   % (max_test_nprocs, tot_nprocs))
            raise testcode2.exceptions.TestCodeError(err)

        # Need to serialize tests that run in the same directory with wildcard
        # patterns in the output file--otherwise we can't figure out which
        # output file belongs to which test.  We might be able to for some
        # wildcards, but let's err on the side of caution.
        wildcards = re.compile('.*(\*|\?|\[.*\]).*')
        serialized_tests = []
        test_store = {}
        for test in tests:
            if test.output and wildcards.match(test.output):
                if test.path in test_store:
                    test_store[test.path].append(test)
                else:
                    test_store[test.path] = [test]
            else:
                serialized_tests.append([test])
        for (key, stests) in test_store.items():
            if (len(stests) > 1) and verbose > 2:
                print('Warning: cannot run tests in %s concurrently.' % stests[0].path)
        serialized_tests += test_store.values()

        semaphore = threading.BoundedSemaphore(tot_nprocs)
        slock = threading.Lock()
        jobs = [threading.Thread(
                    target=run_test_worker,
                    args=(semaphore, slock, test, verbose, cluster_queue,
                          os.getcwd())
                                )
                    for test in serialized_tests]
        for job in jobs:
            # daemonise so thread terminates when master dies
            try:
                job.setDaemon(True)
            except AttributeError:
                job.daemon = True  
            job.start()

        # We avoid .join() which is blocking making it unresponsive to TERM
        while threading.activeCount() > 1:
            time.sleep(0.5)
    else:
        # run straight through, one at a time
        for test in tests:
            test.run_test(verbose, cluster_queue, os.getcwd())


def compare_tests(tests, verbose=1):
    '''Compare tests.

tests: list of tests.
verbose: level of verbosity in output.

Returns:

number of tests not checked due to test output file not existing.
'''

    not_checked = 0

    for test in tests:
        for (inp, args) in test.inputs_args:
            test_file = testcode2.util.testcode_filename(
                    testcode2.FILESTEM['test'],
                    test.test_program.test_id, inp, args
                    )
            test_file = os.path.join(test.path, test_file)
            if os.path.exists(test_file):
                test.verify_job(inp, args, verbose, os.getcwd())
            else:
                if verbose > 0 and verbose <= 2:
                    info_line = testcode2.util.info_line(test.path, inp, args, os.getcwd())
                    print('%sNot checked.' % info_line)
                if verbose > 1:
                    print('Skipping comparison.  '
                          'Test file does not exist: %s.\n' % test_file)
                not_checked += 1

    return not_checked

def recheck_tests(tests, verbose=1, cluster_queue=None, tot_nprocs=0):
    '''Check tests and re-run any failed/skipped tests.

tests: list of tests.
verbose: level of verbosity in output.
cluster_queue: name of cluster system to use.  If None, tests are run locally.
    Currently only PBS is implemented.
tot_nprocs: total number of processors available to run tests on.  As many
    tests (in a LIFO fashion from the tests list) are run at the same time as
    possible without using more processors than this value.  If less than 1 and
    cluster_queue is specified, then all tests are submitted to the cluster at
    the same time.  If less than one and cluster_queue is not set, then
    tot_nprocs is ignored and the tests are run sequentially (default).

Returns:

not_checked: number of tests not checked due to missing test output.
'''

    if verbose == 0:
        sep = ' '
    else:
        sep = '\n\n'

    sys.stdout.write('Comparing tests to benchmarks:'+sep)

    not_checked = compare_tests(tests, verbose)
    end_status(tests, not_checked, verbose, False)

    rerun_tests = []
    skip = testcode2.validation.Status(name='skipped')
    for test in tests:
        stat = test.get_status()
        if sum(stat[key] for key in ('failed', 'unknown')) != 0:
            rerun_tests.append(test)
        elif stat['ran'] != 0:
            # mark tests as skipped using an internal API (naughty!)
            for inp_arg in test.inputs_args:
                test._update_status(skip, inp_arg)

    if verbose > 0:
        print('')
    if rerun_tests:
        sys.stdout.write('Rerunning failed tests:'+sep)
        run_tests(rerun_tests, verbose, cluster_queue, tot_nprocs)

    return not_checked

def diff_tests(tests, diff_program, verbose=1):
    '''Diff tests.

tests: list of tests.
diff_program: diff program to use.
verbose: level of verbosity in output.
'''

    for test in tests:
        cwd = os.getcwd()
        os.chdir(test.path)
        for (inp, args) in test.inputs_args:
            have_benchmark = True
            try:
                benchmark = test.test_program.select_benchmark_file(
                        test.path, inp, args
                        )
            except testcode2.exceptions.TestCodeError:
                err = sys.exc_info()[1]
                have_benchmark = False
            test_file = testcode2.util.testcode_filename(
                    testcode2.FILESTEM['test'],
                    test.test_program.test_id, inp, args
                    )
            if not os.path.exists(test_file):
                if verbose > 0:
                    print('Skipping diff with %s in %s: %s does not exist.'
                            % (benchmark, test.path, test_file))
            elif not have_benchmark:
                if verbose > 0:
                    print('Skipping diff with %s. %s' % (test.path, err))
            else:
                if verbose > 0:
                    print('Diffing %s and %s in %s.' %
                            (benchmark, test_file, test.path))
                diff_cmd = '%s %s %s' % (diff_program, benchmark, test_file)
                diff_popen = subprocess.Popen(diff_cmd, shell=True)
                diff_popen.wait()
        os.chdir(cwd)

def tidy_tests(tests, ndays):
    '''Tidy up test directories.

tests: list of tests.
ndays: test files older than ndays are deleted.
'''

    epoch_time = time.time() - 86400*ndays

    test_globs = ['test.out*','test.err*']

    print(
            'Delete all %s files older than %s days from each job directory?'
                % (' '.join(test_globs), ndays)
         )
    ans = ''
    while ans != 'y' and ans != 'n':
        ans = testcode2.compatibility.compat_input('Confirm [y/n]: ')

    if ans == 'n':
        print('No files deleted.')
    else:
        for test in tests:
            cwd = os.getcwd()
            os.chdir(test.path)
            if test.submit_template:
                file_globs = test_globs + [test.submit_template]
            else:
                file_globs = test_globs
            for file_glob in file_globs:
                for test_file in glob.glob(file_glob):
                    if os.stat(test_file)[-2] < epoch_time:
                        os.remove(test_file)
            os.chdir(cwd)

def make_benchmarks(test_programs, tests, userconfig, copy_files_since,
        insert_id=False):
    '''Make a new set of benchmarks.

test_programs: dictionary of test programs.
tests: list of tests.
userconfig: path to the userconfig file.  This is updated with the new benchmark id.
copy_files_since: files produced since the timestamp (in seconds since the
    epoch) are copied to the testcode_data subdirectory in each test.
insert_id: insert the new benchmark id into the existing list of benchmark ids in
    userconfig if True, otherwise overwrite the existing benchmark ids with the
    new benchmark id (default).
'''

    # All tests passed?
    statuses = [test.get_status() for test in tests]
    npassed = sum(status['passed'] for status in statuses)
    nran = sum(status['ran'] for status in statuses)
    if npassed != nran:
        ans = ''
        print('Not all tests passed.')
        while ans != 'y' and ans != 'n':
            ans = testcode2.compatibility.compat_input(
                                                'Create new benchmarks? [y/n] ')
        if ans != 'y':
            return None

    # Get vcs info.
#    vcs = {}
#    for (key, program) in test_programs.items():
#        if program.vcs and program.vcs.vcs:
#            vcs[key] = program.vcs.get_code_id()
#        else:
#            print('Program not under (known) version control system')
#            vcs[key] = testcode2.compatibility.compat_input(
#                    'Enter revision id for %s: ' % (key))

    # HACK
    code_id = testcode2.compatibility.compat_input(
                    'Enter new revision id : ')

    # Benchmark label from vcs info.
#    if len(vcs) == 1:
#        benchmark = vcs.popitem()[1]
#    else:
#        benchmark = []
#        for (key, code_id) in vcs.items():
#            benchmark.append('%s-%s' % (key, code_id))
#        benchmark = '.'.join(benchmark)

    # HACK
    benchmark = code_id

    # Create benchmarks.
    for test in tests:
        test.create_new_benchmarks(benchmark, copy_files_since)

    # update userconfig file.
    if userconfig:
        config = testcode2.compatibility.configparser.RawConfigParser()
        config.optionxform = str # Case sensitive file.
        config.read(userconfig)
        if insert_id:
            ids = config.get('user', 'benchmark').split()
            if benchmark in ids:
                ids.remove(benchmark)
            ids.insert(0, benchmark)
            benchmark = ' '.join(ids)
        if len(benchmark.split()) > 1:
            print('Setting new benchmarks in userconfig to be: %s.' %
                    (benchmark))
        else:
            print('Setting new benchmark in userconfig to be: %s.' %
                    (benchmark))
        config.set('user', 'benchmark', benchmark)
        userconfig = open(userconfig, 'w')
        config.write(userconfig)
        userconfig.close()

#--- info output ---

def start_status(tests, running, verbose=1):
    '''Print a header containing useful information.

tests: list of tests.
running: true if tests are to be run.
verbose: level of verbosity in output (no output if <1).
'''

    if verbose > 0:
        exes = [test.test_program.exe for test in tests]
        exes = testcode2.compatibility.compat_set(exes)
        if running:
            for exe in exes:
                print('Using executable: %s.' % (exe))
        # All tests use the same test_id and benchmark.
        print('Test id: %s.' % (tests[0].test_program.test_id))
        if len(tests[0].test_program.benchmark) > 1:
            benchmark_ids = ', '.join(tests[0].test_program.benchmark)
            print('Benchmarks: %s.' % (benchmark_ids))
        else:
            print('Benchmark: %s.' % (tests[0].test_program.benchmark[0]))
        print('')

def end_status(tests, not_checked=0, verbose=1, final=True):
    '''Print a footer containing useful information.

tests: list of tests.
not_checked: number of tests not checked (ie not run or compared).
verbose: level of verbosity in output.  A summary footer is produced if greater
    than 0; otherwise a minimal status line is printed out.
final: final call (so print a goodbye messge).
'''

    def pluralise(string, num):
        '''Return plural form (just by adding s) to string if num > 1.'''
        if num > 1:
            string = string+'s'
        return string

    def select_tests(stat_key, tests, statuses):
        '''Select a subset of tests.

        (test.name, test.path) is included if the test object contains at least
        one test of the desired status (stat_key).'''
        test_subset = [(test.name, test.path) for (test, status)
                               in zip(tests, statuses) if status[stat_key] != 0]
        return sorted(test_subset)

    def format_test_subset(subset):
        '''Format each entry in the list returned by select_tests.'''
        subset_fmt = []
        for (name, path) in subset:
            if os.path.abspath(name) == os.path.abspath(path):
                entry = name
            else:
                entry = '%s (test name: %s)' % (path, name)
            if entry not in subset_fmt:
                subset_fmt.append(entry)
        return subset_fmt

    statuses = [test.get_status() for test in tests]
    npassed = sum(status['passed'] for status in statuses)
    nwarning = sum(status['warning'] for status in statuses)
    nfailed = sum(status['failed'] for status in statuses)
    nunknown = sum(status['unknown'] for status in statuses)
    nskipped = sum(status['skipped'] for status in statuses)
    nran = sum(status['ran'] for status in statuses)
    failures = format_test_subset(select_tests('failed', tests, statuses))
    warnings = format_test_subset(select_tests('warning', tests, statuses))
    skipped = format_test_subset(select_tests('skipped', tests, statuses))
    # Treat warnings as passes but add a note about how many warnings.
    npassed += nwarning
    # Treat skipped tests as tests which weren't run.
    nran -= nskipped

    # Pedantic.
    warning = pluralise('warning', nwarning)
    ran_test = pluralise('test', nran)
    failed_test = pluralise('test', nfailed)
    skipped_test = pluralise('test', nskipped)

    add_info_msg = []
    if nwarning != 0:
        add_info_msg.append('%s %s' % (nwarning, warning))
    if nskipped != 0:
        add_info_msg.append('%s skipped' % (nskipped,))
    if nunknown != 0:
        add_info_msg.append('%s unknown' % (nunknown,))
    if not_checked != 0:
        add_info_msg.append('%s not checked' % (not_checked,))
    add_info_msg = ', '.join(add_info_msg)
    if add_info_msg:
        add_info_msg = ' (%s)' % (add_info_msg,)

    if nran == 0:
        print('No tests to run.')
    elif verbose > 0:
        if verbose < 2:
            print('') # Obsessive formatting.
        msg = '%s%s out of %s %s passed%s.'
        if final:
            msg = 'All done. %s' % (msg,)
        if npassed == nran:
            print(msg % ('', npassed, nran, ran_test, add_info_msg))
        else:
            print(msg % ('ERROR: only ', npassed, nran, ran_test, add_info_msg))
        if failures:
            print('Failed %s in:\n\t%s' % (failed_test, '\n\t'.join(failures)))
        if warnings:
            print('%s in:\n\t%s' % (warning.title(), '\n\t'.join(warnings)))
        if skipped:
            print('Skipped %s in:\n\t%s' % (skipped_test, '\n\t'.join(skipped)))
    else:
        print(' [%s/%s%s]'% (npassed, nran, add_info_msg))

    # ternary operator not in python 2.4. :-(
    ret_val = 0
    if nran != npassed:
        ret_val = 1

    return ret_val

#--- main runner ---

def main(args):
    '''main controller procedure.

args: command-line arguments passed to testcode2.
'''

    start_time = time.time()

    (options, actions) = parse_cmdline_args(args)

    # Shortcut names to options used multiple times.
    verbose = options.verbose
    userconfig = options.userconfig
    reuse_id = 'run' not in actions and testcode2.compatibility.compat_any(
                [action in actions for action in ['compare', 'diff', 'recheck']]
                )

    (user_options, test_programs, tests) = init_tests(userconfig,
            options.jobconfig, options.test_id, reuse_id,
            options.executable, options.category, options.nprocs,
            options.benchmark, options.user_option,
            options.job_option)

    ret_val = 0
    if not (len(actions) == 1 and 'tidy' in actions):
        start_status(tests, 'run' in actions, verbose)
    if 'run' in actions:
        run_tests(tests, verbose, options.queue_system, options.tot_nprocs)
        ret_val = end_status(tests, 0, verbose)
    if 'recheck' in actions:
        not_checked = recheck_tests(tests, verbose,
                                        options.queue_system,options.tot_nprocs)
        ret_val = end_status(tests, not_checked, verbose)
    if 'compare' in actions:
        not_checked = compare_tests(tests, verbose)
        ret_val = end_status(tests, not_checked, verbose)
    if 'diff' in actions:
        diff_tests(tests, user_options['diff'], verbose)
    if 'tidy' in actions:
        tidy_tests(tests, options.older_than)
    if 'make-benchmarks' in actions:
        make_benchmarks(test_programs, tests, userconfig, start_time,
                options.insert)

    return ret_val

if __name__ == '__main__':

    try:
        sys.exit(main(sys.argv[1:]))
    except testcode2.exceptions.TestCodeError:
        err = sys.exc_info()[1]
        print(err)
        sys.exit(1)
