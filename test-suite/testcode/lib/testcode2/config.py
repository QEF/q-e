'''
testcode2.config
----------------

Parse jobconfig and userconfig ini files.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import copy
import glob
import os
import subprocess
import time
import warnings

import testcode2
import testcode2.compatibility as compat
import testcode2.exceptions as exceptions
import testcode2.util as util
import testcode2.validation as validation
import testcode2.vcs as vcs

def eval_nested_tuple(string):
    nested_tuple = compat.literal_eval(string)
    if isinstance(nested_tuple[0], (list, tuple)):
        return nested_tuple
    else:
        # Append a comma to the option to ensure literal_eval returns a tuple
        # of tuples, even if the option only contains a single tuple.
        return compat.literal_eval('%s,' % string)

def parse_tolerance_tuple(val):
    '''Parse (abs_tol,rel_tol,name,strict).'''
    if len(val) >= 4:
        strict = val[3]
    else:
        strict = True
    if len(val) >= 3:
        name = val[2]
    else:
        name = None
    if len(val) >= 2:
        rel_tol = val[1]
    else:
        rel_tol = None
    if len(val) >= 1:
        abs_tol = val[0]
    else:
        abs_tol = None
    return (name, validation.Tolerance(name, abs_tol, rel_tol, strict))

def parse_userconfig(config_file, executables=None, test_id=None,
        settings=None):
    '''Parse the user options and job types from the userconfig file.

config_file: location of the userconfig file, either relative or absolute.'''

    if executables is None:
        executables = {}

    if not os.path.exists(config_file):
        raise exceptions.TestCodeError(
                'User configuration file %s does not exist.' % (config_file)
                                      )
    # paths to programs can be specified relative to the config
    # file.
    config_directory = os.path.dirname(os.path.abspath(config_file))

    userconfig = compat.configparser.RawConfigParser()
    userconfig.optionxform = str # Case sensitive file.
    userconfig.read(config_file)

    # Alter config file with additional settings provided.
    if settings:
        for (section_key, section) in settings.items():
            for (option_key, value) in section.items():
                userconfig.set(section_key, option_key, value)

    # Sensible defaults for the user options.
    user_options = dict(benchmark=None, date_fmt='%d%m%Y',
            tolerance='(1.e-10,None)', output_files=None, diff='diff')

    if userconfig.has_section('user'):
        user_options.update(dict(userconfig.items('user')))
        userconfig.remove_section('user')
        user_options['tolerance'] = dict(
                (parse_tolerance_tuple(item)
                     for item in eval_nested_tuple(user_options['tolerance']))
                                        )
        if user_options['benchmark']:
            user_options['benchmark'] = user_options['benchmark'].split()
    else:
        raise exceptions.TestCodeError(
                'user section in userconfig does not exist.'
                                      )

    if not userconfig.sections():
        raise exceptions.TestCodeError(
                'No job types specified in userconfig.'
                                      )

    test_program_options = ('run_cmd_template',
        'launch_parallel', 'ignore_fields', 'data_tag', 'extract_cmd_template',
        'extract_program', 'extract_args', 'extract_fmt', 'verify', 'vcs',
        'skip_program', 'skip_args', 'skip_cmd_template')
    default_test_options = ('inputs_args', 'output', 'nprocs',
        'min_nprocs', 'max_nprocs', 'submit_template',)
    test_programs = {}
    for section in userconfig.sections():
        tp_dict = {}
        tolerances = copy.deepcopy(user_options['tolerance'])
        # Read in possible TestProgram settings.
        for item in test_program_options:
            if userconfig.has_option(section, item):
                tp_dict[item] = userconfig.get(section, item)
        if 'ignore_fields' in tp_dict:
            tp_dict['ignore_fields'] = tp_dict['ignore_fields'].split()
        if section in executables:
            exe = executables[section]
        elif '_tc_all' in executables:
            exe = executables['_tc_all']
        else:
            exe = 'exe'
        if userconfig.has_option(section, exe):
            # exe is set to be a key rather than the path to an executable.
            # Expand.
            exe = userconfig.get(section, exe)
        # Create a default test settings.
        # First, tolerances...
        if userconfig.has_option(section, 'tolerance'):
            for item in (
                    eval_nested_tuple(userconfig.get(section, 'tolerance'))
                        ):
                (name, tol) = parse_tolerance_tuple(item)
                tolerances[name] = tol
        test_dict = dict(
                         default_tolerance=tolerances[None],
                         tolerances=tolerances,
                        )
        # Other settings...
        for item in default_test_options:
            if userconfig.has_option(section, item):
                test_dict[item] = userconfig.get(section, item)
        if userconfig.has_option(section, 'run_concurrent'):
            test_dict['run_concurrent'] = \
                    userconfig.getboolean(section, 'run_concurrent')
        # Programs can be specified relative to the config directory.
        exe = set_program_name(exe, config_directory)
        if 'extract_program' in tp_dict:
            tp_dict['extract_program'] = set_program_name(
                                tp_dict['extract_program'], config_directory)
        if 'skip_program' in tp_dict:
            tp_dict['skip_program'] = set_program_name(
                                tp_dict['skip_program'], config_directory)
        if 'submit_template' in test_dict:
            test_dict['submit_template'] = os.path.join(config_directory,
                                                   test_dict['submit_template'])
        for key in ('nprocs', 'max_nprocs', 'min_nprocs'):
            if key in test_dict:
                test_dict[key] = int(test_dict[key])
        if 'inputs_args' in test_dict:
            # format: (input, arg), (input, arg)'
            test_dict['inputs_args'] = (
                    eval_nested_tuple(test_dict['inputs_args']))
        # Create a default test.
        tp_dict['default_test_settings'] = testcode2.Test(None, None, None,
                **test_dict)
        if 'vcs' in tp_dict:
            tp_dict['vcs'] = vcs.VCSRepository(tp_dict['vcs'],
                    os.path.dirname(exe))
        program = testcode2.TestProgram(section, exe, test_id,
            user_options['benchmark'], **tp_dict)
        test_programs[section] = program

        if len(test_programs) == 1:
            # only one program; set default program which helpfully is the most
            # recent value of section from the previous loop.
            user_options['default_program'] = section

    return (user_options, test_programs)

def parse_jobconfig(config_file, user_options, test_programs, settings=None):
    '''Parse the test configurations from the jobconfig file.

config_file: location of the jobconfig file, either relative or absolute.'''

    if not os.path.exists(config_file):
        raise exceptions.TestCodeError(
                'Job configuration file %s does not exist.' % (config_file)
                                      )

    # paths to the test directories can be specified relative to the config
    # file.
    config_directory = os.path.dirname(os.path.abspath(config_file))

    jobconfig = compat.configparser.RawConfigParser()
    jobconfig.optionxform = str # Case sensitive file.
    jobconfig.read(config_file)

    # Alter config file with additional settings provided.
    if settings:
        for (section_key, section) in settings.items():
            for (option_key, value) in section.items():
                jobconfig.set(section_key, option_key, value)

    # Parse job categories.
    # Just store as list of test names for now.
    if jobconfig.has_section('categories'):
        test_categories = dict(jobconfig.items('categories'))
        for (key, val) in test_categories.items():
            test_categories[key] = val.split()
        jobconfig.remove_section('categories')
    else:
        test_categories = {}

    # Parse individual sections for tests.
    # Note that sections/paths may contain globs and hence correspond to
    # multiple tests.
    # First, find out the tests each section corresponds to.
    test_sections = []
    for section in jobconfig.sections():
        # Expand any globs in the path/section name and create individual Test
        # objects for each one.
        if jobconfig.has_option(section, 'path'):
            path = os.path.join(config_directory,
                                jobconfig.get(section, 'path'))
            jobconfig.remove_option(section, 'path')
            globbed_tests = [(section, os.path.abspath(test_path))
                                            for test_path in glob.glob(path)]
        else:
            path = os.path.join(config_directory, section)
            globbed_tests = [(test_path, os.path.abspath(test_path))
                                            for test_path in glob.glob(path)]
        test_sections.append((section, globbed_tests))
    test_sections.sort(key=lambda sec_info: len(sec_info[1]), reverse=True)
    test_info = {}
    for (section, globbed_tests) in test_sections:
        test_dict = {}
        # test program
        if jobconfig.has_option(section, 'program'):
            test_program = test_programs[jobconfig.get(section, 'program')]
        else:
            test_program = test_programs[user_options['default_program']]
        # tolerances
        if jobconfig.has_option(section, 'tolerance'):
            test_dict['tolerances'] = {}
            for item in (
                    eval_nested_tuple(jobconfig.get(section,'tolerance'))
                        ):
                (name, tol) = parse_tolerance_tuple(item)
                test_dict['tolerances'][name] = tol
            jobconfig.remove_option(section, 'tolerance')
            if None in test_dict['tolerances']:
                test_dict['default_tolerance'] = test_dict['tolerances'][None]
        # inputs and arguments
        if jobconfig.has_option(section, 'inputs_args'):
            # format: (input, arg), (input, arg)'
            test_dict['inputs_args'] = (
                    eval_nested_tuple(jobconfig.get(section, 'inputs_args')))
            jobconfig.remove_option(section, 'inputs_args')
        if jobconfig.has_option(section, 'run_concurrent'):
            test_dict['run_concurrent'] = \
                    jobconfig.getboolean(section, 'run_concurrent')
            jobconfig.remove_option(section, 'run_concurrent')
        # Other options.
        for option in jobconfig.options(section):
            test_dict[option] = jobconfig.get(section, option)
        for key in ('nprocs', 'max_nprocs', 'min_nprocs'):
            if key in test_dict:
                test_dict[key] = int(test_dict[key])
        if 'submit_template' in test_dict:
            test_dict['submit_template'] = os.path.join(config_directory,
                                                   test_dict['submit_template'])
        for (name, path) in globbed_tests:
            # Need to take care with tolerances: want to *update* existing
            # tolerance dictionary rather than overwrite it.
            # This means we can't just use test_dict to update the relevant
            # dictionary in test_info.
            tol = None
            if (name, path) in test_info:
                # Just update existing info.
                test = test_info[(name, path)]
                if  'tolerances' in test_dict:
                    test[1]['tolerances'].update(test_dict['tolerances'])
                    tol = test_dict.pop('tolerances')
                test[0] = test_program
                test[1].update(test_dict)
                if tol:
                    test_dict['tolerances'] = tol
            else:
                # Create new test_info value.
                # Merge with default values.
                # Default test options.
                default_test = test_program.default_test_settings
                test = dict(
                        inputs_args=default_test.inputs_args,
                        output=default_test.output,
                        default_tolerance=default_test.default_tolerance,
                        tolerances = copy.deepcopy(default_test.tolerances),
                        nprocs=default_test.nprocs,
                        min_nprocs=default_test.min_nprocs,
                        max_nprocs=default_test.max_nprocs,
                        run_concurrent=default_test.run_concurrent,
                        submit_template=default_test.submit_template,
                    )
                if  'tolerances' in test_dict:
                    test['tolerances'].update(test_dict['tolerances'])
                    tol = test_dict.pop('tolerances')
                test.update(test_dict)
                # restore tolerances for next test in the glob.
                if tol:
                    test_dict['tolerances'] = tol
                test_info[(name, path)] = [test_program, copy.deepcopy(test)]

    # Now create the tests (after finding out what the input files are).
    tests = []
    for ((name, path), (test_program, test_dict)) in test_info.items():
        old_dir = os.getcwd()
        os.chdir(path)
        # Expand any globs in the input files.
        inputs_args = []
        for input_arg in test_dict['inputs_args']:
            # Be a little forgiving for the input_args config option.
            # If we're given ('input'), then clearly the user meant for the
            # args option to be empty.  However, literal_eval returns
            # a string rather than a tuple in such cases, which causes
            # problems.
            if isinstance(input_arg, str):
                inp = input_arg
                arg = ''
            elif len(input_arg) == 2:
                inp = input_arg[0]
                arg = input_arg[1]
            else:
                inp = input_arg[0]
                arg = ''
            if inp:
                # the test, error and benchmark filenames contain the input
                # filename, so we need to filter them out.
                inp_files = sorted(glob.glob(inp))
                if not inp_files:
                    err = 'Cannot find input file %s in %s.' % (inp, path)
                    warnings.warn(err)
                    continue
                # We use a glob for the input argument to avoid the
                # case where the argument is empty and hence a pattern
                # such as *.inp also matches files like
                # test.out.test_id.inp=x.inp and hence considering
                # previous output files to actually be an input file in
                # their own right.
                test_files = [
                     util.testcode_filename(stem[1], '*', '*', arg)
                     for stem in testcode2._FILESTEM_TUPLE
                             ]
                testcode_files = []
                for tc_file in test_files:
                    testcode_files.extend(glob.glob(tc_file))
                for inp_file in inp_files:
                    if inp_file not in testcode_files:
                        inputs_args.append((inp_file, arg))
            else:
                inputs_args.append((inp, arg))
        test_dict['inputs_args'] = tuple(inputs_args)
        os.chdir(old_dir)
        # Create test.
        if test_dict['run_concurrent']:
            for input_arg in test_dict['inputs_args']:
                test_dict['inputs_args'] = (input_arg,)
                tests.append(testcode2.Test(name, test_program, path,
                                            **test_dict))
        else:
            tests.append(testcode2.Test(name, test_program, path, **test_dict))

    return (tests, test_categories)

def get_unique_test_id(tests, reuse_id=False, date_fmt='%d%m%Y'):
    '''Find a unique test id based upon the date and previously run tests.'''
    todays_id = time.strftime(date_fmt)
    newest_file = None
    test_id = '0'*len(todays_id)
    for test in tests:
        test_globs = glob.glob('%s*' %
                os.path.join(test.path, testcode2.FILESTEM['test'])
                              )
        for test_file in test_globs:
            if (not newest_file or
                    os.stat(test_file)[-2] > os.stat(newest_file)[-2]):
                newest_file = test_file
                # keep track of the latest file with today's test_id (in case
                # the most recent test was run with a user-specified test_id).
                newest_test_id = util.testcode_file_id(
                                 newest_file, testcode2.FILESTEM['test']
                                                 )
                if newest_test_id[:len(todays_id)] == todays_id:
                    test_id = newest_test_id
    if reuse_id:
        # Want test_id to be the most recent set of tests.
        if not newest_file:
            err = 'Cannot find any previous test outputs.'
            raise exceptions.TestCodeError(err)
        test_id = util.testcode_file_id(newest_file, testcode2.FILESTEM['test'])
    elif test_id[:len(todays_id)] == todays_id:
        # Have run at more than one test today already.  Create unique id.
        if len(test_id) == len(todays_id):
            test_id = 1
        else:
            test_id = int(test_id[len(todays_id)+1:]) + 1
        test_id = '%s-%s' % (todays_id, test_id)
    else:
        # First test of the day!
        test_id = todays_id
    return test_id

def select_tests(all_tests, test_categories, selected_categories, prefix=''):
    '''Return the set of tests contained by the selected test categories.'''
    test_categories['_all_'] = [test.path for test in all_tests]
    if ('_default_' in selected_categories
            and '_default_' not in test_categories):
        selected_categories = ['_all_']
    # Recursively expand job categories.
    while compat.compat_any(
                    cat in test_categories for cat in selected_categories
                           ):
        tmp = []
        for cat in selected_categories:
            if cat in test_categories:
                tmp.extend(test_categories[cat])
            else:
                # cat has been fully expanded and now refers to a test
                # contained within the directory named cat.
                tmp.append(cat)
        selected_categories = tmp
    # Select tests to run.
    tests = []
    parent = lambda pdir, cdir: \
            not os.path.relpath(cdir, start=pdir).startswith(os.pardir)
    for cat in selected_categories:
        # test paths are relative to the config directory but absolute paths
        # are stored .
        found = False
        cat_paths = glob.glob(os.path.join(prefix, cat))
        for test in all_tests:
            if cat == test.name:
                found = True
                tests.append(test)
            elif compat.compat_any(os.path.exists(path) and
                    os.path.samefile(path, test.path) for path in cat_paths):
                found = True
                tests.append(test)
            elif compat.compat_any(parent(path, test.path)
                    for path in cat_paths):
                # test contained within a subdirectory of a cat_path.
                found = True
                tests.append(test)
        if not found:
            print('WARNING: %s test/category not found.\n' % cat)
    # Only want to run each test once.
    tests = list(compat.compat_set(tests))
    return tests

def set_program_name(program, relative_path):
    '''Set a full path to the given program.

If the program exists on PATH, then return the full path to that program.
Otherwise, assume program is given relative to relative_path and hence return
the full path.
'''
    program_path = os.path.join(relative_path, program)
    program_path = os.path.expandvars(program_path)
    if not os.path.exists(program_path):
        # Program not supplied as a relative or full path.
        # Does program exist on the user's path?
        which_popen = subprocess.Popen(['which', program],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        which_popen.wait()
        if which_popen.returncode == 0:
            # Program is on user's path.
            # Return full path to program.
            program_path = which_popen.communicate()[0].decode('utf-8').strip()
        else:
            # Cannot find program.
            # This still allows us to manipulate previously run tests, just not
            # run new ones...
            print('WARNING: cannot find program: %s.' % (program))
            # Allow things to proceed with the original path -- the user might
            # know what they're doing and the above tests are not always
            # sufficient (e.g. if using cygwin but using an MPI implementation
            # which requires a Windows-based path).
            program_path = program

    return program_path
