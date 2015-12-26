.. _config:

Configuration files
===================

For convenience, tests can be specified via configuration files rather than
using the testcode API directly.  These configuration files are required for
work with the command-line interface.

The two configuration files are, by default, :ref:`jobconfig` and
:ref:`userconfig` in the working directory.  Different names and/or paths can
be specified if required.

Both configuration files take options in the ini format (as understood by
Python's `configparser <http://docs.python.org/library/configparser.html>`_ module).  For example::

    [section_1]
    a = 2
    b = test_option

    [section_2]
    v = 4.5
    hello = world

defines an ini file with two sections (named 'section_1' and 'section_2'), each
with two variables set.

.. note::

    Any paths can either be absolute or relative to the directory containing
    the configuration file.  The full path need not be given for any program
    which exists on the user's PATH.  Environment variables in **program** names will be expanded.
