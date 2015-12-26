'''
testcode2.compatibility
-----------------------

Functions for compatibility with different versions of python.

testcode2 is developed using python 3.2; these statements exist to enable
testcode to function transparently (i.e. without using 2to3) on python 2.4
onwards.

Rather than using conditional statements in the main source code, instead place
the required statements in this module and then use (e.g.)

import testcode2.compatibility as compat

var = compat.compat_set([1,2,3,1])

in the main source code.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import sys

### python 2.4 ###

# Import from the sets module if sets are not part of the language.
try:
    compat_set = set
except NameError:
    from sets import Set as compat_set

# Any and all don't exist in python <2.5. Define our own in pure python.
try:
    compat_all = all
except NameError:
    def compat_all(iterable):
        '''all(iterable) -> bool

Return True if bool(x) is True for all values x in the iterable.
'''
        for val in iterable:
            if not val:
                return False
        return True
try:
    compat_any = any
except NameError:
    def compat_any(iterable):
        '''any(iterable) -> bool

Return True if bool(x) is True for any x in the iterable.
'''
        for val in iterable:
            if val:
                return True

try:
    import functools
except ImportError:
    import testcode2._functools_dummy as functools

### python 2.4, python 2.5 ###

# math.isnan was introduced in python 2.6, so need a workaround for 2.4 and 2.5.
try:
    from math import isnan
except ImportError:
    def isnan(val):
        '''Return True if x is a NaN (not a number), and False otherwise.

:param float val: number.

Replacement for math.isnan for python <2.6.
This is not guaranteed to be portable, but does work under Linux.
'''
        return type(val) is float and val != val

try:
    # python >=2.6
    from ast import literal_eval
except ImportError:
    # python 2.4, 2.5
    from compiler import parse
    from compiler import ast
    def literal_eval(node_or_string):
        """Safely evaluate a node/string containing a Python expression.

Thestring or node provided may only consist of the following Python literal
structures: strings, numbers, tuples, lists, dicts,  booleans, and None.

Essentially a backport of the literal_eval function in python 2.6 onwards.
From: http://mail.python.org/pipermail/python-list/2009-September/1219992.html
        """
        _safe_names = {'None': None, 'True': True, 'False': False}
        if isinstance(node_or_string, basestring):
            node_or_string = parse(node_or_string, mode='eval')
        if isinstance(node_or_string, ast.Expression):
            node_or_string = node_or_string.node
        def _convert(node):
            '''Convert node/string to expression.'''
            if isinstance(node, ast.Const) and isinstance(node.value,
                    (basestring, int, float, long, complex)):
                return node.value
            elif isinstance(node, ast.Tuple):
                return tuple(_convert(element) for element in node.nodes)
            elif isinstance(node, ast.List):
                return list(_convert(element) for element in node.nodes)
            elif isinstance(node, ast.Dict):
                return dict((_convert(k), _convert(v)) for k, v
                            in node.items)
            elif isinstance(node, ast.Name):
                if node.name in _safe_names:
                    return _safe_names[node.name]
            elif isinstance(node, ast.UnarySub):
                return -_convert(node.expr)
            raise ValueError('malformed string')
        return _convert(node_or_string)

# os.path.relpath was introduced in python 2.6.
try:
    from os.path import relpath
except ImportError:
    import os.path
    def relpath(path, start=os.path.curdir):
        """Return a relative version of a path"""

        if not path:
            raise ValueError("no path specified")

        filter_null = lambda lll: [x for x in lll if x]

        start_list = filter_null(os.path.abspath(start).split(os.path.sep))
        path_list = filter_null(os.path.abspath(path).split(os.path.sep))

        common = len(os.path.commonprefix([start_list, path_list]))

        rel_list = [os.pardir] * (len(start_list)-common) + path_list[common:]
        if not rel_list:
            return os.path.curdir
        return os.path.join(*rel_list)

### python 2 ###

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

try:
    compat_input = raw_input
except NameError:
    compat_input = input

try:
    maxint = sys.maxint
except AttributeError:
    maxint = sys.maxsize
