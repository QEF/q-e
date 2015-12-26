'''
testcode2._dummy_functools
--------------------------

Dummy stub functions of required functools objects used.

This means that we can use python 2.4 and advanced features in later versions
of python.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

def wraps(func1):
    '''Upgrade from python 2.4 to use functools.wraps.'''
    def wrapper(func2):
        '''Upgrade from python 2.4 to use functools.wraps.'''
        def decorated_func(*args, **kwargs):
            '''Upgrade from python 2.4 to use functools.wraps.'''
            return func2(*args, **kwargs)
        return decorated_func
    return wrapper
