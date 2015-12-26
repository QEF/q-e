'''
testcode2.exceptions
--------------------

Custom exceptions.  Initialise signal handler for the interrupt signal.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import signal
import sys

def signal_handler(sig, frame):
    '''Capture signal and leave quietly.'''
    print('Signal: %s has been caught.  Bye!' % (sig))
    sys.exit(1)


class RunError(Exception):
    '''Exception used for errors running test jobs.'''
    pass


class AnalysisError(Exception):
    '''Exception used for errors running test jobs.'''
    pass


class TestCodeError(Exception):
    '''Top level exception for testcode errors.'''
    pass

signal.signal(signal.SIGINT, signal_handler) # Listen out for Ctrl-C. 
