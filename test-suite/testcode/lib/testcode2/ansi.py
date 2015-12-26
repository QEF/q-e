'''
testcode2.ansi
--------------

(A subset of) ANSI codes and functions to wrap codes around strings.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import sys

ANSI_START = '\033['
ANSI_END = '\033[0m'

ANSI_SGR = dict(
            bold=1,
        )

ANSI_INTENSITY = dict(
            normal=0,
            bright=60,
        )

ANSI_COLOUR = dict(
            black=30,
            red=31,
            green=32,
            yellow=33,
            blue=34,
            magenta=35,
            cyan=36,
            white=37,
        )

def ansi_format(string, colour='black', intensity='normal', style=None,
                override=False):
    '''Return string wrapped in the appropriate ANSI codes.

Note: if not writing to a true terminal and override is false, then the string
is not modified.
'''

    if sys.stdout.isatty() or override:
        code = str(ANSI_COLOUR[colour]+ANSI_INTENSITY[intensity])
        if style:
            code += ';%s' % (ANSI_SGR[style])
        code += 'm'
        return ANSI_START + code + string + ANSI_END + ANSI_END
    else:
        return string
