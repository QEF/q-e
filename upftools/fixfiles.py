#!/usr/bin/env python3
"""
simple script to mark as CDATA all text contained in PP_INPUTFILE element. It
avoids any issue related to the presence on XML reserved characters in that
section.  
"""
from subprocess import Popen, PIPE
from xml.etree import cElementTree as eT
import os


def fix_file(f):
    """
    Reads file f: looks for the INPUTFILE element, corrects it and returns 1.  
    If is doesn't find it prints a message and returns 0.
    raise an exception if the PP_INPUTFILE is not correctly terminated
    :f: is the file name
    """ 
    import os
    with open(f, 'r') as fin:
        lin = fin.readlines()
    lin = [l.strip() for l in lin]
    try:
        start = lin.index('<PP_INPUTFILE>')
    except ValueError:
        print ( 'PP_INPUTFILE not present in %s' % f ) 
        return 0
    try:
        end = lin.index('</PP_INPUTFILE>')
    except ValueError:
        print ('PP_INPUTFILE is not correctly closed in %s' % f)
        raise

    lout = lin[:start+1]+['<![CDATA[']+lin[start+2:end]+[']]>']+lin[end:]
    os.rename(f, '%s_orig' % f)
    with open(f, 'w') as fout:
        for l in lout:
            fout.write(l+'\n')
    return 1


with Popen('ls *.UPF', shell=True, stdout=PIPE) as p:
    out = p.communicate()[0]
upf_files = out.decode().split('\n')

upf_v2 = []
upf_v1 = []

for _ in upf_files:
    with Popen(['grep', '<UPF', _], stdout=PIPE) as p:
        out = p.communicate()[0]
    if out.decode() == '':
        upf_v1.append(_)
    if out.decode() != '':
        upf_v2.append(_)

print('In this directory you  have  {} files UPF v1'.format(len(upf_v1)))
print('In this directory you have {} files UPF v2'.format(len(upf_v2)))

not_well_formed = []
for p in upf_v2:
    try:
        _ = eT.parse(p)
    except eT.ParseError:
        not_well_formed.append(p)


if len(not_well_formed) == 0:
    print ('All UPF v2 files in this directory can  be parsed')
else:
    print(' {} files UPF v2 in this directory need to be corrected'.format(len(not_well_formed)))



count = 0
with open('log', 'a') as log:
    for p in not_well_formed:
        fix_file(p)
        try:
            _  = eT.parse(p)
            log.write("%s has been fixed \n" % p)
        except eT.ParseError:
            log.write("could not fix %s \n" % p)
            try:
                os.rename('%s_orig' % p, p)
            except FileNotFoundError:
                pass

            count += 1 

if count != 0:
    print ('Fixing failed for {} files. See the log !!!'.format(count))
