#!/usr/bin/env python3
#
# Post-processing script from of PH data in format used by EPW
# 14/07/2015 - Creation of the script - Samuel Ponce
# 14/03/2018 - Automatically reads the number of q-points - Michael Waters
# 14/03/2018 - Detect if SOC is included in the calculation - Samuel Ponce
# 05/06/2019 - Removed SOC for xml detection instead - Felix Goudreault
#
from __future__ import print_function
try:
    from builtins import input
except ImportError:
    print('Install future. e.g. "pip install --user future"')
# import numpy as np

import os
import re
from xml.dom import minidom


# Return the number of q-points in the IBZ
def get_nqpt(prefix):
    fname = '_ph0/' + prefix + '.phsave/control_ph.xml'

    fid = open(fname, 'r')
    lines = fid.readlines()
    # these files are relatively small so reading the whole thing shouldn't
    # be an issue
    fid.close()

    line_number_of_nqpt = 0
    while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]:
        # increment to line of interest
        line_number_of_nqpt += 1
    line_number_of_nqpt += 1  # its on the next line after that text

    nqpt = int(lines[line_number_of_nqpt])

    return nqpt


# Check if the calculation include SOC
def hasSOC(prefix):
    fname = prefix+'.save/data-file-schema.xml'

    xmldoc = minidom.parse(fname)
    item = xmldoc.getElementsByTagName('spinorbit')[0]
    lSOC = item.childNodes[0].data

    return lSOC


# Check if the calculation includes PAW
def hasPAW(prefix):
    fname = prefix+'.save/data-file-schema.xml'

    xmldoc = minidom.parse(fname)
    item = xmldoc.getElementsByTagName('paw')[0]
    lPAW = (item.childNodes[0].data == 'true')

    return lPAW


# Check if the calculation used .fc or .fc.xml files
def hasfc(prefix):
    fname = str(prefix)+'.fc.xml'
    if (os.path.isfile(fname)):
        lfc = True
    else:
        fname_no_xml = re.sub('\.xml$', '', fname)
        if (os.path.isfile(fname_no_xml)):
            lfc = True
        else:
            lfc = False

    return lfc


# check if calculation used xml files (irrelevant of presence of SOC)
def hasXML(prefix):
    # check for a file named prefix.dyn1.xml
    # if it exists => return True else return False
    fname = os.path.join(prefix + ".dyn1.xml")
    if os.path.isfile(fname):
        return True
    # check if the other without .xml extension exists
    # if not raise an error
    fname_no_xml = re.sub('\.xml$', '', fname)

    class FileNotFoundError(Exception):
        pass
    if not os.path.isfile(fname_no_xml):
        raise FileNotFoundError(
                "No dyn0 file found cannot tell if xml format was used.")
    return False


# Check if the calculation was done in sequential
def isSEQ(prefix):
    fname = '_ph0/'+str(prefix)+'.dvscf'
    if (os.path.isfile(fname)):
        lseq = True
    else:
        lseq = False

    return lseq


# Enter the number of irr. q-points
user_input = input(
        'Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = str(user_input)

# # Test if SOC
# SOC = hasSOC(prefix)
# Test if '.xml' files are used
XML = hasXML(prefix)

# Test if PAW
PAW = hasPAW(prefix)

# Test if fc
fc = hasfc(prefix)

# Test if seq. or parallel run
SEQ = isSEQ(prefix)

if True:  # this gets the nqpt from the outputfiles
    nqpt = get_nqpt(prefix)

else:
    # Enter the number of irr. q-points
    user_input = input(
            'Enter the number of irreducible q-points\n')
    nqpt = user_input
    try:
        nqpt = int(user_input)
    except ValueError:
        raise Exception('The value you enter is not an integer!')

os.system('mkdir save 2>/dev/null')

for iqpt in range(1, nqpt+1):
    label = str(iqpt)

    # Case calculation in seq.
    if SEQ:
        # Case with XML files
        if XML:
            os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
            os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix
                      + '.dyn_q'+label+'.xml')
            if (iqpt == 1):
                os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'
                          + label)
                os.system('cp -r _ph0/'+prefix+'.phsave save/')
                if fc:
                    os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.dvscf_paw* save/'+prefix +
                              '.dvscf_paw_q'+label)
            else:
                os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                          '.dvscf* save/'+prefix+'.dvscf_q'+label)
                os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf_paw* save/'+prefix+'.dvscf_paw_q'+label)
        # Case without XML files
        else:
            os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q' +
                      label)
            if (iqpt == 1):
                os.system('cp _ph0/'+prefix+'.dvscf save/'+prefix+'.dvscf_q' +
                          label)
                os.system('cp -r _ph0/'+prefix+'.phsave save/')
                if fc:
                    os.system('cp '+prefix+'.fc save/ifc.q2r')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.dvscf_paw save/'+prefix +
                              '.dvscf_paw_q'+label)
            else:
                os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                          '.dvscf save/'+prefix+'.dvscf_q'+label)
                os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf_paw save/'+prefix+'.dvscf_paw_q'+label)
    else:
        # Case with XML format
        if XML:
            os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
            os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix +
                      '.dyn_q'+label+'.xml')
            if (iqpt == 1):
                os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q' +
                          label)
                os.system('cp -r _ph0/'+prefix+'.phsave save/')
                if fc:
                    os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.dvscf_paw1 save/'+prefix +
                              '.dvscf_paw_q'+label)
            else:
                os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                          '.dvscf1 save/'+prefix+'.dvscf_q'+label)
                os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)
        # Case without XML format
        else:
            os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q' +
                      label)
            if (iqpt == 1):
                os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q' +
                          label)
                os.system('cp -r _ph0/'+prefix+'.phsave save/')
                if fc:
                    os.system('cp '+prefix+'.fc save/ifc.q2r')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.dvscf_paw1 save/'+prefix +
                              '.dvscf_paw_q'+label)
            else:
                os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                          '.dvscf1 save/'+prefix+'.dvscf_q'+label)
                os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                if PAW:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)
