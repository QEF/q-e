'''
testcode2.vcs
-------------

Lightweight access to required version control system functions.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import os
import subprocess

class VCSRepository(object):
    '''Handle information about a version control repository.

vcs: version control system used.  Currently git, mercurial and subversion are supported.
repository: (local) directory containing a checked-out version of the repository.
remote_repository: remote location of the repository.
'''
    def __init__(self, vcs, repository, remote_repository=None):
        if vcs in ['svn', 'git', 'hg']:
            self.vcs = vcs
        else:
            self.vcs = None
        self.repository = repository
        if remote_repository:
            self.remote_repository = remote_repository

    def get_code_id(self):
        '''Return the id (i.e. version number or hash) of the VCS repository.'''
        old_dir = os.getcwd()
        os.chdir(self.repository)
        code_id = 'UNKNOWN'
        id_popen = None
        if self.vcs == 'svn':
            id_popen = subprocess.Popen(['svnversion', '.'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif self.vcs == 'git':
            id_popen = subprocess.Popen(['git', 'rev-parse', '--short', 'HEAD'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif self.vcs == 'hg':
            id_popen = subprocess.Popen(['hg', 'id', '-i'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if id_popen:
            id_popen.wait()
            code_id = id_popen.communicate()[0].decode('utf-8').strip()
        os.chdir(old_dir)
        return (code_id)
