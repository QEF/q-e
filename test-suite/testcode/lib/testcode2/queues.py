'''
testcode.queues
---------------

Access to external queueing systems.

:copyright: (c) 2012 James Spencer.
:license: modified BSD; see LICENSE for more details.
'''

import os.path
import subprocess
import sys
import time

import testcode2.exceptions as exceptions

class ClusterQueueJob:
    '''Interface to external queueing system.

:param string submit_file: filename of submit script to be submitted to the
    queueing system.
:param string system: name of queueing system.  Currently only an interface to
    PBS is implemented.
'''
    def __init__(self, submit_file, system='PBS'):
        self.job_id = None
        self.submit_file = submit_file
        self.system = system
        if self.system == 'PBS':
            self.submit_cmd = 'qsub'
            self.queue_cmd = 'qstat'
            self.job_id_column = 0
            self.status_column = 4
            self.finished_status = 'C'
        else:
            err = 'Queueing system not implemented: %s' % self.system
            raise exceptions.RunError(err)
    def create_submit_file(self, pattern, string, template):
        '''Create a submit file.
        
Replace pattern in the template file with string and place the result in
self.submit_file.

:param string pattern: string in template to be replaced.
:param string string: string to replace pattern in template.
:param string template: filename of file containing the template submit script.
'''
        # get template
        if not os.path.exists(template):
            err = 'Submit file template does not exist: %s.' % (template,)
            raise exceptions.RunError(err)
        ftemplate = open(template)
        submit = ftemplate.read()
        ftemplate.close()
        # replace marker with our commands
        submit = submit.replace(pattern, string)
        # write to submit script
        fsubmit = open(self.submit_file, 'w')
        fsubmit.write(submit)
        fsubmit.close()
    def start_job(self):
        '''Submit job to cluster queue.'''
        submit_cmd = [self.submit_cmd, self.submit_file]
        try:
            submit_popen = subprocess.Popen(submit_cmd, stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT)
            submit_popen.wait()
            self.job_id = submit_popen.communicate()[0].strip().decode('utf-8')
        except OSError:
            # 'odd' syntax so exceptions work with python 2.5 and python 2.6/3.
            err = 'Error submitting job: %s' % (sys.exc_info()[1],)
            raise exceptions.RunError(err)
    def wait(self):
        '''Returns when job has finished running on the cluster.'''
        running = True
        # Don't ask the queueing system for the job itself but rather parse the
        # output from all current jobs and look  gor the job in question. 
        # This works around the problem where the job_id is not a sufficient
        # handle to query the system directly (e.g. on the CMTH cluster).
        qstat_cmd = [self.queue_cmd]
        while running:
            time.sleep(15)
            qstat_popen = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
            qstat_popen.wait()
            if qstat_popen.returncode != 0:
                err = ('Error inspecting queue system: %s' %
                                                      qstat_popen.communicate())
                raise exceptions.RunError(err)
            qstat_out = qstat_popen.communicate()[0]
            # Assume job has finished unless it appears in the qstat output.
            running = False
            for line in qstat_out.splitlines():
                words = line.split()
                if words[self.job_id_column] == self.job_id:
                    running = words[self.status_column] != self.finished_status
                    break
