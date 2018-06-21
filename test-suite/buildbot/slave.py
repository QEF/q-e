#
# 2016-2018 : Samuel Ponce' and Martin Schlipf
# 
# Setup used by the different Buildbot slaves. 
# 

from buildbot.plugins import steps
from buildbot.steps.shell import ShellCommand
from buildbot.locks import SlaveLock

class Steps:

  def __init__(self,Environ):
    # Max number of running builds
    build_lock = SlaveLock('build',
         maxCount = 1,
         maxCountForSlave = {
             'farmer-slave1': 1,
    })
    
    # All repo
    all_repos = {
        'quantum_espresso': {
            'repository': 'https://gitlab.com/QEF/q-e.git',
            'branch': 'develop',
        },
        'sternheimer_gw': {
            'repository': 'https://github.com/mmdg-oxford/SternheimerGW.git',
            'branch': 'develop',
        },
        'wannier90': {
            'repository': 'https://github.com/wannier-developers/wannier90.git',
            'branch': 'develop',
        },
    }

############################################################################
# QE code
############################################################################
  
    self.checkout_qe = [steps.Git(
                 name="checkout_qe",
                 method="copy",
                 repourl=all_repos["quantum_espresso"]["repository"],
                 branch=all_repos["quantum_espresso"]["branch"],
                 haltOnFailure = True,
                 alwaysUseLatest = True,
             )]
  
    self.configure_qe = [ShellCommand(
                   name="configure_qe",
                   command=["./configure"], 
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe"]
               )]
    
    self.dep_qe     = [ShellCommand(
                   name="dep_qe",
                   #command=["make","depend"], 
                   # DBSP: Temporary until QE 6.2 is release 
                   command=["ls"], 
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["dep_qe"]
               )]
    
    self.make_pw    = [ShellCommand(
                   name="make_pw",
                   command=["make","pwall","cp","ld1","upf"], 
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True, descriptionDone=["make_pw"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.make_ph    = [ShellCommand(
                   name="make_ph",
                   command=["make","ph"], 
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True, descriptionDone=["make_ph"],
                   locks=[build_lock.access('counting')]
                )]
    
    
    self.make_epw0   = [ShellCommand(
                   name="make_epw0",
                   command=["make"], 
                   env=Environ,
                   workdir="build/EPW/src/",
                   haltOnFailure=True, descriptionDone=["make_epw"],
                   locks=[build_lock.access('counting')]
                )]
    
    
    self.make_epw   = [ShellCommand(
                   name="make_epw",
                   command=["make","epw"], 
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True, descriptionDone=["make_epw"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.make_lr    = [ShellCommand(
                   name="make_lr",
                   command=["make","-j","8","lrmods"],
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True,
                   descriptionDone=["make_lr"],
                   locks=[build_lock.access('counting')],
                )]
    
    self.test_clean = [ShellCommand(
                  name="test_clean",
                  command=["make", "clean"],
                  env=Environ,
                  workdir="build/test-suite",
                  descriptionDone = ["test_clean"],
                  locks=[build_lock.access('counting')],
                )]

    self.clean  = [ShellCommand(
                   command=["make", "veryclean"],
                   alwaysRun=True,
                   flunkOnFailure = False,
                   workdir="build"
               )]
    
    self.test0      = [ShellCommand(
                   name="test_prolog",
                   command=["make","prolog"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["make prolog"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.test_para_PW = [ShellCommand(
                   name="PW_para",
                   command=["make","run-tests-pw-parallel"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PW para tests"],
                   locks=[build_lock.access('counting')]
                )]
    self.test_serial_PW = [ShellCommand(
                   name="PW_serial",
                   command=["make","run-tests-pw-serial"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PW serial tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_para_CP = [ShellCommand(
                   name="CP_para",
                   command=["make","run-tests-cp-parallel"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["CP para tests"],
                   locks=[build_lock.access('counting')]
                )]
    self.test_serial_CP = [ShellCommand(
                   name="CP_serial",
                   command=["make","run-tests-cp-serial"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["CP serial tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_para_PH = [ShellCommand(
                   name="PH_para",
                   command=["make","run-tests-ph-parallel"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PH para tests"],
                   locks=[build_lock.access('counting')]
                )]
    self.test_serial_PH = [ShellCommand(
                   name="PH_serial",
                   command=["make","run-tests-ph-serial"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PH serial tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_para_EPW  = [ShellCommand(
                   name="EPW_para",
                   command=["make","run-tests-epw-parallel"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["EPW para tests"],
                   locks=[build_lock.access('counting')]
                )]
    self.test_serial_EPW  = [ShellCommand(
                   name="EPW_serial",
                   command=["make","run-tests-epw-serial"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["EPW serial tests"],
                   locks=[build_lock.access('counting')]
                )]


############################################################################
# SGW code
############################################################################
    self.configure_qe2 = [ShellCommand(
                   name="configure_qe",
                   command=["./configure"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe"]
               )]

    self.make_pw2   = [ShellCommand(
                   name="make_pw",
                   command=["make","pw","lrmods"], 
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True, descriptionDone=["make_pw"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.checkout_sgw = [steps.Git(
                   name="checkout_sgw",
                   repourl=all_repos["sternheimer_gw"]["repository"],
                   branch=all_repos["sternheimer_gw"]["branch"],
                   workdir="build/SGW",
                   haltOnFailure = True,
                   alwaysUseLatest = True,
                )]

    self.make_sgw   = [ShellCommand(
                  name="make_sgw",
                  command=["make"],
                  env=Environ,
                  workdir="build/SGW",
                  haltOnFailure = True,
                  descriptionDone = ["make_sgw"],
                  locks=[build_lock.access('counting')],
                )]

    self.test_sgw   = [ShellCommand(
                  name="test_sgw",
                  command=["make", "run-tests"],
                  env=Environ,
                  workdir="build/SGW/test-suite",
                  haltOnFailure = True,
                  descriptionDone = ["test_sgw"],
                  locks=[build_lock.access('counting')],
                )]

    self.test_clean_sgw = [ShellCommand(
                  name="test_clean",
                  command=["make", "clean"],
                  env=Environ,
                  workdir="build/SGW/test-suite",
                  descriptionDone = ["test_clean"],
                  locks=[build_lock.access('counting')],
                )]

    
############################################################################
# Wannier code
############################################################################
    
    self.checkout_wannier = [steps.Git(
                   name="checkout_wannier",
                   method="copy",
                   workdir="build/WAN",
                   repourl=all_repos["wannier90"]["repository"],
                   branch=all_repos["wannier90"]["branch"],
                   haltOnFailure = True,
                   alwaysUseLatest = True,
                )]
    
    self.cpconfig    = [ShellCommand(
                   name="cp_config",
                   command=["cp","test-suite/config/EPW_testfarm/farmer_gcc485.inc","make.inc"], 
                   env=Environ,
                   workdir="build/WAN",
                   haltOnFailure=True, descriptionDone=["cp_config"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.clean_wannier   = [ShellCommand(
                  name="clean_wannier",
                  command=["python","clean_tests"],
                  env=Environ,
                  workdir="build/WAN/test-suite",
                  haltOnFailure = True, 
                  descriptionDone = ["clean_wannier"],
                  locks=[build_lock.access('counting')],
                )]
    
    self.make_wannier   = [ShellCommand(
                  name="make_wannier",
                  command=["make"],
                  env=Environ,
                  workdir="build/WAN",
                  haltOnFailure = True, 
                  descriptionDone = ["make_wannier"],
                  locks=[build_lock.access('counting')],
                )]
    
    self.make_wannier2   = [ShellCommand(
                  name="make_wannier2",
                  command=["make","default","w90chk2chk"],
                  env=Environ,
                  workdir="build/WAN",
                  haltOnFailure = True, 
                  descriptionDone = ["make_wannier2"],
                  locks=[build_lock.access('counting')],
                )]
    
    self.test_wannier_serial   = [ShellCommand(
                  name="test_wannier_seq",
                  command=["./run_tests","--category=default"],
                  env=Environ,
                  workdir="build/WAN/test-suite",
                  haltOnFailure = True, 
                  descriptionDone = ["test_wannier_seq"],
                  locks=[build_lock.access('counting')],
                )]
    
    self.test_wannier_para   = [ShellCommand(
                  name="test_wannier_para",
                  command=["./run_tests","--category=default", "--numprocs=4"],
                  env=Environ,
                  workdir="build/WAN/test-suite",
                  haltOnFailure = True, 
                  descriptionDone = ["test_wannier_para"],
                  locks=[build_lock.access('counting')],
                )]
  

