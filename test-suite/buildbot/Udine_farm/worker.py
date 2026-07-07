from buildbot.plugins import steps
from buildbot.steps.shell import ShellCommand
from buildbot.locks import WorkerLock
from buildbot.process.properties import Interpolate


class Steps:

  def __init__(self,Environ):
    # Max number of running builds
    build_lock = WorkerLock('build',
         maxCount = 2,
         maxCountForWorker = {
             'farmer-worker1': 2,
    })
    
    # All repo
    all_repos = {
        'quantum_espresso': {
            'repository': 'https://gitlab.com/QEF/q-e.git',
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
                 env=Environ,
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

    self.configure_qe_hdf5 = [ShellCommand(
                   name="configure_qe_hdf5",
                   command=["./configure","--with-hdf5=/home/buildbot/local/hdf5-146-gcc102"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe_hdf5"]
               )]    

    self.configure_qe_libxc = [ShellCommand(
                   name="configure_qe_libxc",
                   command=["./configure", "--with-libxc=yes", "--with-libxc-prefix=/home/buildbot/local/libxc"],
                   env=Environ,
                   workdir="build", 
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True, 
                   descriptionDone=["configure_qe_libxc"] 
               )] 

    self.configure_qe_serial = [ShellCommand(
                   name="configure_qe_serial",
                   command=["./configure","--disable-parallel","--enable-openmp"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe_serial"]
               )]

    self.configure_qe_mp = [ShellCommand(
                   name="configure_qe_mp",
                   command=["./configure","--enable-openmp","--enable-parallel"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe_mp"]
               )]

    self.configure_qe_GPU = [ShellCommand(
                   name="configure_qe_GPU",
                   command=["./configure","--with-cuda=$CUDA_HOME","--with-cuda-runtime=11.4","--with-cuda-cc=70","--disable-parallel"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe_GPU"]
               )]

    self.configure_qe_GPU2 = [ShellCommand(
                   name="configure_qe_GPU2",
                   command=["./configure","--with-cuda=$CUDA_HOME","--with-cuda-runtime=11.7","--with-cuda-cc=70","--enable-parallel","--with-cuda-mpi=yes"],
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["configure_qe_GPU2"]
               )]

    self.debug_flags  = [ShellCommand(
                   name="debug_flags",
                   command=Interpolate('sed -i "s/FFLAGS         = -O3 -g/FFLAGS         = -g -Wall -fbounds-check -frange-check -finit-integer=987654321 -finit-real=nan -finit-logical=true -finit-character=64/g" make.inc'), 
                   env=Environ,
                   workdir="build",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,descriptionDone=["debug_flags"]
               )]

    self.env_omp3    = [ShellCommand(
                   name="env_qe1",
                   command=Interpolate('sed -i "s/OMP_NUM_THREADS=1/OMP_NUM_THREADS=3/g" ENVIRONMENT'),
                   env=Environ,
                   workdir="build/test-suite/",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,
                   descriptionDone=["env_qe1"]
               )]

    self.env_omp2    = [ShellCommand(
                   name="env_qe1",
                   command=Interpolate('sed -i "s/OMP_NUM_THREADS=1/OMP_NUM_THREADS=2/g" ENVIRONMENT'),
                   env=Environ,
                   workdir="build/test-suite/",
                   locks=[build_lock.access('counting')],
                   haltOnFailure = True,
                   descriptionDone=["env_qe1"]
               )]

    self.make_pw    = [ShellCommand(
                   name="make_pw",
                   command=["make","-j","4","pwall","cp","ld1","hp"], 
                   env=Environ,
                   workdir="build",
                   haltOnFailure=True, descriptionDone=["make_pw"],
                   locks=[build_lock.access('counting')]
                )]

    self.make_pw_GPU = [ShellCommand(
                   name="make_pw_GPU",
                   command=["make","-j","4","pw"],
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
    
    self.test_PW = [ShellCommand(
                   name="PW_test",
                   command=["make","run-tests-pw"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PW para tests"],
                   locks=[build_lock.access('counting')]
                )]
    
    self.test_PP = [ShellCommand(
                   name="PP_test",
                   command=["make","run-tests-pp"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PP para tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_CP = [ShellCommand(
                   name="CP_test",
                   command=["make","run-tests-cp"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["CP para tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_PH = [ShellCommand(
                   name="PH_test",
                   command=["make","run-tests-ph"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["PH para tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_para_IMAGE = [ShellCommand(
                   name="IMAGE_para",
                   command=["make","run-tests-image","NPROCS=4"],
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["IMAGE para tests"],
                   locks=[build_lock.access('counting')]
                )]    

    self.test_EPW  = [ShellCommand(
                   name="EPW_test",
                   command=["make","run-tests-epw"], 
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["EPW para tests"],
                   locks=[build_lock.access('counting')]
                )]

    self.test_HP  = [ShellCommand(
                   name="HP_test",
                   command=["make","run-tests-hp"],
                   env=Environ,
                   workdir="build/test-suite",
                   haltOnFailure=False, descriptionDone=["HP para tests"],
                   locks=[build_lock.access('counting')]
                )]

