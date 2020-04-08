
import os

from buildslave.bot import BuildSlave
from twisted.application import service

basedir = '.'
rotateLength = 512000
maxRotatedFiles = 3

# if this is a relocatable tac file, get the directory containing the TAC
if basedir == '.':
    import os.path
    basedir = os.path.abspath(os.path.dirname(__file__))

# note: this line is matched against to check that this is a buildslave
# directory; do not edit it.
application = service.Application('buildslave')

try:
  from twisted.python.logfile import LogFile
  from twisted.python.log import ILogObserver, FileLogObserver
  logfile = LogFile.fromFullPath(os.path.join(basedir, "twistd.log"), rotateLength=rotateLength,
                                 maxRotatedFiles=maxRotatedFiles)
  application.setComponent(ILogObserver, FileLogObserver(logfile).emit)
except ImportError:
  # probably not yet twisted 8.2.0 and beyond, can't set log yet
  pass

buildmaster_host = 'localhost'
port = 9989
slavename = 'farmer-slave2'
passwd = 'XXXXXXXXXXXXXXXXXXXXX'
keepalive = 120
usepty = 0
umask = 027
maxdelay = 30
allow_shutdown = None

s = BuildSlave(buildmaster_host, port, slavename, passwd, basedir,
               keepalive, usepty, umask=umask, maxdelay=maxdelay,
               allow_shutdown=allow_shutdown)
s.setServiceParent(application)

