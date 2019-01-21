#!/usr/bin/python
#
# Post-processing script from of PH data in format used by EPW
# 14/07/2015 - Creation of the script - Samuel Ponce
# 14/03/2018 - Automatically reads the number of q-points - Michael Waters
# 14/03/2018 - Detect if SOC is included in the calculation - Samuel Ponce 
# 13/11/2018 - Write dyn files in xml format for SOC case - Shunhong Zhang (USTC)
# 
import numpy as np
import os
from xml.dom import minidom

# Convert the dyn files to the xml form, for SOC case - Shunhong Zhang (USTC)
def dyn2xml(prefix):
    ndyn=int(os.popen('head -2 {0}.dyn0|tail -1'.format(prefix)).read())
    for idyn in range(1,ndyn+1):
        print '{0}.dyn{1} to {0}.dyn_q{1}.xml'.format(prefix,idyn)
        dynmat=dyn(prefix,idyn)
        dynmat._write_xml()
def get_geom_info():
    if os.path.isfile('ph.out')==False:
       print 'cannot extract geometry info from ph.out'
       return 1
    else:
       volm=float(os.popen('grep -a volume ph.out 2>/dev/null|tail -1').readline().split()[-2])
       get_at=os.popen('grep -a -A 3 "crystal axes" ph.out 2>/dev/null|tail -3').readlines()
       at=np.array([[float(item) for item in line.split()[3:6]] for line in get_at])
       get_bg=os.popen('grep -a -A 3 "reciprocal axes" ph.out 2>/dev/null|tail -3').readlines()
       bg=np.array([[float(item) for item in line.split()[3:6]] for line in get_bg])
       return volm,at,bg

class dyn(object):
    def __init__(self,prefix,idyn):
        self._prefix=prefix
        self._idyn=idyn
        fil='{0}.dyn{1}'.format(prefix,idyn)
        f=open(fil)
        self._comment=f.readline()
        f.readline()
        line=f.readline().split()
        self._ntype=int(line[0])
        self._natom=int(line[1])
        self._ibrav=int(line[2])
        self._nspin=1
        self._cell_dim=np.array([float(ii) for ii in line[3:]])
        self._volm=0
        self._at=np.zeros((3,3),float)
        self._bg=np.zeros((3,3),float)
        try: self._volm,self._at,self._bg = get_geom_info()
        except: print 'warning: lattice info not found'
        self._species=[];
        self._mass=[]
        for i in range(self._ntype):
            line=f.readline().split()
            self._species.append(line[1].strip("'"))
            self._mass.append(float(line[-1])/911.4442)  # normalize to atomic mass
        self._atom_type=np.zeros(self._natom,int)
        self._pos=np.zeros((self._natom,3),float)
        for i in range(self._natom):
            line=f.readline().split()
            self._atom_type[i]=int(line[1])
            for j in range(3): self._pos[i,j]=float(line[j+2])
        self._nqpt=int(os.popen('grep -c "Dynamical  Matrix" {0}'.format(fil)).read().split()[0])
        self._qpt=[]
        self._dynmat=np.zeros((self._nqpt,self._natom,self._natom,3,3,2),float)
        f.readline()
        for iqpt in range(self._nqpt):
            f.readline();
            f.readline()
            line=f.readline().split()
            self._qpt.append(np.array([float(item) for item in line[3:6]]))
            f.readline()
            for i in range(self._natom):
                for j in range(self._natom):
                    f.readline()
                    data=np.fromfile(f,sep=' ',count=18,dtype=float).reshape(3,3,2)
                    self._dynmat[iqpt,i,j]=data
        self._qpt=np.array(self._qpt)
        for i in range(5): f.readline()
        self._freq=np.zeros((self._natom*3,2),float)
        self._disp=np.zeros((self._natom*3,self._natom,3,2),float)
        for i in range(self._natom*3):
            line=f.readline().split()
            self._freq[i,0]=float(line[4])
            self._freq[i,1]=float(line[7])
            for j in range(self._natom):
                line=f.readline().split()[1:-1]
                data=np.array([float(item) for item in line]).reshape(3,2)
                self._disp[i,j]=data

    def _write_xml(self):
        doc=minidom.Document()
        root = doc.createElement('Root')
        doc.appendChild(root)
        geom_info=doc.createElement('GEOMETRY_INFO')
        tags=('NUMBER_OF_TYPES','NUMBER_OF_ATOMS','BRAVAIS_LATTICE_INDEX','SPIN_COMPONENTS')
        numbers=(self._ntype,self._natom,self._ibrav,self._nspin)
        for i,(tag,num) in enumerate(zip(tags,numbers)):
            inode=doc.createElement(tag)
            inode.setAttribute('type','integer')
            inode.setAttribute('size','1')
            inode.text=num
            inode.appendChild(doc.createTextNode(str(num)))
            geom_info.appendChild(inode)
        cell_dim=doc.createElement('CELL_DIMENSIONS')
        cell_dim.setAttribute('type','real')
        cell_dim.setAttribute('size','6')
        for i in range(6):
            cell_dim.appendChild(doc.createTextNode('{0:16.10f}'.format(self._cell_dim[i])))
        geom_info.appendChild(cell_dim)
        tags=['AT','BG']
        for tag,lat in zip(tags,(self._at,self._bg)):
            inode=doc.createElement(tag)
            inode.setAttribute('type','real')
            inode.setAttribute('size','9')
            inode.setAttribute('columns','3')
            for i in range(3):
                text=' '.join(['{0:16.10f}'.format(item) for item in lat[i]])
                inode.appendChild(doc.createTextNode(text))
            geom_info.appendChild(inode)
        volm=doc.createElement('UNIT_CELL_VOLUME_AU')
        volm.setAttribute('type','real')
        volm.setAttribute('size','1')
        volm.appendChild(doc.createTextNode('{0:16.10f}'.format(self._volm)))
        geom_info.appendChild(volm)
        for itype in range(self._ntype):
            nt=doc.createElement('TYPE_NAME.{0}'.format(itype+1))
            nt.setAttribute('type','character')
            nt.setAttribute('size','1')
            nt.setAttribute('len','3')
            nt.appendChild(doc.createTextNode('{0}'.format(self._species[itype])))
            na=doc.createElement('MASS.{0}'.format(itype+1))
            na.setAttribute('type','real')
            na.setAttribute('size','1')
            na.appendChild(doc.createTextNode('{0:16.10f}'.format(self._mass[itype])))
            geom_info.appendChild(nt)
            geom_info.appendChild(na)
        for iat in range(self._natom):
            at=doc.createElement('ATOM.{0}'.format(iat+1))
            at.setAttribute('SPECIES','{0}'.format(self._species[self._atom_type[iat]-1]))
            at.setAttribute('INDEX',str(iat+1))
            pos=' '.join(['{0:16.10f}'.format(item) for item in self._pos[iat]])
            at.setAttribute('TAU',pos)
            geom_info.appendChild(at)
        nqpt=doc.createElement('NUMBER_OF_Q')
        nqpt.setAttribute('type','integer')
        nqpt.setAttribute('size','1')
        nqpt.appendChild(doc.createTextNode(str(self._nqpt)))
        geom_info.appendChild(nqpt)
        root.appendChild(geom_info)
        for iqpt in range(self._nqpt):
            dynmat=doc.createElement('DYNAMICAL_MAT_.{0}'.format(iqpt+1))
            qpt=doc.createElement('Q_POINT')
            qpt.setAttribute('type','real')
            qpt.setAttribute('size','3')
            qpt.setAttribute('columns','3')
            tnode=doc.createTextNode(' '.join(['{0:16.10f}'.format(item) for item in self._qpt[iqpt]]))
            qpt.appendChild(tnode)
            dynmat.appendChild(qpt)
            for iat in range(self._natom):
                for jat in range(self._natom):
                    ph=doc.createElement('PHI.{0}.{1}'.format(iat+1,jat+1))
                    ph.setAttribute('type','complex')
                    ph.setAttribute('size','9')
                    ph.setAttribute('columns','3')
                    for i in range(3):
                        for j in range(3):
                            text='{0:16.10f} {1:16.10f}'.format(self._dynmat[iqpt,iat,jat,i,j,0],self._dynmat[iqpt,iat,jat,i,j,1])
                            ph.appendChild(doc.createTextNode(text))
                    dynmat.appendChild(ph)
            root.appendChild(dynmat)
        mode=doc.createElement('FREQUENCIES_THZ_CMM1')
        for iomega in range(self._natom*3):
            inode=doc.createElement('OMEGA.{0}'.format(iomega+1))
            inode.setAttribute('type','real')
            inode.setAttribute('size','2')
            inode.setAttribute('columns','2')
            inode.appendChild(doc.createTextNode('{0:16.10f} {1:16.10f}'.format(self._freq[iomega,0],self._freq[iomega,1])))
            idisp=doc.createElement('DISPLACEMENT.{0}'.format(iomega+1))
            idisp.setAttribute('tpye','complex')
            idisp.setAttribute('size','3')
            for iat in range(self._natom):
                for j in range(3):
                    tnode=doc.createTextNode('{0:16.10f} {1:16.10f}'.format(self._disp[iomega,iat,j,0],self._disp[iomega,iat,j,1]))
                    idisp.appendChild(tnode)
            mode.appendChild(inode)
            mode.appendChild(idisp)
        root.appendChild(mode)
        fp = open('{0}.dyn_q{1}.xml'.format(self._prefix,self._idyn), 'w')
        doc.writexml(fp, addindent='  ', newl='\n')

# Return the number of q-points in the IBZ
def get_nqpt(prefix):
  fname = '_ph0/' +prefix+'.phsave/control_ph.xml'

  fid = open(fname,'r')
  lines = fid.readlines() # these files are relatively small so reading the whole thing shouldn't be an issue
  fid.close()

  line_number_of_nqpt = 0
  while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]: # increment to line of interest
    line_number_of_nqpt +=1
  line_number_of_nqpt +=1 # its on the next line after that text

  nqpt = int(lines[line_number_of_nqpt])

  return nqpt

# Check if the calculation include SOC
def hasSOC(prefix):
  fname = prefix+'.save/data-file-schema.xml'

  xmldoc = minidom.parse(fname)
  item = xmldoc.getElementsByTagName('spinorbit')[0]
  lSOC = item.childNodes[0].data
  
  return lSOC

# Check if the calculation was done in sequential
def isSEQ(prefix):
  fname = '_ph0/'+str(prefix)+'.dvscf'
  if (os.path.isfile(fname)):
    lseq = True
  else:
    lseq = False
 
  return lseq
    
# Enter the number of irr. q-points
user_input = raw_input('Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = str(user_input)

# Test if SOC
SOC = hasSOC(prefix)

# If SOC detected, but dyn is not in XML and we want to convert it
if SOC=='true':
  user_input = raw_input('Calculation with SOC detected. Do you want to convert dyn in XML format [y/n]?\n')
  if str(user_input) == 'y':
    dyn2xml(prefix)
    os.system('mv {0}.dyn*.xml save'.format(prefix))

# If no SOC detected, do you want to convert into XML format 
if SOC=='false':
  user_input = raw_input('Calculation without SOC detected. Do you want to convert to xml anyway [y/n]?\n')
  if str(user_input) == 'y':
    SOC = 'true'
    dyn2xml(prefix)
    os.system('mv {0}.dyn*.xml save'.format(prefix))

# Test if seq. or parallel run
SEQ = isSEQ(prefix)

if True: # this gets the nqpt from the outputfiles
  nqpt =  get_nqpt(prefix)

else:
  # Enter the number of irr. q-points
  user_input = raw_input('Enter the number of irreducible q-points\n')
  nqpt = user_input
  try:
    nqpt = int(user_input)
  except ValueError:
    raise Exception('The value you enter is not an integer!')

os.system('mkdir save 2>/dev/null')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  # Case calculation in seq.
  if SEQ:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
 
  else:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )

