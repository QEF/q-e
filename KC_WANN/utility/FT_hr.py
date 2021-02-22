import numpy as np
import sys
from datetime import datetime as dt


#-------------------------------------------------
#--------------- INPUT PARAMETERS ----------------
#-------------------------------------------------

prefix = 'Si'           # prefix as in KC input
num_wann_occ = 4
num_wann_emp = 4
alat = 10.263101844     # lattice parameter in BOHR
kmesh = 4,4,4           # k-mesh from the KC calculation

avec = np.array(([-0.5, 0.0, 0.5],
                 [ 0.0, 0.5, 0.5],
                 [-0.5, 0.5, 0.0]))	# primitive vectors in alat

# K_PATH for band structure
k_path = np.array([[0.  ,0.  ,0.  ],
                   [0.05,0.  ,0.  ],
                   [0.1 ,0.  ,0.  ],
                   [0.15,0.  ,0.  ],
                   [0.2 ,0.  ,0.  ],
                   [0.25,0.  ,0.  ],
                   [0.3 ,0.  ,0.  ],
                   [0.35,0.  ,0.  ],
                   [0.4 ,0.  ,0.  ],
                   [0.45,0.  ,0.  ],
                   [0.5 ,0.  ,0.  ],
                   [0.45,0.  ,0.  ],
                   [0.4 ,0.  ,0.  ],
                   [0.35,0.  ,0.  ],
                   [0.3 ,0.  ,0.  ],
                   [0.25,0.  ,0.  ],
                   [0.2 ,0.  ,0.  ],
                   [0.15,0.  ,0.  ],
                   [0.1 ,0.  ,0.  ],
                   [0.05,0.  ,0.  ],
                   [0.  ,0.  ,0.  ]])

#-------------------------------------------------
#--------------- END INPUT SECTION ---------------
#-------------------------------------------------
#
#
#
#-------------------------------------------------
#------------------- MODULES ---------------------
#-------------------------------------------------

def read_hr(prefix):
    
    ifile = open(prefix+'.kc_hr.dat','r')
    lines = ifile.readlines()
    ifile.close()

    num_wann = int(lines[1].split()[5])
    nrtot = int(lines[1].split()[2])

    if ( len(lines)-2 != nrtot*num_wann*num_wann ): sys.exit('\nInconsistency between nrtot, num_wann and length of hr file\n')

    hr = np.zeros((nrtot,num_wann,num_wann),dtype=complex)
    Rvec = np.zeros((nrtot,3),dtype=int)

    for n in range(2,len(lines)):
        ir = int((n-2)/num_wann**2)
        if ( n%(num_wann**2) ):
            Rvec[ir,:] = np.array((lines[n].split()[:3]),dtype=int)
        iband = int(lines[n].split()[3]) - 1
        jband = int(lines[n].split()[4]) - 1
        hr[ir,iband,jband] = float(lines[n].split()[5]) + 1.j*float(lines[n].split()[6])

    return hr,Rvec,num_wann,nrtot


# This routine reads the centers from .xyz file created by Wannier90
# and fold them within the R=0 primitive cell. The output is in cartesian 
# (divided by alat) units.
def read_centers(prefix,num_wann_occ,num_wann_emp,alat):
    
    centers = np.zeros((num_wann_occ+num_wann_emp,3))
    counter = 0

    for typ in ['occ','emp']:
        if (typ == 'occ'):
            ifile = open(prefix+'_centres.xyz','r')
        else:
            ifile = open(prefix+'_emp_centres.xyz','r')

        lines = ifile.readlines()
        ifile.close()

        for n in range(len(lines)):
            if ( lines[n].split()[0] == 'X' ):
                centers[counter,:] = np.array((lines[n].split()[1:]),dtype=float)
                counter = counter + 1

    centers = centers / (alat*0.52917720859)

    return centers


def corr_phase(mp1,mp2,mp3,avec,kvec,rvect,ws_dist):
    dist_min = np.linalg.norm(ws_dist + rvect)
    Tnn = []                # List of degenerate T-vectors, corresponding to equidistant nn WFs
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                tvect = np.array((i*mp1,j*mp2,k*mp3))       # T is a lattice vector of the SC
                Tvec = crys_to_cart(tvect,avec,+1)
                dist = np.linalg.norm(ws_dist + rvect + Tvec)
                if abs(dist - dist_min) < 1.e-3:
                    Tnn.append(tvect)
                elif (dist < dist_min):   # A new nearest neighbor is found
                    dist_min = dist
                    Tnn = [tvect]
                else:
                    continue
    
    phase = sum(np.exp(1j*2*np.pi*np.dot(Tnn,kvec))) / len(Tnn)
    return phase
    


# Function to transform the vector coordinates from crystal to cartesian (divided by 
# alat), or viceversa, as it is done in QE: typ=+1 for crystal-to-cartesian, typ=-1 
# for cartesian-to-crystal.
# For a real space vector trmat in input must be avec if typ=+1 and bvec if typ=-1.
# For a k-space vector trmat in input must be bvec if typ=+1 and avec if typ=-1.
def crys_to_cart(vec,trmat,typ):    
    if ( typ == +1 ):                   # crystal-to-cartesian conversion
        vec_tr = np.dot(trmat.transpose(),vec)
    elif ( typ == -1 ):                 # cartesian-to-crystal conversion
        vec_tr = np.dot(trmat,vec)
    else:
        sys.exit('\ntyp in crys_to_cart call must be either +1 or -1\n')
    return vec_tr



def MP_mesh(mp1,mp2,mp3):
    k_mesh = []
    for i in range(mp1):
        for j in range(mp2):
            for k in range(mp3):
                k_mesh.append(np.array((i,j,k),dtype=float))
    k_mesh = np.array(k_mesh) / np.array((mp1,mp2,mp3))
    return k_mesh


def path_to_plot(k_path,bvec):
    k_path_tmp = []
    x = 0.
    for kn in range(len(k_path)):
        k_path[kn] = np.dot(k_path[kn],bvec)
        if kn > 0:      x = x + np.linalg.norm(k_path[kn]-k_path[kn-1])
        k_path_tmp.append(x)
     
    return k_path_tmp


def write_output_eigk(eig_k,k_path):
    ofile = open('bands_python.dat','w')
    ofile.write('# Written on %d-%d-%d at %d:%d:%02d' %(dt.now().day,dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
    for n in range(eig_k.shape[1]):             # Loop over the band index
        for kn in range(eig_k.shape[0]):        # Loop over the k-points (of the path)
            ofile.write('\n%10.4f%10.4f' %(k_path[kn],eig_k[kn,n]))
        ofile.write('\n')
    ofile.close()



#-------------------------------------------------
#----------------- MAIN PROGRAM ------------------
#-------------------------------------------------

hr,Rvec,num_wann,nrtot = read_hr(prefix)
if ( num_wann != num_wann_occ+num_wann_emp ):   sys.exit('\nnum_wann not matching num_wann_occ+num_wann_emp\n')

bvec = np.zeros((3,3))      # primitive vectors of the reciprocal lattice (in 2*pi/alat)
bvec[0] = np.linalg.solve(avec,np.array(([1.,0.,0.])))
bvec[1] = np.linalg.solve(avec,np.array(([0.,1.,0.])))
bvec[2] = np.linalg.solve(avec,np.array(([0.,0.,1.])))

centers = read_centers(prefix,num_wann_occ,num_wann_emp,alat)

eig_k = np.zeros((len(k_path),num_wann))

for kn in range(len(k_path)):
    hk = np.zeros((num_wann,num_wann),dtype=complex)
    for iband in range(num_wann):
        for jband in range(num_wann):
            ws_dist = centers[jband] - centers[iband]
            for rn in range(nrtot):
                rvect = crys_to_cart(Rvec[rn],avec,+1)
                phase = corr_phase(kmesh[0],kmesh[1],kmesh[2],avec,k_path[kn],rvect,ws_dist)
                hk[iband,jband] = hk[iband,jband] + np.exp( 1.j * 2 * np.pi * np.dot(k_path[kn],Rvec[rn]) ) * phase * hr[rn,iband,jband]

    eig_k[kn] = np.linalg.eigvalsh(hk)

eig_k = eig_k
k_path = path_to_plot(k_path,bvec)
write_output_eigk(eig_k,k_path)

#-------------------------------------------------
#--------------- END MAIN PROGRAM ----------------
#-------------------------------------------------
