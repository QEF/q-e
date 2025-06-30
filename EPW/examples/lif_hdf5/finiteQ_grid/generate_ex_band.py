import os
import shutil

def read_klist(filename):
    klist = []
    file = open(filename,'r')

    for line in file:
        if 'K_POINTS' in line:
            numK = int(file.readline().strip())
            for ik in range(numK):
                temp = file.readline().split()
                kpts = []
                for ix in range(3):
                    kpts.append(float(temp[ix]))
                klist.append(kpts)
            break

    file.close()
    
    return klist
        
def generate_inputs(kpt, filepath):
    towrite = 'exciton_Q_shift 2 ';
    for i in range(3):
        towrite = towrite + '{:15.10f}'.format(-kpt[i]);
    
    print(towrite)
    file = open(filepath + '13-kernel/kernel.inp','a')
    file.write(towrite + '\n');
    file.close();
    
    file = open(filepath + '14-absorption/absorption.inp','a')
    file.write(towrite + '\n');
    file.close();

def generate_all_input(klist, seed):
    numK = len(klist);
    for i in range(numK):
        filepath = './Q' + str(i + 1) + '/';
        os.makedirs(filepath, exist_ok=True)
        os.system('cp -r ' + seed + '/13-kernel ' + filepath + '13-kernel/')
        os.system('cp -r ' + seed + '/14-absorption ' + filepath + '14-absorption/')
        generate_inputs(klist[i][:], filepath)
        os.system('cp link.sh ' + filepath)
        #os.system('bash ' + filepath + '/link.sh ')
        
if __name__ == '__main__':
    klist = read_klist('../espresso/05-Wfn_fi/wfn.in') 
    print(klist)
    generate_all_input(klist, 'seed')

