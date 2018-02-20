# real and complex never tested! Only 8 lengh 
types=['INTEGER', 'REAL', 'REAL(8)', 'COMPLEX', 'COMPLEX(8)']
typeconv={'INTEGER': 'INT', 'REAL': 'REAL', 'REAL(8)': 'DBLE', 'COMPLEX': 'CMPLX', 'COMPLEX(8)': 'DCMPLX'}

datasize='10'
ranks  ={'1': '', 'v': '(datasize)', 'm': '(datasize,datasize)', 't': '(datasize,datasize,datasize)'}
nextr={'1': 'v', 'v': 'm', 'm': 't', 't': ''}

import os, re

mkfile = open('./autotest.inc', 'w')
mkfile.write('override SRCS += ')

for file in os.listdir("."):
  if file.endswith(".tmpl"):
    
    with open(file,'r') as f:
      data = f.read()
      
      # ! Implemented: i1, iv, rv, iv, rm, im, cv
      mtch = re.findall(r'Implemented\s*:\s*((\w\w|\*)(,\s*\w\w)*)',data)
      if mtch:
          impl_interfaces = mtch[0][0].replace(' ','').split(',')
      else:
          impl_interfaces = ['*']
          
      for t in types:
        for k,s in ranks.items():
          tmp = os.path.splitext(os.path.basename(file))[0]
          tname = t[0].lower()+k
          ofile = 'auto' + tmp + '_' + tname + '.f90'
          
          if not ( ('*' in impl_interfaces) or (tname in impl_interfaces) ):
              continue
          
          with open(ofile,'w') as fo:
            if k == '1':
              allf = ''; sumf = ''
            else:
              allf = 'ALL'; sumf = 'SUM'

            fo.write(data.format(  vname=(t[0].lower()+k), \
                                    datasize=datasize, \
                                    type=t, size=s, \
                                    all=s.replace('datasize',':'), \
                                    sizep1=ranks.get(nextr[k],'invalid'), \
                                    allp1=ranks.get(nextr[k],'invalid').replace('datasize',':'), \
                                    allf=allf, sumf=sumf, typeconv=typeconv[t]))
          mkfile.write(ofile + '\\\n')
