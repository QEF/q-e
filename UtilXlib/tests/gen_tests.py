# real and complex never tested! Only 8 lengh 
types=['INTEGER', 'REAL', 'REAL(8)', 'COMPLEX', 'COMPLEX(8)']
typeconv={'INTEGER': 'INT', 'REAL': 'REAL', 'REAL(8)': 'DBLE', 'COMPLEX': 'CMPLX', 'COMPLEX(8)': 'DCMPLX'}

datasize='10'
ranks  ={'1': '', 'v': '(datasize)', 'm': '(datasize,datasize)', 't': '(datasize,datasize,datasize)'}
alldata={'1': '', 'v': '(:)', 'm': '(:,:)', 't': '(:,:,:)'}
import os
for file in os.listdir("."):
  if file.endswith(".tmpl"):
    
    with open(file,'r') as f:
      data = f.read()
      for t in types:
        for k,s in ranks.items():
          tmp = os.path.splitext(os.path.basename(file))[0]
          ofile = 'auto' + tmp +'_'+t[0].lower()+k+'.f90'
          with open(ofile,'w') as fo:
            if k == '1':
              allf = ''; sumf = ''
            else:
              allf = 'ALL'; sumf = 'SUM'
            
            fo.write(data.format(  vname=(t[0].lower()+k), \
                                    datasize=datasize, \
                                    type=t, size=s, \
                                    all=alldata[k], \
                                    allf=allf, sumf=sumf, typeconv=typeconv[t]))
          print (ofile + '\\')
