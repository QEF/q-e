types=['INTEGER', 'REAL', 'REAL(8)', 'COMPLEX', 'COMPLEX(8)']
ranks={'1': '', 'v': '(10)', 'm': '(10,10)', 't': '(10,10,10)'}

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
            allf = 'ALL'
            if k == '1':
              allf = ''
            fo.write(data.format(vname=(t[0].lower()+k), type=t, size=s, all=s.replace('10',':'), allf=allf))
          print (ofile + '\\')
