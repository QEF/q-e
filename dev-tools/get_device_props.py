import sys
import subprocess as sp

have_yaml = True
try:
    import yaml
except:
    have_yaml = False
    
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO

print("\n\n This is a helper tool to check the details of your GPU before configuring QE.\n\n")
print("""Remeber to load CUDA environemt and run this on the COMPUTE NODE 
if you are not sure that the frontend node shares 
the same configuration as the backend nodes\n\n""")


compilation = sp.Popen("nvcc device_props.c -o device_props".split(), stdout=sp.PIPE)
streamdata = compilation.communicate()[0]
if compilation.returncode != 0:
    print("Compilation with nvcc failed. Did you load CUDA?")
    sys.exit()

compilation = sp.Popen("./device_props", stdout=sp.PIPE)
streamdata = compilation.communicate()[0]
if compilation.returncode != 0:
    print("\nDetails acquisition failed.")
else:
    yaml_data = StringIO(streamdata)
    conf_cc=""; conf_rt=""
    if have_yaml:
        data = yaml.load(yaml_data)
        for dev in data['devices']:
            print("Compute capabilities for dev {}: {}.{}".format(dev['devId'], dev['major'],dev['minor']))
            conf_cc=str(dev['major'])+str(5 if dev['minor']>=5 else 0)
        conf_rt = data['system']['runtimeVersion']
        conf_rt = str(int(conf_rt/1000))+'.'+str(int(conf_rt/10)%10)
    else:
        print("Yaml not available. Printing minimal info.\n")
        minor = ""; major = ""; devnum=0
        lines = yaml_data.readlines()
        for line in lines:
            if 'runtimeVersion' in line:
                _, conf_rt = line.split(':')
                conf_rt = conf_rt.strip()
                conf_rt = conf_rt[0]+'.'+conf_rt[2]
            if 'minor' in line:
                _, minor = line.split(':')
                minor = str(5 if int(minor)>=5 else 0)
            if 'major' in line:
                _, major = line.split(':')
            
            if minor != "" and major != "":
                print("Compute capabilities for dev {}: {}.{}".format(str(devnum), major.strip(),minor.strip()))
                conf_cc=major.strip()+minor.strip()
                minor = ""; major = ""; devnum += 1
                
    print("\n If all compute capabilities match, configure QE with:")
    print("./configure --with-cuda=$CUDA_HOME --with-cuda-cc={} --with-cuda-runtime={}\n".format(conf_cc, conf_rt))

    
