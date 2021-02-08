import sys, os, jinja2

def render(tpl_path, context):
    path, filename = os.path.split(tpl_path)
    return jinja2.Environment(undefined=jinja2.StrictUndefined,
        loader=jinja2.FileSystemLoader(path or './')
    ).get_template(filename).render(context)


context = {'module_name': sys.argv[1], 'vars': [dict(zip(['name','type','dim'],x.split(':'))) for x in sys.argv[2:]]}

with open(sys.argv[1]+'_gpu.f90','w') as f:
    f.write(render('duplicated_module.jf90', context))
