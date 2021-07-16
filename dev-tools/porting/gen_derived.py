# Usage: python ./gen_derived.py name_of_original_module
#
# a json file named 'name_of_original_module.json' must be present in the current directory.

import sys, os, json, jinja2

def render(tpl_path, context):
    path, filename = os.path.split(tpl_path)
    return jinja2.Environment(undefined=jinja2.StrictUndefined,
        loader=jinja2.FileSystemLoader(path or './')
    ).get_template(filename).render(context)


with open(sys.argv[1]+'.json','r') as f:
    context = json.load(f)

with open(sys.argv[1]+'_gpu.f90','w') as f:
    f.write(render('derived_type_duplicated_module.jf90', context))
