#!/bin/bash

# This script allows to compile FoX with PGI v.19.10 Community Edition
# on Windows 10 - configure works in general but fails for FoX

set -x

rm -rf FoX

mkdir -p FoX
cd FoX

tar xvzf ../archive/fox.tgz
cd fox

myflags="-fast -Mcache_align -Mlarge_arrays -mp "
mydefs="-Mpreprocess -DPGF90 -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DFC_ABORT_ARG -DFC_EOR_LF"

mkdir -p objs/lib objs/finclude

cd fsys
pgfortran $myflags -c   $mydefs  fox_m_fsys_abort_flush.F90 -o fox_m_fsys_abort_flush.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_array_str.F90 -o fox_m_fsys_array_str.o
pgfortran $myflags -c     fox_m_fsys_realtypes.f90 -o fox_m_fsys_realtypes.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_format.F90 -o fox_m_fsys_format.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_parse_input.F90 -o fox_m_fsys_parse_input.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_count_parse_input.F90 -o fox_m_fsys_count_parse_input.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_string.F90 -o fox_m_fsys_string.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_string_list.F90 -o fox_m_fsys_string_list.o
pgfortran $myflags -c   $mydefs  fox_m_fsys_varstr.F90 -o fox_m_fsys_varstr.o
ar ruv libFoX_fsys.a fox_m_fsys_abort_flush.o fox_m_fsys_array_str.o fox_m_fsys_format.o fox_m_fsys_parse_input.o fox_m_fsys_count_parse_input.o fox_m_fsys_string.o fox_m_fsys_string_list.o fox_m_fsys_realtypes.o fox_m_fsys_varstr.o
cp -p libFoX_fsys.a ../objs/lib/libFoX_fsys.lib 
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd utils
pgfortran $myflags -c  -I../objs/finclude $mydefs  fox_m_utils_mtprng.F90 -o fox_m_utils_mtprng.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  fox_m_utils_uuid.F90 -o fox_m_utils_uuid.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  fox_m_utils_uri.F90 -o fox_m_utils_uri.o
pgfortran $myflags -c  -I../objs/finclude   FoX_utils.f90 -o FoX_utils.o
ar ruv libFoX_utils.a FoX_utils.o fox_m_utils_mtprng.o fox_m_utils_uuid.o fox_m_utils_uri.o
cp -p libFoX_utils.a ../objs/lib/libFoX_utils.lib 
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd common
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_charset.F90 -o m_common_charset.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_content_model.F90 -o m_common_content_model.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_error.F90 -o m_common_error.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_namecheck.F90  -o m_common_namecheck.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_element.F90 -o m_common_element.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_attrs.F90 -o m_common_attrs.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_buffer.F90 -o m_common_buffer.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_entities.F90 -o m_common_entities.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_notations.F90 -o m_common_notations.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_struct.F90 -o m_common_struct.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_namespaces.F90 -o m_common_namespaces.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_elstack.F90 -o m_common_elstack.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_io.F90 -o m_common_io.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  FoX_common.F90 -o FoX_common.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_common_entity_expand.F90 -o m_common_entity_expand.o
ar ruv libFoX_common.a m_common_attrs.o m_common_buffer.o m_common_charset.o m_common_namespaces.o m_common_error.o m_common_elstack.o m_common_io.o FoX_common.o m_common_namecheck.o m_common_entities.o m_common_notations.o m_common_element.o m_common_struct.o m_common_entity_expand.o m_common_content_model.o
cp -p libFoX_common.a ../objs/lib/libFoX_common.lib 
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd wxml
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wxml_escape.F90  -o m_wxml_escape.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wxml_core.F90 -o m_wxml_core.o
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wxml_overloads.F90 -o m_wxml_overloads.o
pgfortran $myflags -c  -I../objs/finclude   FoX_wxml.f90 -o FoX_wxml.o
ar  ruv libFoX_wxml.a m_wxml_escape.o m_wxml_core.o m_wxml_overloads.o FoX_wxml.o
cp -p libFoX_wxml.a ../objs/lib/libFoX_wxml.lib
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd wcml
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_stml.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_coma.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_metadata.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_core.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_geometry.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_lattice.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_lists.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_molecule.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_parameter.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_property.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wcml_inputdec.F90 
pgfortran $myflags -c  -I../objs/finclude   FoX_wcml.f90
ar ruv libFoX_wcml.a FoX_wcml.obj m_wcml_coma.obj m_wcml_core.obj m_wcml_stml.obj m_wcml_parameter.obj m_wcml_property.obj m_wcml_metadata.obj m_wcml_lattice.obj m_wcml_geometry.obj m_wcml_molecule.obj m_wcml_lists.obj m_wcml_inputdec.obj
cp -p libFoX_wcml.a ../objs/lib/libFoX_wcml.lib ;
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd wkml
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_color_def.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_color.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_contours.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_lowlevel.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_styling.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_core.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_chart.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_features.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_contours.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_wkml_coverage.F90 
pgfortran $myflags -c  -I../objs/finclude FoX_wkml.f90
ar  ruv libFoX_wkml.a  FoX_wkml.obj m_wkml_lowlevel.obj m_wkml_color.obj m_wkml_styling.obj m_wkml_features.obj m_wkml_coverage.obj m_wkml_core.obj m_wkml_contours.obj m_contours.obj m_wkml_color_def.obj m_wkml_chart.obj
cp -p libFoX_wkml.a ../objs/lib/libFoX_wkml.lib ;
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

cd sax
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_xml_source.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_reader.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_types.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_tokenizer.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_parser.F90 
pgfortran $myflags -c  -I../objs/finclude $mydefs  m_sax_operate.F90 
pgfortran $myflags -c  -I../objs/finclude   FoX_sax.f90
ar  ruv libFoX_sax.a m_sax_types.obj m_sax_tokenizer.obj m_sax_reader.obj m_sax_parser.obj m_sax_operate.obj m_sax_xml_source.obj FoX_sax.obj
cp -p libFoX_sax.a ../objs/lib/libFoX_sax.lib 
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..


cd dom
pgfortran $myflags -c  -I../objs/finclude m_dom_error.f90 -o m_dom_error.o
pgfortran $myflags -c  -I../objs/finclude $mydefs m_dom_dom.F90 -o m_dom_dom.o
pgfortran $myflags -c  -I../objs/finclude m_dom_parse.f90 -o m_dom_parse.o
pgfortran $myflags -c  -I../objs/finclude m_dom_utils.f90 -o m_dom_utils.o
pgfortran $myflags -c  -I../objs/finclude $mydefs m_dom_extras.F90 -o m_dom_extras.o
pgfortran $myflags -c  -I../objs/finclude FoX_dom.f90 -o FoX_dom.o
ar ruv libFoX_dom.a m_dom_error.o m_dom_parse.o m_dom_utils.o m_dom_extras.o m_dom_dom.o FoX_dom.o
cp -p libFoX_dom.a ../objs/lib/libFoX_dom.lib
for i in *.mod
do
   cp -p $i ../objs/finclude
done
cd ..

chmod 777 objs/lib/*
chmod 777 objs/finclude/*

cd ..

mkdir -p lib finclude bin
cp -p fox/objs/lib/* lib
cp -p fox/objs/finclude/* finclude

cd ..
