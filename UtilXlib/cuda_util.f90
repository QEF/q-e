!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contains a CUF KERNEL to perform the selected operation
!
MODULE cuda_util
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  PUBLIC :: cuf_memcpy, cuf_memset
  !
  INTERFACE cuf_memcpy
    MODULE PROCEDURE &
      cuf_memcpy_r1d, &
      cuf_memcpy_r2d, &
      cuf_memcpy_r3d, &
      cuf_memcpy_c1d, &
      cuf_memcpy_c2d, &
      cuf_memcpy_c3d
  END INTERFACE
  !
  INTERFACE cuf_memset
    MODULE PROCEDURE &
      cuf_memset_r1d, &
      cuf_memset_r2d, &
      cuf_memset_r3d, &
      cuf_memset_c1d, &
      cuf_memset_c2d, &
      cuf_memset_c3d
  END INTERFACE
  !
  CONTAINS
  !
  SUBROUTINE cuf_memcpy_r1d(array_out, array_in, range1 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN) :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = array_in(i1 )
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r1d
  !
  SUBROUTINE cuf_memcpy_r2d(array_out, array_in, range1, range2 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN) :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = array_in(i1,i2 )
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r2d
  !
  SUBROUTINE cuf_memcpy_r3d(array_out, array_in, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN) :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r3d
  !
  SUBROUTINE cuf_memcpy_c1d(array_out, array_in, range1 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN) :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = array_in(i1 )
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c1d
  !
  SUBROUTINE cuf_memcpy_c2d(array_out, array_in, range1, range2 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN) :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = array_in(i1,i2 )
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c2d
  !
  SUBROUTINE cuf_memcpy_c3d(array_out, array_in, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN) :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c3d
  !
  !
  SUBROUTINE cuf_memset_r1d(array_out, val, range1 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = val
    ENDDO
    !
  END SUBROUTINE cuf_memset_r1d
  !
  SUBROUTINE cuf_memset_r2d(array_out, val, range1, range2 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = val
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_r2d
  !
  SUBROUTINE cuf_memset_r3d(array_out, val, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = val
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_r3d
  !
  SUBROUTINE cuf_memset_c1d(array_out, val, range1 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = val
    ENDDO
    !
  END SUBROUTINE cuf_memset_c1d
  !
  SUBROUTINE cuf_memset_c2d(array_out, val, range1, range2 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = val
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_c2d
  !
  SUBROUTINE cuf_memset_c3d(array_out, val, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = val
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_c3d
  !
END MODULE cuda_util
!
!
! === TEMPLATE USED TO GENERATE THIS FILE ===
!
!    MODULE cuda_util
!      !
!      USE util_param,   ONLY : DP
!      !
!      IMPLICIT NONE
!      !
!      PUBLIC :: cuf_memcpy, cuf_memset
!      !
!      INTERFACE cuf_memcpy
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          cuf_memcpy_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!      END INTERFACE
!      !
!      INTERFACE cuf_memset
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          cuf_memset_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!      END INTERFACE
!      !
!      CONTAINS
!      !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE cuf_memcpy_{{t[0]|lower}}{{d}}d(array_out, array_in,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN) :: array_in({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(2)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_out, array_in
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: i{{dd+1}}, d{{dd+1}}s, d{{dd+1}}e
!    {%- endfor %}
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}s = range{{dd+1}}(1)
!        d{{dd+1}}e = range{{dd+1}}(2)
!    {%- endfor %}
!        !
!        !$cuf kernel do({{d}})
!    {%- for dd in range(d,0,-1) %}
!        DO i{{dd}} = d{{dd}}s, d{{dd}}e
!    {%- endfor %}
!           array_out( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} ) = array_in( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} )
!    {%- for dd in range(d) %}
!        ENDDO
!    {%- endfor %}
!        !
!      END SUBROUTINE cuf_memcpy_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!      !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE cuf_memset_{{t[0]|lower}}{{d}}d(array_out, val,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN) :: val
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(2)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_out
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: i{{dd+1}}, d{{dd+1}}s, d{{dd+1}}e
!    {%- endfor %}
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}s = range{{dd+1}}(1)
!        d{{dd+1}}e = range{{dd+1}}(2)
!    {%- endfor %}
!        !
!        !$cuf kernel do({{d}})
!    {%- for dd in range(d,0,-1) %}
!        DO i{{dd}} = d{{dd}}s, d{{dd}}e
!    {%- endfor %}
!           array_out( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} ) = val
!    {%- for dd in range(d) %}
!        ENDDO
!    {%- endfor %}
!        !
!      END SUBROUTINE cuf_memset_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!    END MODULE cuda_util
!
! === CODE TO GENERATE THE f90 FILE ===
!
!import sys, os, jinja2
!
!def render(tpl_path, context):
!    path, filename = os.path.split(tpl_path)
!    return jinja2.Environment(undefined=jinja2.StrictUndefined,
!        loader=jinja2.FileSystemLoader(path or './')
!    ).get_template(filename).render(context)
!with open('cuda_util.f90', 'w') as f: f.write(render('cuda_util.jf90', {'types': ['REAL', 'COMPLEX'], 'dimensions': 3}))
!