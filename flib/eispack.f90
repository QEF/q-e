#include "f_defs.h"
     SUBROUTINE RS(nm, n, a, w, matz, z, fv1, fv2, ierr)

        INTEGER n,nm,ierr,matz
        REAL*8 :: a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
        CHARACTER(LEN=1) :: jobz, uplo
        INTEGER            info, lda, lwork
        !verificare! cosi' non va bene
        REAL*8   work( 3*n-1 ), ztmp(nm,n)

        ztmp(1:N,1:N) = a(1:N,1:N)
        if ( matz == 0) then
            jobz = 'N'
        else
            jobz = 'V'
        end if

        lda = nm
        uplo = 'U'
        lwork = 3*n-1


        CALL DSYEV(jobz, uplo, n, ztmp, lda, w, work, lwork, info)


        ierr = info

        !non posso fare ierr = INFO  perche' non hanno gli stessi codici d'errore
        ! a meno che non ce ne freghi niente quando e' maggiore di 0 (perche' non si controlla il codice)

        !if (INFO.eq.0) then
        !    ierr = 0
        !end if

        if ( matz /= 0) then
          z(1:N,1:N) = ztmp(1:N,1:N)
        end if

        RETURN
      END SUBROUTINE
