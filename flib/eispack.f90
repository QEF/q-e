     SUBROUTINE RS(nm, n, a, w, matz, z, fv1, fv2, ierr)

        INTEGER n,nm,ierr,matz
        DOUBLE PRECISION a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
        CHARACTER          jobz, uplo
        INTEGER            info, lda, lwork
        !verificare! cosi' non va bene
        DOUBLE PRECISION   work( 3*n-1 )

        z = a
        if (matz.eq.0) then
            jobz = 'N'
        else
            jobz = 'V'
        end if

        lda = nm
        uplo = 'U'
        lwork = 3*n-1

        CALL DSYEV(jobz, uplo, n, z, lda, w, work, lwork, info)

        ierr = info

        !non posso fare ierr = INFO  perche' non hanno gli stessi codici d'errore
        ! a meno che non ce ne freghi niente quando e' maggiore di 0 (perche' non si controlla il codice)

        !if (INFO.eq.0) then
        !    ierr = 0
        !end if

        RETURN
      END SUBROUTINE

