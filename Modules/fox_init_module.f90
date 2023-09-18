MODULE fox_init_module
#if defined(__fox)
USE mp,   ONLY: mp_bcast, mp_barrier 
USE io_global, ONLY: ionode, ionode_id
USE mp_images, ONLY: intra_image_comm
USE m_common_io, ONLY: setup_io, io_err, io_eor, io_eof
#endif
IMPLICIT NONE 
PRIVATE
PUBLIC     :: fox_init 

CONTAINS 
   SUBROUTINE fox_init() 
      INTEGER   :: errcodes(3)
#if defined(__fox)
      IF (ionode) THEN
         call setup_io() 
         errcodes(1) = io_err
         errcodes(2) = io_eor
         errcodes(3) = io_eof
      END IF
      CALL mp_barrier(intra_image_comm) 
      CALL mp_bcast(errcodes, ionode_id, intra_image_comm) 
      CALL setup_io(ERR_CODE = errcodes(1), EOR_CODE = errcodes(2), EOF_CODE = errcodes(3))
#else
         errcodes(:) = 0
#endif
   END SUBROUTINE fox_init
END MODULE fox_init_module
