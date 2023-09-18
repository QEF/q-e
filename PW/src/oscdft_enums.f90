MODULE oscdft_enums
#if defined (__OSCDFT)
   IMPLICIT NONE
   SAVE
   INTEGER, PARAMETER :: CONV_MULTIPLIERS = 0,&
                         CONV_GRADIENT = 1,&
                         CONV_ENERGY = 2,&
                         CONV_ALWAYS_FALSE = -1,&
                         CONV_ALWAYS_TRUE = -2,&
                         OPT_GRADIENT_DESCENT = 0,&
                         OPT_GRADIENT_DESCENT2 = 1,&
                         CONV_FUNC_MAXVAL = 0,&
                         CONV_FUNC_NORM = 1,&
                         CONV_FUNC_NORM_AVERAGE = 2,&
                         CONSTR_FALSE = 0,&
                         CONSTR_TRUE = 1,&
                         CONSTR_LE = 2,&
                         CONSTR_GE = 3,&
                         CONSTR_LE2 = 4,&
                         CONSTR_GE2 = 5,&
                         CONSTR_LE3 = 6,&
                         CONSTR_GE3 = 7,&
                         CONSTR_D1 = 8,&
                         ITER_MULTIPLIERS_RHO = 0,&
                         ITER_RHO_MULTIPLIERS = 1,&
                         OCCUP_TRACE = -1,&
                         OCCUP_SUM = -2,&
                         OSCDFT_NONE = 0,&
                         ! OSCDFT_FREEZE = 1,&
                         OSCDFT_PERMUTE = 2,&
                         OSCDFT_DUMMY = 99
#endif
END MODULE oscdft_enums
