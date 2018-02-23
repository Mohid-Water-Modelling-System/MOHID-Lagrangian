module tracer_precision

    implicit none 

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.d0)

    ! Precision used here
    integer,  parameter :: prec      = sp 
    integer,  parameter :: prec_time = sp 
    integer,  parameter :: prec_wrt  = sp 

    ! Missing value aliases 
    real(prec), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(prec), parameter :: MV     = MISSING_VALUE_DEFAULT
    real(prec), parameter :: MV_INT = int(MISSING_VALUE_DEFAULT)

    ! Error values
    real(prec), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter   :: ERR_IND  = -1 
    
end module tracer_precision 