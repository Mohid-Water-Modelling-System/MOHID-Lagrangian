!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracer_precision
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module to control the precision of the variables trough the project.
!------------------------------------------------------------------------------
    
module tracer_precision

    implicit none
    private
    public :: prec, prec_time, prec_wrt
    public :: MISSING_VALUE_DEFAULT, MV, MV_INT
    public :: ERR_DIST, ERR_IND

    integer,  parameter :: sp  = kind(1.0)  !> Simple precision definition switch
    integer,  parameter :: dp  = kind(1.d0) !> Double precision definition switch

    ! Precision used throughout is define here. Change at will.
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