!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : simulation_precision
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

module simulation_precision_mod

    use penf

    implicit none
    private
    public :: prec, prec_time, prec_wrt
    public :: MISSING_VALUE_DEFAULT, MV, MV_INT
    public :: ERR_DIST, ERR_IND
    public :: CHAR_LEN

    integer,  parameter :: sps  = kind(1._R4P)   !< Simple precision definition switch
    integer,  parameter :: dps  = kind(1._R8P)   !< Double precision definition switch

    ! Precision used throughout is define here. Change at will.
    integer,  parameter :: prec      = sps
    integer,  parameter :: prec_time = sps
    integer,  parameter :: prec_wrt  = sps

    ! Missing value aliases
    real(prec), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dps
    real(prec), parameter :: MV     = MISSING_VALUE_DEFAULT
    real(prec), parameter :: MV_INT = int(MISSING_VALUE_DEFAULT)

    ! Error values
    real(prec), parameter :: ERR_DIST = 1E8_dps
    integer,  parameter   :: ERR_IND  = -1
    
    ! char handling
    integer, parameter :: CHAR_LEN = 99

end module simulation_precision_mod
