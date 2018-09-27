    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : utilities
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that provides useful back-end routines.
    !------------------------------------------------------------------------------

    module utilities_mod

    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    use simulation_precision_mod
    use simulation_logger_mod

    implicit none
    private

    !public routines
    public :: get_closest_twopow, int2str, geo2m, m2geo

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public function that returns a vector in meters given an array in
    !> geographical coordinates (lon, lat, z) and a lattitude
    !> @param[in] geovec, lat
    !---------------------------------------------------------------------------
    type(vector) function geo2m(geovec, lat) result(res)
    type(vector), intent(in) :: geovec
    real(prec), intent(in) :: lat
    integer :: R
    real(prec) :: pi
    R = 6378137 !earth radius in meters
    pi = 3.1415926
    res = geovec
    res%x = res%x*R*cos(pi*lat/180.0)
    res%y = res%y*R
    end function geo2m

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public function that returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters and a lattitude
    !> @param[in] mvec, lat
    !---------------------------------------------------------------------------
    type(vector) function m2geo(mvec, lat) result(res)
    type(vector), intent(in) :: mvec
    real(prec), intent(in) :: lat
    integer :: R
    real(prec) :: pi
    R = 6378137 !earth radius in meters
    pi = 3.1415926
    res = mvec
    res%x = res%x/(R*cos(pi*lat/180.0))
    res%y = res%y/R
    end function m2geo

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public function that returns a zero paded string from an integer number
    !> and a format descriptor
    !> @param[in] fmt, i
    !---------------------------------------------------------------------------
    function int2str(fmt, i) result(res)
    character(:), allocatable :: res
    character(len=6), intent(in) :: fmt ! format descriptor
    integer, intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp, fmt) i
    res = trim(tmp)
    end function

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public function that returns the closest power of 2 or a given real number
    !> @param[in] num
    !---------------------------------------------------------------------------
    function get_closest_twopow(num) result(twopow)
    implicit none
    real(prec), intent(in) :: num
    real(prec) :: twopow
    integer :: i
    real(prec) :: dist1, dist2
    do i=-4, 10
        twopow = 2.0**i
        if (num < twopow) then
            dist1 = sqrt(twopow-num)
            dist2 = sqrt(num-2.0**(i-1))
            if (dist2 < dist1) then
                twopow = 2.0**(i-1)
                exit
            endif
            exit
        endif
    enddo
    end function get_closest_twopow

    end module utilities_mod
