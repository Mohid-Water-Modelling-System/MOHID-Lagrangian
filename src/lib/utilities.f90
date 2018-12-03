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
    
    type :: utils_class
    contains
    procedure :: find_str
    procedure :: geo2m
    procedure :: m2geo_vec, m2geo_comp
    generic   :: m2geo => m2geo_vec, m2geo_comp
    procedure :: int2str
    procedure :: get_closest_twopow
    end type utils_class
    
    type(utils_class) :: Utils
    
    !public objects
    public :: Utils

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns the index of a given string in an array of strings. 
    !> Has optional mandatory flag.
    !> @param[in] self, str_array, str, mandatory
    !---------------------------------------------------------------------------
    integer function find_str(self, str_array, str, mandatory)
    class(utils_class), intent(in) :: self
    type(string), dimension(:), intent(in) :: str_array
    type(string), intent(in) :: str
    logical, optional, intent(in) :: mandatory
    type(string) :: outext
    do find_str=1, size(str_array)
        if (str == str_array(find_str)) return
    end do
    if(present(mandatory)) then
        if (mandatory) then
            outext = '[Utils::find_str]: string "'// str //'" not found on list, stopping'
            call Log%put(outext)
            stop
        end if
    end if
    end function find_str
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in meters given an array in
    !> geographical coordinates (lon, lat, z) and a lattitude
    !> @param[in] geovec, lat
    !---------------------------------------------------------------------------
    type(vector) function geo2m(self, geovec, lat) result(res)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: geovec
    real(prec), intent(in) :: lat
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    !pi = 3.1415926
    res = geovec
    res%x = res%x*R*cos(pi*lat/180.0)
    res%y = res%y*R
    end function geo2m

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters and a lattitude
    !> @param[in] self, mvec, lat
    !---------------------------------------------------------------------------
    type(vector) function m2geo_vec(self, mvec, lat) result(res)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: mvec    
    real(prec), intent(in) :: lat
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    !pi = 3.1415926
    res = mvec
    res%x = res%x/(R*cos(pi*lat/180.0))
    res%y = res%y/R
    end function m2geo_vec
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters and a lattitude
    !> @param[in] self, mvec, lat
    !---------------------------------------------------------------------------
    function m2geo_comp(self, mvec, lat, component)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: mvec    
    real(prec), dimension(:), intent(in) :: lat
    logical, intent(in) :: component
    real(prec), dimension(size(mvec)) :: m2geo_comp
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    m2geo_comp = mvec
    if (component) then
        m2geo_comp = m2geo_comp/R
    else
        m2geo_comp = m2geo_comp/(R*cos(pi*lat/180.0))
    end if
    end function m2geo_comp

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a zero paded string from an integer number
    !> and a format descriptor
    !> @param[in] self, fmt, i
    !---------------------------------------------------------------------------
    function int2str(self, fmt, i) result(res)
    class(utils_class), intent(in) :: self
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
    !> Returns the closest power of 2 or a given real number
    !> @param[in] self, num
    !---------------------------------------------------------------------------
    function get_closest_twopow(self, num) result(twopow)
    class(utils_class), intent(in) :: self
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
