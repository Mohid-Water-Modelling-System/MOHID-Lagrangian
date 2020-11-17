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
    use datetime_module
    use simulationPrecision_mod
    use simulationLogger_mod

    implicit none
    private

    type :: utils_class
    contains
    procedure :: getDateFromISOString
    procedure :: getRelativeTimeFromString
    procedure :: getDateTimeFromDate
    procedure :: find_str
    procedure :: geo2cart_vec, geo2cart_comp, geo2cart_compSingle
    generic   :: geo2cart => geo2cart_vec, geo2cart_comp, geo2cart_compSingle
    procedure :: cart2geo_vec, cart2geo_comp
    generic   :: cart2geo => cart2geo_vec, cart2geo_comp
    procedure :: geo2m_vec, geo2m_comp, geo2m_compSingle, geo2m_vecFull
    generic   :: geo2m => geo2m_vec, geo2m_comp, geo2m_compSingle, geo2m_vecFull
    procedure :: m2geo_vec, m2geo_comp, m2geo_vecFull
    generic   :: m2geo => m2geo_vec, m2geo_comp, m2geo_vecFull
    procedure :: int2str
    procedure :: real2str
    procedure :: get_closest_twopow
    procedure :: isBoundedSingle, isBoundedArray
    generic :: isBounded => isBoundedSingle, isBoundedArray
    procedure :: appendArraysUniqueReal
    end type utils_class

    type(utils_class) :: Utils

    !public objects
    public :: Utils

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Function that returns a real with a time in seconds of a given date
    !> string, counting from a given date.
    !> @param[in] self, dateStr, RefDate
    !---------------------------------------------------------------------------
    real(prec) function getRelativeTimeFromString(self, dateStr, RefDate)
    class(utils_class), intent(in) :: self
    type(string), intent(in) :: dateStr
    type(datetime), intent(in) :: RefDate
    type(string), allocatable :: dc(:)
    integer, dimension(6) :: AbsoluteDateStr
    type(datetime) :: AbsoluteDate
    type(timedelta) :: delta

    call dateStr%split(tokens=dc, sep=' ')
    if (size(dc) > 1) then
        AbsoluteDateStr = self%getDateFromISOString(dateStr)
        AbsoluteDate = self%getDateTimeFromDate(AbsoluteDateStr)
        delta = AbsoluteDate - RefDate
        getRelativeTimeFromString = delta%total_seconds()
    else
        getRelativeTimeFromString = dateStr%to_number(kind=1._R8P)
    end if

    end function getRelativeTimeFromString

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Function that returns a real with a time in seconds of a given date
    !> string, counting from a given date.
    !> @param[in] self, dateArray
    !---------------------------------------------------------------------------
    type(datetime) function getDateTimeFromDate(self, dateArray)
    class(utils_class), intent(in) :: self
    integer, dimension(6), intent(in) :: dateArray
    type(string) :: outext
    if (size(dateArray) == 6) then
        getDateTimeFromDate = datetime(dateArray(1), dateArray(2), dateArray(3), dateArray(4), dateArray(5), dateArray(6))
    else
        outext = '[Utils::getDateTimeFromDate] Array is not the correct size, check iso date specs. Stopping'
        call Log%put(outext)
        stop
    end if
    end function getDateTimeFromDate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Function that returns an integer array of type (year, month, day, hour,
    !> minute, second) from an ISO date string.
    !> @param[in] self, dateStr
    !---------------------------------------------------------------------------
    function getDateFromISOString(self, dateStr)
    class(utils_class), intent(in) :: self
    type(string), intent(in) :: dateStr
    integer, dimension(6) :: getDateFromISOString
    integer :: i
    type(string), allocatable :: dc(:)
    type(string) :: outext
    call dateStr%split(tokens=dc, sep=' ')
    if (size(dc) == 6) then
        do i=1, size(dc)
            getDateFromISOString(i) = dc(i)%to_number(kind=1._R8P)
        end do
    else
        outext = '[Utils::getDateFromISOString] Date '// dateStr //' not in correct format. Eg. "2009 03 01 00 00 00"'
        call Log%put(outext)
        stop
    end if
    end function getDateFromISOString

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
    find_str = MV_INT
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
    !> geographical coordinates (lon, lat, z), using a refecence point
    !> @param[in] self, geovec
    !---------------------------------------------------------------------------
    elemental type(vector) function geo2cart_vec(self, geovec, ref)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: geovec
    type(vector), intent(in) :: ref
    integer :: R
    real(prec) :: pi
    pi = 4*atan(1.0)
    R = 6378000 !earth radius in meters
    geo2cart_vec%x = (geovec%x - ref%x)*(R*pi/180.0)*cos(pi*geovec%y/180.0)
    geo2cart_vec%y = (geovec%y - ref%y)*(R*pi/180.0)
    end function geo2cart_vec
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns an array of coordinates converted to degrees, using a refecence 
    !> point
    !> @param[in] self, geovec, lat, isLat
    !---------------------------------------------------------------------------
    function geo2cart_comp(self, geovec, lat, ref, isLat)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: geovec
    real(prec), dimension(:), intent(in) :: lat
    type(vector), intent(in) :: ref
    logical, intent(in) :: isLat
    real(prec), dimension(size(geovec)) :: geo2cart_comp
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378000 !earth radius in meters
    if (isLat) then
        geo2cart_comp = (geovec - ref%y)*(R*pi/180.0)
    else
        geo2cart_comp = (geovec - ref%x)*((R*pi/180.0)*cos(pi*lat/180.0))
    end if
    end function geo2cart_comp
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a coordinate converted to degrees, using a refecence 
    !> point
    !> @param[in] self, geovec, lat, isLat
    !---------------------------------------------------------------------------
    function geo2cart_compSingle(self, geovec, lat, ref, isLat)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: geovec
    real(prec), intent(in) :: lat
    type(vector), intent(in) :: ref
    logical, intent(in) :: isLat
    real(prec), dimension(size(geovec)) :: geo2cart_compSingle
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378000 !earth radius in meters
    if (isLat) then
        geo2cart_compSingle = (geovec - ref%y)*(R*pi/180.0)
    else
        geo2cart_compSingle = (geovec - ref%x)*((R*pi/180.0)*cos(pi*lat/180.0))
    end if
    end function geo2cart_compSingle
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters, using a refecence 
    !> point
    !> @param[in] self, mvec
    !---------------------------------------------------------------------------
    elemental type(vector) function cart2geo_vec(self, mvec, ref)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: mvec
    type(vector), intent(in) :: ref
    real(prec) :: R
    real(prec) :: pi
    pi = 4*atan(1.0)
    R = 6378000 !earth radius in meters
    cart2geo_vec%y = mvec%y/(R*pi/180.0) + ref%y
    cart2geo_vec%x = mvec%x/((R*pi/180.0)*cos(pi*mvec%y/180.0)) + ref%x
    end function cart2geo_vec
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters, using a refecence 
    !> point
    !> @param[in] self, mvec, lat, component
    !---------------------------------------------------------------------------
    function cart2geo_comp(self, mvec, lat, ref, isLat)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: mvec
    real(prec), dimension(:), intent(in) :: lat
    type(vector), intent(in) :: ref
    logical, intent(in) :: isLat
    real(prec), dimension(size(mvec)) :: cart2geo_comp
    real(prec) :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378000 !earth radius in meters
    if (isLat) then
        cart2geo_comp = mvec/(R*pi/180.0)+ ref%y
    else
        cart2geo_comp = mvec/((R*pi/180.0)*cos(pi*lat/180.0)) + ref%x
    end if
    end function cart2geo_comp
    
    
    
    

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in meters given an array in
    !> geographical coordinates (lon, lat, z) and a lattitude
    !> @param[in] self, geovec, lat
    !---------------------------------------------------------------------------
    type(vector) function geo2m_vec(self, geovec, lat) result(res)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: geovec
    real(prec), intent(in) :: lat
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    !pi = 3.1415926
    res = geovec
    res%x = res%x*(R*pi/180.0)*cos(pi*lat/180.0)
    res%y = res%y*(R*pi/180.0)
    end function geo2m_vec
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in meters given an array in
    !> geographical coordinates (lon, lat, z) and a lattitude
    !> @param[in] self, geovec
    !---------------------------------------------------------------------------
    elemental type(vector) function geo2m_vecFull(self, geovec)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: geovec
    integer :: R
    real(prec) :: pi
    pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    !pi = 3.1415926
    geo2m_vecFull = geovec
    geo2m_vecFull%x = geo2m_vecFull%x*(R*pi/180.0)*cos(pi*geo2m_vecFull%y/180.0)
    geo2m_vecFull%y = geo2m_vecFull%y*(R*pi/180.0)
    end function geo2m_vecFull    

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns an array of coordinates converted to degrees. 
    !> @param[in] self, geovec, lat, isLat
    !---------------------------------------------------------------------------
    function geo2m_comp(self, geovec, lat, isLat)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: geovec
    real(prec), dimension(:), intent(in) :: lat
    logical, intent(in) :: isLat
    real(prec), dimension(size(geovec)) :: geo2m_comp
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    if (isLat) then
        geo2m_comp = geovec*(R*pi/180.0)
    else
        geo2m_comp = geovec*((R*pi/180.0)*cos(pi*lat/180.0))
    end if
    end function geo2m_comp
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns an array of coordinates converted to degrees. 
    !> @param[in] self, geovec, lat, isLat
    !---------------------------------------------------------------------------
    function geo2m_compSingle(self, geovec, lat, isLat)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: geovec
    real(prec), intent(in) :: lat
    logical, intent(in) :: isLat
    real(prec), dimension(size(geovec)) :: geo2m_compSingle
    integer :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137 !earth radius in meters
    if (isLat) then
        geo2m_compSingle = geovec*(R*pi/180.0)
    else
        geo2m_compSingle = geovec*((R*pi/180.0)*cos(pi*lat/180.0))
    end if
    end function geo2m_compSingle
    
    
    
    

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
    real(prec) :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137.0 !earth radius in meters
    res = mvec
    res%y = res%y/(R*pi/180.0)
    res%x = res%x/((R*pi/180.0)*cos(pi*lat/180.0))
    end function m2geo_vec

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters and a lattitude
    !> @param[in] self, mvec
    !---------------------------------------------------------------------------
    elemental type(vector) function m2geo_vecFull(self, mvec)
    class(utils_class), intent(in) :: self
    type(vector), intent(in) :: mvec
    real(prec) :: R
    real(prec) :: pi
    pi = 4*atan(1.0)
    R = 6378137.0 !earth radius in meters
    m2geo_vecFull=mvec
    m2geo_vecFull%y = mvec%y/(R*pi/180.0)
    m2geo_vecFull%x = mvec%x/((R*pi/180.0)*cos(pi*m2geo_vecFull%y/180.0))
    end function m2geo_vecFull
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a vector in geographical coordinates
    !> (lon, lat, z) given an array in meters and a lattitude
    !> @param[in] self, mvec, lat, component
    !---------------------------------------------------------------------------
    function m2geo_comp(self, mvec, lat, component)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: mvec
    real(prec), dimension(:), intent(in) :: lat
    logical, intent(in) :: component
    real(prec), dimension(size(mvec)) :: m2geo_comp
    real(prec) :: R
    real(prec) :: pi = 4*atan(1.0)
    R = 6378137.0 !earth radius in meters
    m2geo_comp = mvec
    if (component) then
        m2geo_comp = m2geo_comp/(R*pi/180.0)
    else
        m2geo_comp = m2geo_comp/((R*pi/180.0)*cos(pi*lat/180.0))
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
    !> Returns a string from an real number
    !> and a format descriptor
    !> @param[in] self, fmt, i
    !---------------------------------------------------------------------------
    function real2str(self, fmt, i) result(res)
    class(utils_class), intent(in) :: self
    character(:), allocatable :: res
    character(len=6), intent(in) :: fmt ! format descriptor
    real(prec), intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp, fmt) i
    res = trim(tmp)
    end function real2str

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

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Logical function that checks if a number is bounded between 2 values
    !> @param[in] self, nums, minBound, maxBound, eta
    !---------------------------------------------------------------------------
    logical function isBoundedSingle(self, nums, minBound, maxBound, eta)
    class(utils_class), intent(in) :: self
    real(prec), intent(in) :: nums
    real(prec), intent(in) :: minBound
    real(prec), intent(in) :: maxBound
    real(prec), optional, intent(in) :: eta
    real(prec) :: ieta
    ieta = 0.0
    if (present(eta)) ieta = eta
    isBoundedSingle = (nums >= minBound).and.(nums <= maxBound+ieta)
    end function isBoundedSingle

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Logical function that checks if a set of numbers are bounded between 2 values
    !> @param[in] self, nums, minBound, maxBound, eta
    !---------------------------------------------------------------------------
    logical function isBoundedArray(self, nums, minBound, maxBound, eta)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), intent(in) :: nums
    real(prec), intent(in) :: minBound
    real(prec), intent(in) :: maxBound
    real(prec), optional, intent(in) :: eta
    real(prec) :: ieta
    ieta = 0.0
    if (present(eta)) ieta = eta
    isBoundedArray = all(nums >= minBound).and.all(nums <= maxBound+ieta)
    end function isBoundedArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> receives 2 real arrays and appends the new values of the second to the first.
    !> also returns an optional logical array with the used values of the second array
    !> @param[in] self, baseArray, newArray, usedArray
    !---------------------------------------------------------------------------
    subroutine appendArraysUniqueReal(self, baseArray, newArray, usedArray)
    class(utils_class), intent(in) :: self
    real(prec), dimension(:), allocatable, intent(inout) :: baseArray
    real(prec), dimension(:), allocatable, intent(in) :: newArray
    logical, dimension(:), allocatable, intent(out), optional :: usedArray
    logical, dimension(:), allocatable :: usedTemp
    real(prec), dimension(:), allocatable :: newBaseArray
    integer :: i, j, k

    allocate(usedTemp(size(newArray)))
    usedTemp = .false.
    do i=1, size(newArray)
        if (.not.any(baseArray == newArray(i))) then
            usedTemp(i) = .true.
        end if
    end do
    allocate(newBaseArray(size(baseArray) + count(usedTemp)))
    newBaseArray(1:size(baseArray)) = baseArray
    k = 1
    do i= size(baseArray)+1, size(baseArray) + count(usedTemp)
        do j=k, size(newArray)
            if (usedTemp(j)) then
                newBaseArray(i) =  newArray(j)
                k = k+1
                exit
            end if
            k = k+1
        end do
    end do
    deallocate(baseArray)
    allocate(baseArray, source=newBaseArray)

    if (present(usedArray)) allocate(usedArray, source = usedTemp)
    end subroutine appendArraysUniqueReal

    end module utilities_mod