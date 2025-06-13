    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : geometry
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines geometry classes and related methods.
    !------------------------------------------------------------------------------

    module geometry_mod

    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    use ModuleEnterData
    use ModuleGlobalData

    use simulationPrecision_mod
    use simulationLogger_mod
    !use simulationGlobals_mod
    use utilities_mod
    use xmlParser_mod

    implicit none
    private

    !Update with any new additions as they are added
    type :: geometry_class
        type(string), allocatable, dimension(:) :: list !< String list with the name of possible geometry types.
    contains
    procedure :: initialize => allocatelist !<Builds the geometry list, possible geometry types (new types must be manually added)
    procedure :: inlist                     !<checks if a given geometry is defined as a derived type (new types must be manually added)
    procedure :: allocateShape              !< Returns allocated Shape with a specific shape given a name
    procedure :: fillSize                   !<Gets the number of points that fill a geometry
    procedure, public :: getFillPoints      !< returns a list of points that fill a given shape
    procedure :: getCenter                  !<Function that retuns the shape baricenter
    procedure :: getPoints                  !<Function that retuns the points (vertexes) that define the geometrical shape
    procedure :: getnumPoints               !<Function that returns the number of points (vertexes) that define the geometrical shape
    procedure, public :: setPolygon         !<assembles a polygon from a file
    procedure, public :: isPolygon          !<returns true if type is polygon
    procedure :: print => printGeometry     !<prints the geometry type and contents
    end type geometry_class

    type :: shape                      !<Type - extendable shape class
        type(vector) :: pt              !< Coordinates of a point
    end type

    type, extends(shape) :: point   !<Type - point class
    end type

    type, extends(shape) :: line    !<Type - line class
        type(vector) :: last            !< Coordinates of the end point
    end type
    
    type, extends(shape) :: polyline    !<Type - line class
        type(vector), dimension(:), allocatable :: point            !< point array
    end type

    type, extends(shape) :: sphere  !<Type - sphere class
        real(prec) :: radius            !< Sphere radius (m)
    end type

    type, extends(shape) :: box     !<Type - box class
        type(vector) :: size            !< Box size (m)
    end type

    type, extends(shape) :: polygon     !<Type - polygon class
        type(vector), dimension(:), allocatable :: vertex            !< vertex array
        type(vector) :: bbMin, bbMax
    contains
    procedure, public :: setBoundingBox
    procedure, public :: getMetricPolygon
    end type
    
    type(string) :: notRead

    type(geometry_class) :: Geometry

    public :: shape, point, line, sphere, box, polygon, polyline
    public :: Geometry

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns allocated Shape with a specific shape given a name
    !> @param[in] self, name, shape_obj
    !---------------------------------------------------------------------------
    subroutine allocateShape(self, name, shape_obj)
    class(geometry_class), intent(in) :: self
    type(string), intent(in) :: name
    class(shape), allocatable, intent(inout) :: shape_obj
    type(string) :: outext
    select case (name%chars())
    case ('point')
        allocate(point::shape_obj)
    case ('sphere')
        allocate(sphere::shape_obj)
    case ('box')
        allocate(box::shape_obj)
    case ('line')
        allocate(line::shape_obj)
    case ('polygon')
        allocate(polygon::shape_obj)
    case ('polyline')
        allocate(polyline::shape_obj)
        case default
        outext='[Geometry::allocateShape]: unexpected type for geometry object "'//name//'", stoping'
        call Log%put(outext)
        stop
    end select
    end subroutine allocateShape

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public routine to allocate the possible geometry name list
    !---------------------------------------------------------------------------
    subroutine allocatelist(self)
    class(geometry_class), intent(inout) :: self
    allocate(self%list(6))
    self%list(1) ='point'
    self%list(2) ='line'
    self%list(3) ='box'
    self%list(4) ='sphere'
    self%list(5) ='polygon'
    self%list(6) ='polyline'
    
    notRead = 'notRead'
    end subroutine allocatelist

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public function that returns a logical if the input geometry name is valid
    !> @param[in] self, geomname
    !---------------------------------------------------------------------------
    logical function inlist(self, geomname) result(tf)
    implicit none
    class(geometry_class), intent(in) :: self
    type(string), intent(in) :: geomname
    integer :: i
    tf = .false.
    do i=1, size(self%list)
        if (geomname == self%list(i)) then
            tf = .true.
        endif
    enddo
    end function inlist
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Public function that returns a logical if the input geometry name is a polygon
    !> @param[in] self, geomname
    !---------------------------------------------------------------------------
    logical function isPolygon(self, geomname) result(tf)
    implicit none
    class(geometry_class), intent(in) :: self
    type(string), intent(in) :: geomname
    integer :: i
    tf = .false.

    if (geomname == 'polygon') then
        tf = .true.
    endif
    end function isPolygon

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to get the number of points that fill a given geometry
    !> @param[in] self, shapetype, dp
    !---------------------------------------------------------------------------
    function fillSize(self, shapetype, dp)
    class(geometry_class), intent(in) :: self
    class(shape), intent(in) :: shapetype
    type(vector), intent(in) :: dp
    integer :: fillSize
    type(vector) :: temp
    type(string) :: outext
    select type (shapetype)
    type is (shape)
    class is (box)
        fillSize = max((int(shapetype%size%x/dp%x)+1)*(int(shapetype%size%y/dp%y)+1)*(int(shapetype%size%z/dp%z)+1),1)
    class is (point)
        fillSize = 1
    class is (line)
        fillSize = lineNpCount(dp, shapetype%pt, shapetype%last)
    class is (polyline)
        fillSize = polylineNpCount(dp, shapetype)
    class is (sphere)
        fillSize = sphere_np_count(dp, shapetype%radius)
    class is (polygon)
        fillSize = polygonNpCount(dp, shapetype)
        class default
        outext='[geometry::fillSize] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    end function fillSize

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns an array of points that fill a given shape
    !> @param[in] self, shapetype, dp
    !---------------------------------------------------------------------------
    function getFillPoints(self, shapetype, dp)
    class(geometry_class), intent(in) :: self
    class(shape) :: shapetype
    type(vector), intent(in) :: dp
    type(vector), dimension(:), allocatable :: getFillPoints
    integer :: fillSize
    type(string) :: outext

    fillSize = self%fillSize(shapetype, dp)
    allocate(getFillPoints(fillSize))
    select type (shapetype)
    type is (shape)
    class is (box)
        call box_grid(dp, shapetype%size, fillSize, getFillPoints)
    class is (point)
        getFillPoints(1)=0.0
    class is (line)
        call line_grid(Utils%geo2m(shapetype%last-shapetype%pt, shapetype%pt%y), fillSize, getFillPoints)
    class is (polyline)
        call polylineGrid(dp, shapetype, getFillPoints)
    class is (sphere)
        call sphere_grid(dp, shapetype%pt, shapetype%radius, fillSize, getFillPoints)
    class is (polygon)
        call polygonGrid(dp, shapetype, getFillPoints)
        class default
        outext='[geometry::getFillPoints] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    end function getFillPoints

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to get the baricenter of a given geometry
    !> @param[in] self, shapetype
    !---------------------------------------------------------------------------
    function getCenter(self, shapetype) result(center)
    class(geometry_class), intent(in) :: self
    class(shape), intent(in) :: shapetype
    type(vector) :: center
    type(string) :: outext
    select type (shapetype)
    type is (shape)
    class is (box)
        center = shapetype%pt + Utils%m2geo(shapetype%size, shapetype%pt%y)/2.0
    class is (point)
        center = shapetype%pt
    class is (line)
        center = shapetype%pt + (shapetype%last-shapetype%pt)/2.0
    class is (polyline)
        center = polylineCenter(shapetype)
    class is (sphere)
        center = shapetype%pt
    class is (polygon)
        center = (shapetype%bbMin + shapetype%bbMax)/2.0
        class default
        outext='[geometry::getCenter] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    end function getCenter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that returns the points defining a given geometry
    !> @param[in] self, shapetype
    !---------------------------------------------------------------------------
    function getPoints(self, shapetype) result(pts)
    class(geometry_class), intent(in) :: self
    class(shape), intent(in) :: shapetype
    type(vector), dimension(:), allocatable :: pts
    type(string) :: outext
    integer :: n
    type(vector) :: temp
    select type (shapetype)
    type is (shape)
    class is (box)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        temp = shapetype%size
        pts(1) = shapetype%pt
        pts(2) = shapetype%pt + temp%y*ey
        pts(3) = pts(2) + temp%z*ez
        pts(4) = shapetype%pt + temp%z*ez
        pts(5) = shapetype%pt + temp%x*ex
        pts(6) = pts(5) + temp%y*ey
        pts(7) = shapetype%pt + temp
        pts(8) = pts(5) + temp%z*ez
    class is (point)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        pts(1) = shapetype%pt
    class is (line)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        pts(1) = shapetype%pt
        pts(2) = shapetype%last
    class is (polyline)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        pts(1) = shapetype%pt
        pts(2:) = shapetype%point
    class is (sphere)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        pts(1) = shapetype%pt
    class is (polygon)
        n = self%getnumPoints(shapetype)
        allocate(pts(n))
        pts = shapetype%vertex
        class default
        outext='[geometry::getPoints] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    end function getPoints

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method the points defining a given geometry
    !> @param[in] self, shapetype
    !---------------------------------------------------------------------------
    integer function getnumPoints(self, shapetype)
    class(geometry_class), intent(in) :: self
    class(shape), intent(in) :: shapetype
    type(string) :: outext
    select type (shapetype)
    type is (shape)
    class is (box)
        getnumPoints = 8
    class is (point)
        getnumPoints = 1
    class is (line)
        getnumPoints = 2
    class is (polyline)
        getnumPoints = size(shapetype%point) + 1
    class is (sphere)
        getnumPoints = 1
    class is (polygon)
        getnumPoints  = size(shapetype%vertex)
        class default
        outext='[geometry::getnumPoints] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    end function getnumPoints

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to set the bounding box of a polygon
    !> @param[in] self, zMin, zMax
    !---------------------------------------------------------------------------
    subroutine setBoundingBox(self, zMin, zMax)
    class(polygon), intent(inout) :: self
    real(prec), intent(in), optional :: zMin
    real(prec), intent(in), optional :: zMax
    integer :: i
    self%bbMin = self%vertex(1)
    self%bbMax = self%vertex(1)
    do i= 2, size(self%vertex)
        self%bbMin%x = min(self%bbMin%x, self%vertex(i)%x)
        self%bbMin%y = min(self%bbMin%y, self%vertex(i)%y)
        self%bbMin%z = min(self%bbMin%z, self%vertex(i)%z)
        self%bbMax%x = max(self%bbMax%x, self%vertex(i)%x)
        self%bbMax%y = max(self%bbMax%y, self%vertex(i)%y)
        self%bbMax%z = max(self%bbMax%z, self%vertex(i)%z)
    end do
    if (present(zMin)) then
        if (present(zMax)) then
            self%bbMin%z = zMin
            self%bbMax%z = zMax
        end if
    end if
    self%pt = Geometry%getCenter(self)
    end subroutine setBoundingBox

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns a polygon in metric units, centered at the origin
    !> @param[in] self
    !---------------------------------------------------------------------------
    type(polygon) function getMetricPolygon(self)
    class(polygon), intent(in) :: self

    allocate(getMetricPolygon%vertex(size(self%vertex)))
    getMetricPolygon%vertex = Utils%geo2m(self%vertex - self%pt)
    call getMetricPolygon%setBoundingBox(self%bbMin%z-self%bbMin%z/2.0, self%bbMax%z-self%bbMin%z/2.0)
    getMetricPolygon%pt = self%pt
    end function getMetricPolygon

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to print the details of a given geometry
    !> @param[in] self, shapetype
    !---------------------------------------------------------------------------
    subroutine printGeometry(self, shapetype)
    class(geometry_class), intent(in) :: self
    class(shape) :: shapetype
    integer :: i
    type(vector) :: temp(2)
    type(string) :: temp_str(6)
    type(string) :: outext

    temp_str(1) = shapetype%pt%x
    temp_str(2) = shapetype%pt%y
    temp_str(3) = shapetype%pt%z
    select type (shapetype)
    type is (shape)
    class is (box)
        temp_str(4) = shapetype%size%x
        temp_str(5) = shapetype%size%y
        temp_str(6) = shapetype%size%z
        outext='      Box at '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')//&
            '       with '//temp_str(4)//' X '//temp_str(5)//' X '//temp_str(6)
    class is (point)
        outext='      Point at '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)
    class is (line)
        temp_str(4) = shapetype%last%x
        temp_str(5) = shapetype%last%y
        temp_str(6) = shapetype%last%z
        outext='      Line from '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')//&
            '       to '//temp_str(4)//' '//temp_str(5)//' '//temp_str(6)
    class is (polyline)
        i = size(shapetype%point)
        temp_str(4) = shapetype%point(i)%x
        temp_str(5) = shapetype%point(i)%y
        temp_str(6) = shapetype%point(i)%z
        outext='      Polyline from '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')//&
            '       to '//temp_str(4)//' '//temp_str(5)//' '//temp_str(6)
    class is (sphere)
        temp_str(4) = shapetype%radius
        outext='      Sphere at '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')//&
            '       with radius '//temp_str(4)
    class is (polygon)
        temp_str(1) = shapetype%bbMin%x
        temp_str(2) = shapetype%bbMin%y
        temp_str(3) = shapetype%bbMin%z
        temp_str(4) = shapetype%bbMax%x
        temp_str(5) = shapetype%bbMax%y
        temp_str(6) = shapetype%bbMax%z
        outext='      polygon contained in min['//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//']'//new_line('a')//&
            '       max['//temp_str(4)//' '//temp_str(5)//' '//temp_str(6)//']'
        class default
        outext='[geometry::print] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select
    call Log%put(outext,.false.)

    end subroutine printGeometry

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Sets a polygon from a given file
    !> @param[in] self, poly, fileName, zMin_str, zMax_str
    !---------------------------------------------------------------------------
    subroutine setPolygon(self, poly, fileName, zMin_str, zMax_str)
    class(geometry_class), intent(in) :: self
    class(polygon), intent(inout) :: poly
    type(string), intent(in) :: fileName, zMin_str, zMax_str
    integer :: PolygonsFile, STAT_CALL, ClientNumber, FromBlock, iflag
    integer :: StartLine, EndLine
    integer :: i
    logical :: found, vertBounds
    character(len = line_length) :: FullBufferLine
    type(string) :: outext!, tag
    real, dimension(:), pointer :: PointCoordinates
    real(prec) :: zMin, zMax
    type(xmlparser_class) :: xmlReader

    vertBounds = .false.
    if (zMin_str /= notRead) then
        if(zMax_str /= notRead) then
            zMin = zMin_str%to_number(kind=1._R8P)
            zMax = zMax_str%to_number(kind=1._R8P)
            vertBounds = .true.
        end if
    end if
    
    outext='-> Reading polygon file '//fileName
    call Log%put(outext)
    if (fileName%extension() == '.xy') then
        PolygonsFile = 0
        call ConstructEnterData(PolygonsFile, fileName%chars(), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) then
            outext='[geometry::setPolygon] : error reading polygon file '//fileName//', stoping'
            call Log%put(outext)
            stop
        end if
        call ExtractBlockFromBuffer(PolygonsFile, ClientNumber, '<beginpolygon>', '<endpolygon>', found, STAT = STAT_CALL)
        if (STAT_CALL .EQ. SUCCESS_) then
            if (found) then
                call GetExtractType(FromBlock = FromBlock)
                call GetBlockSize(PolygonsFile, ClientNumber, StartLine, EndLine, FromBlock, STAT = STAT_CALL)
                allocate(poly%vertex(EndLine - StartLine - 1))
                allocate(PointCoordinates (1:2))
                do i=1, size(poly%vertex)
                    call GetFullBufferLine(PolygonsFile, i+1, FullBufferLine, STAT = STAT_CALL)
                    call GetExtractType(FromBlock = FromBlock)
                    call GetData(PointCoordinates, PolygonsFile, flag = iflag, SearchType = FromBlock, Buffer_Line = i+1, STAT = STAT_CALL)
                    poly%vertex(i)%x = PointCoordinates(1)
                    poly%vertex(i)%y = PointCoordinates(2)
                end do
                deallocate(PointCoordinates)
            else
                call Block_Unlock(PolygonsFile, ClientNumber, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) then
                    outext='[geometry::setPolygon] : error reading block from polygon file '//fileName//', stoping'
                    call Log%put(outext)
                    stop
                end if
            end if
        else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
            if(STAT_CALL .ne. SUCCESS_) then
                outext='[geometry::setPolygon] : error extracting block from polygon file '//fileName//', stoping'
                call Log%put(outext)
                stop
            end if
        end if
        call KillEnterData(PolygonsFile, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) then
            outext='[geometry::setPolygon] : error destroying reader for polygon file '//fileName//', stoping'
            call Log%put(outext)
            stop
        end if
        if (vertBounds) then
            call poly%setBoundingBox(zMin, zMax)
        else
            outext='[geometry::setPolygon] : .xy polygon files demand "verticalBoundingBox" field on setup file, stoping'
            call Log%put(outext)
            stop
        end if
     else if (fileName%extension() == '.kmz') then
         call xmlReader%getPolygonFromKMZFile(fileName,poly%vertex)
         if (vertBounds) then
            call poly%setBoundingBox(zMin, zMax)
        else
            call poly%setBoundingBox()
        end if
    else
        outext='-> Format '//fileName%extension()//' not suported, stoping'
        call Log%put(outext)
        stop
    end if

    end subroutine setPolygon

    !orphan routines from now on

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private function that returns the number of points distributed on a grid
    !> with spacing dp inside a sphere
    !> @param[in] dp, r
    !---------------------------------------------------------------------------
    function sphere_np_count(dp, r) result(np)
    type(vector), intent(in) :: dp
    real(prec), intent(in) :: r
    integer :: np
    integer :: i, j, k, nx, ny, nz
    type(vector) :: pts
    np=0
    nx=int(3*r/dp%x)
    ny=int(3*r/dp%y)
    nz=int(3*r/dp%z)
    do i=1, nx
        do j=1, ny
            do k=1, nz
                pts = (ex*(i-1)*dp%x +ey*(j-1)*dp%y +ez*(k-1)*dp%z) - r*(ex+ey+ez)
                if (pts%normL2() <= r) then
                    np=np+1
                end if
            end do
        end do
    end do
    if (np == 0) then !Just the center point
        np=1
    end if
    end function sphere_np_count
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private function that returns the number of points distributed allong a
    !> line with spacing dp
    !> @param[in] dp, firstPoint, lastPoint
    !---------------------------------------------------------------------------
    integer function lineNpCount(dp, firstPoint, lastPoint)
    type(vector), intent(in) :: dp
    type(vector), intent(in) :: firstPoint
    type(vector), intent(in) :: lastPoint
    type(vector) :: temp
    
    temp = firstPoint - lastPoint
    temp = Utils%geo2m(temp, firstPoint%y)
    lineNpCount = max(int(temp%normL2()/dp%normL2()),1)
    if (lineNpCount == 0) then !Just the center point
        lineNpCount=1
    end if

    end function lineNpCount
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private function that returns the number of points distributed allong a
    !> polyline with spacing dp
    !> @param[in] dp, polylineObj
    !---------------------------------------------------------------------------
    integer function polylineNpCount(dp, polylineObj)
    type(vector), intent(in) :: dp
    type(polyline), intent(in) :: polylineObj    
    integer :: i
    
    polylineNpCount = lineNpCount(dp, polylineObj%pt, polylineObj%point(1))
    do i=1, size(polylineObj%point)-1
        polylineNpCount = polylineNpCount + lineNpCount(dp, polylineObj%point(i), polylineObj%point(i+1))
    end do    
    
    if (polylineNpCount == 0) then !Just the center point
        polylineNpCount=1
    end if

    end function polylineNpCount

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private function that returns the number of points distributed on a grid
    !> with spacing dp inside a polygon
    !> @param[in] dp, polyIn
    !---------------------------------------------------------------------------
    integer function polygonNpCount(dp, polyIn)
    type(vector), intent(in) :: dp
    type(polygon), intent(in) :: polyIn
    type(polygon) :: poly
    integer :: i, j, k, nx, ny, nz
    type(vector) :: pts

    poly = getMetricPolygon(polyIn)
    polygonNpCount = 0
    nx=max(int(abs(poly%bbMax%x - poly%bbMin%x)/dp%x+1),1)
    ny=max(int(abs(poly%bbMax%y - poly%bbMin%y)/dp%y+1),1)
    nz=max(int(abs(poly%bbMax%z - poly%bbMin%z)/dp%z),1)
    do i=1, nx
        do j=1, ny
            do k=1, nz
                pts = poly%bbMin + (ex*(i-1)*dp%x + ey*(j-1)*dp%y + ez*(k-1)*dp%z)
                if (pointInPolygon(pts, poly%vertex, poly%bbMin%z, poly%bbMax%z)) then
                    polygonNpCount = polygonNpCount+1
                end if
            end do
        end do
    end do
    if (polygonNpCount == 0) then !Just the center point
        polygonNpCount=1
    end if

    end function polygonNpCount

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp inside a polygon
    !> @param[in] dp, polyIn, ptlist
    !---------------------------------------------------------------------------
    subroutine polygonGrid(dp, polyIn, ptlist)
    type(vector), intent(in) :: dp
    type(polygon), intent(in) :: polyIn
    type(vector), dimension(:), intent(out) :: ptlist
    type(polygon) :: poly
    integer :: i, j, k, nx, ny, nz, p
    type(vector) :: pts

    poly = getMetricPolygon(polyIn)
    if (size(ptlist) == 1) then
        ptlist(1) = poly%pt
        return
    end if
    nx=max(int(abs(poly%bbMax%x - poly%bbMin%x)/dp%x+1),1)
    ny=max(int(abs(poly%bbMax%y - poly%bbMin%y)/dp%y+1),1)
    nz=max(int(abs(poly%bbMax%z - poly%bbMin%z)/dp%z),1)
    p=0
    do i=1, nx
        do j=1, ny
            do k=1, nz
                pts = poly%bbMin + (ex*(i-1)*dp%x + ey*(j-1)*dp%y + ez*(k-1)*dp%z)                
                if (pointInPolygon(pts, poly%vertex, poly%bbMin%z, poly%bbMax%z)) then
                    p=p+1
                    ptlist(p)=pts
                end if
            end do
        end do
    end do

    end subroutine polygonGrid

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp inside a sphere
    !> @param[in] dp, centerPt, r, np, ptlist)
    !---------------------------------------------------------------------------
    subroutine sphere_grid(dp, centerPt, r, np, ptlist)
    type(vector), intent(in) :: dp
    type(vector), intent(in) :: centerPt
    real(prec), intent(in) :: r
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    integer :: i, j, k, p, nx, ny, nz
    type(vector) :: pts

    if (np == 1) then !Just the center point
        ptlist(1)= centerPt
        return
    end if
    nx=int(3*r/dp%x)
    ny=int(3*r/dp%y)
    nz=int(3*r/dp%z)
    p=0
    do i=1, nx
        do j=1, ny
            do k=1, nz
                pts = (ex*(i-1)*dp%x +ey*(j-1)*dp%y +ez*(k-1)*dp%z) - r*(ex+ey+ez)
                if (pts%normL2() <= r) then
                    p=p+1
                    ptlist(p)=pts
                end if
            end do
        end do
    end do

    end subroutine sphere_grid

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp inside a box
    !> @param[in] dp, size, np, ptlist
    !---------------------------------------------------------------------------
    subroutine box_grid(dp, size, np, ptlist)
    type(vector), intent(in) :: dp
    type(vector), intent(in) :: size
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    integer :: i, j, k, p
    p=0
    do i=1, int(size%x/dp%x)+1
        do j=1, int(size%y/dp%y)+1
            do k=1, int(size%z/dp%z)+1
                p=p+1
                ptlist(p) = (ex*(i-1)*dp%x + ey*(j-1)*dp%y + ez*(k-1)*dp%z)
            end do
        end do
    end do
    if (np == 1) then !Just the origin
        ptlist(1)= 0*ex + 0*ey +0*ez
    end if
    end subroutine box_grid
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the points distributed on a polyline
    !> with spacing dp
    !> @param[in] dp, polylineObj, ptlist
    !---------------------------------------------------------------------------
    subroutine polylineGrid(dp, polylineObj, ptlist)
    type(vector), intent(in) :: dp
    type(polyline), intent(in) :: polylineObj
    type(vector), dimension(:), intent(out) :: ptlist
    type(vector) :: offset
    integer :: i, nptsSegm, accm
    
    if (size(ptlist) == 1) then
        ptlist(1) = polylineObj%pt
        return
    end if
    nptsSegm = lineNpCount(dp, polylineObj%pt, polylineObj%point(1))
    call line_grid(Utils%geo2m(polylineObj%point(1)-polylineObj%pt, polylineObj%pt%y), nptsSegm, ptlist(1:nptsSegm))
    accm = nptsSegm    
    do i=1, size(polylineObj%point)-1
        offset = ptlist(accm)
        nptsSegm = lineNpCount(dp, polylineObj%point(i), polylineObj%point(i+1))
        call line_grid(Utils%geo2m(polylineObj%point(i+1)-polylineObj%point(i), polylineObj%point(i)%y), nptsSegm, ptlist(accm + 1 : accm + nptsSegm), offset)
        accm = accm + nptsSegm
    end do

    end subroutine polylineGrid

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp along a line
    !> @param[in] dist, np, ptlist, offset
    !---------------------------------------------------------------------------
    subroutine line_grid(dist, np, ptlist, offset)
    type(vector), intent(in) :: dist
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    type(vector), intent(in), optional :: offset
    type(vector) :: offset_
    integer :: i, j, k, p
    offset_ = 0*ex + 0*ey +0*ez
    if (present(offset)) offset_ = offset
    do p=1, np
        ptlist(p) = dist*(p-1)/np + offset_
    end do
    if (np == 1) then !Just the origin
        ptlist(1)= 0*ex + 0*ey +0*ez + offset_
    end if
    end subroutine line_grid
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private routine that returns the center of a polyline
    !> @param[in] polylineObj
    !---------------------------------------------------------------------------
    type(vector) function polylineCenter(polylineObj)
    type(polyline), intent(in) :: polylineObj
    integer :: i, num
    
    polylineCenter = polylineObj%pt
    num = 1
    do i=1, size(polylineObj%point)
        polylineCenter = polylineCenter + polylineObj%point(i)
        num = num + 1
    end do
    polylineCenter = polylineCenter/num

    end function polylineCenter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> private function that returns true for a point in a polygon, assumed
    !> regular on the vertical coordinate
    !> @param[in] pt, poly, zmin, zmax
    !---------------------------------------------------------------------------
    logical function pointInPolygon(pt, poly, zmin, zmax)
    type(vector), intent(in) :: pt
    type(vector), dimension(:), intent(in) :: poly
    real(prec), intent(in) :: zmin
    real(prec), intent(in) :: zmax
    integer :: i, ip1
    real(prec) :: t

    pointInPolygon = .false.
    if (pt%z <= zmax ) then
        if (pt%z >= zmin ) then
            do i = 1, size(poly)
                ip1 = mod ( i, size(poly) ) + 1
                if ( poly(ip1)%y < pt%y .eqv. pt%y <= poly(i)%y ) then
                    t = pt%x - poly(i)%x - ( pt%y - poly(i)%y ) * ( poly(ip1)%x - poly(i)%x ) / ( poly(ip1)%y - poly(i)%y )
                    if ( t < 0.0D+00 ) then
                        pointInPolygon = .not.pointInPolygon
                    end if
                end if
            end do
        end if
    end if

    end function pointInPolygon

    end module geometry_mod