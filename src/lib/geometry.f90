    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : source
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines geometry classes and related methods.
    !------------------------------------------------------------------------------

    module geometry

    use vecfor
    use simulation_precision
    use stringifor
    use simulation_logger

    implicit none
    private

    !Update with any new additions as they are added
    type(string), allocatable, dimension(:) :: GeomList !< String list (array) with the name of possible geometry types.

    !>Type - extendable shape class
    type shape
        !> Coordinates of a point
        type(vector) :: pt
    contains
    procedure :: getnp
    procedure :: getpointdistribution
    end type

    !>Type - point class
    type, extends(shape) :: point
    end type

    !>Type - line class
    type, extends(shape) :: line
        !> Coordinates of the end point
        type(vector) :: last
    end type

    !>Type - sphere class
    type, extends(shape) :: sphere
        !> Sphere radius
        real(prec) :: radius
    end type

    !>Type - point class
    type, extends(shape) :: box
        ! Coordinates of the lower left corner point are defined by shape class
        !> Box size
        type(vector) :: size
    end type

    public :: shape, point, line, sphere, box
    public :: GeomList
    public :: AllocateGeomList, IsValidGeom

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to get the number of points that fill a given geometry
    !
    !> @param[in] self, np
    !---------------------------------------------------------------------------
    subroutine getnp(self,np,dp)
    implicit none
    class(shape) :: self
    integer, intent(out) :: np
    real(prec), intent(in) :: dp
    type(vector) :: temp
    type(string) :: outext

    select type (self)
    type is (shape)
    class is (box)
        np = max(int(self%size%x/dp)*int(self%size%y/dp)*int(self%size%z/dp),1)
    class is (point)
        np = 1
    class is (line)
        temp = self%pt-self%last !simply the difference
        np = max(int(temp%normL2()/dp),1)
    class is (sphere)
        call sphere_np_count(dp, self%radius, self%pt, np)
        class default
        outext='geometry::getnp : unexpected type for geometry object!'
        call ToLog(outext)
        stop
    end select

    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to get the number of points that fill a given geometry
    !
    !> @param[in] self, np
    !---------------------------------------------------------------------------
    subroutine getpointdistribution(self,np,dp,ptlist)
    implicit none
    class(shape) :: self
    integer, intent(in) :: np
    real(prec), intent(in) :: dp
    type(vector), intent(inout) :: ptlist(np)
    type(vector) :: temp
    type(string) :: outext

    select type (self)
    type is (shape)
    class is (box)
        
    class is (point)
        
    class is (line)
        
    class is (sphere)
        call sphere_grid(dp, self%radius, self%pt, np, ptlist)
        class default
        outext='geometry::getpointdistribution : unexpected type for geometry object!'
        call ToLog(outext)
        stop
    end select

    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public routine to allocate the possible geometry name list
    !---------------------------------------------------------------------------
    subroutine AllocateGeomList
    implicit none
    allocate(GeomList(4))
    GeomList(1) ='point'
    GeomList(2) ='line'
    GeomList(3) ='box'
    GeomList(4) ='sphere'
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public function that returns a logical if the input geometry name is valid
    !
    !> @param[in] geomname
    !---------------------------------------------------------------------------
    logical function IsValidGeom(geomname) result(tf)
    implicit none
    type(string), intent(in) :: geomname
    integer :: i
    tf = .false.
    do i=1, size(GeomList)
        if (geomname == GeomList(i)) then
            tf = .true.
        endif
    enddo
    end function IsValidGeom

    !---------------------------------------------------------------------------
    !> @John Burkardt
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the number of points distributed on a grid
    !> inside a sphere, centered on its center.
    !    The grid is defined by specifying the radius and center of the sphere,
    !    and the number of subintervals N into which the horizontal radius
    !    should be divided.  Thus, a value of N = 2 will result in 5 points
    !    along that horizontal line.
    !
    !> @param[in] dp, r, c, ng
    !Input, real dp, point spacing
    !Input, real r, the radius of the sphere.
    !Input, real c(3), the coordinates of the center of the sphere.
    !Output, integer ng, the number of grid points inside the sphere.
    !---------------------------------------------------------------------------
    subroutine sphere_np_count(dp, r, center, ng)
    implicit none
    real(prec), intent(in) :: dp
    real(prec), intent(in) :: r
    type(vector), intent(in) :: center
    integer, intent(out) :: ng
    integer :: i, j, k, n
    real(prec) :: x, y, z, c(3)

    c(1)=center%x
    c(2)=center%y
    c(3)=center%z
    n = int(r/dp)
    !print*, "N = ", n
    ng = 0
    do i = 0, n
        x = c(1) + r * real ( 2 * i, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
        do j = 0, n
            y = c(2) + r * real ( 2 * j, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
            do k = 0, n
                z = c(3) + r * real ( 2 * k, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
                if ( r * r < ( x - c(1) )**2 &
                    + ( y - c(2) )**2 &
                    + ( z - c(3) )**2 ) then
                exit
                end if
                ng = ng + 1
                if ( 0 < i ) then
                    ng = ng + 1
                end if
                if ( 0 < j ) then
                    ng = ng + 1
                end if
                if ( 0 < k ) then
                    ng = ng + 1
                end if
                if ( 0 < i .and. 0 < j ) then
                    ng = ng + 1
                end if
                if ( 0 < i .and. 0 < k ) then
                    ng = ng + 1
                end if
                if ( 0 < j .and. 0 < k ) then
                    ng = ng + 1
                end if
                if ( 0 < i .and. 0 < j .and. 0 < k ) then
                    ng = ng + 1
                end if
            end do
        end do
    end do

    return
    end

    
    !---------------------------------------------------------------------------
    !> @John Burkardt
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns points distributed on a grid
    !> inside a sphere, centered on its center.
    !    The grid is defined by specifying the radius and center of the sphere,
    !    and the number of subintervals N into which the horizontal radius
    !    should be divided.  Thus, a value of N = 2 will result in 5 points
    !    along that horizontal line.
    !
    !> @param[in] dp, r, center, ng, ptlist
    !Input, real dp, point spacing
    !Input, real r, the radius of the sphere.
    !Input, real c(3), the coordinates of the center of the sphere.
    !Input, integer ng, the number of grid points inside the sphere.
    !---------------------------------------------------------------------------
    subroutine sphere_grid(dp, r, center, ng, ptlist)   
    implicit none
    real(prec), intent(in) :: dp
    real(prec), intent(in) :: r
    type(vector), intent(in) :: center
    integer, intent(in)::  ng
    type(vector), intent(out) :: ptlist(ng)
    real(prec) :: bg(3,ng)
    integer :: i, j, k, p, n
    real(prec) :: x, y, z, c(3)

    c(1)=center%x
    c(2)=center%y
    c(3)=center%z
    n = int(r/dp)
    p = 0
    do i = 0, n
        x = c(1) + r * real ( 2 * i, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
        do j = 0, n
            y = c(2) + r * real ( 2 * j, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
            do k = 0, n
                z = c(3) + r * real ( 2 * k, kind = 8 ) / real ( 2 * n + 1, kind = 8 )
                if ( r * r < ( x - c(1) )**2 &
                    + ( y - c(2) )**2 &
                    + ( z - c(3) )**2 ) then
                exit
                end if
                p = p + 1
                bg(1,p) = x
                bg(2,p) = y
                bg(3,p) = z
                if ( 0 < i ) then
                    p = p + 1
                    bg(1,p) = 2.0D+00 * c(1) - x
                    bg(2,p) = y
                    bg(3,p) = z
                end if
                if ( 0 < j ) then
                    p = p + 1
                    bg(1,p) = x
                    bg(2,p) = 2.0D+00 * c(2) - y
                    bg(3,p) = z
                end if
                if ( 0 < k ) then
                    p = p + 1
                    bg(1,p) = x
                    bg(2,p) = y
                    bg(3,p) = 2.0D+00 * c(3) - z
                end if
                if ( 0 < i .and. 0 < j ) then
                    p = p + 1
                    bg(1,p) = 2.0D+00 * c(1) - x
                    bg(2,p) = 2.0D+00 * c(2) - y
                    bg(3,p) = z
                end if
                if ( 0 < i .and. 0 < k ) then
                    p = p + 1
                    bg(1,p) = 2.0D+00 * c(1) - x
                    bg(2,p) = y
                    bg(3,p) = 2.0D+00 * c(3) - z
                end if
                if ( 0 < j .and. 0 < k ) then
                    p = p + 1
                    bg(1,p) = x
                    bg(2,p) = 2.0D+00 * c(2) - y
                    bg(3,p) = 2.0D+00 * c(3) - z
                end if
                if ( 0 < i .and. 0 < j .and. 0 < k ) then
                    p = p + 1
                    bg(1,p) = 2.0D+00 * c(1) - x
                    bg(2,p) = 2.0D+00 * c(2) - y
                    bg(3,p) = 2.0D+00 * c(3) - z
                end if
            end do
        end do
    end do
    
    do p=1, ng
        ptlist(p)%x=bg(1,p)
        ptlist(p)%y=bg(2,p)
        ptlist(p)%z=bg(3,p)
    end do

    return
    end

    end module geometry
