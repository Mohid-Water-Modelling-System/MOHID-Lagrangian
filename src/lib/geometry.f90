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
    use stringifor
    use simulation_precision    
    use simulation_logger

    implicit none
    private

    !Update with any new additions as they are added
    type(string), allocatable, dimension(:) :: GeomList !< String list (array) with the name of possible geometry types.
        
    type, abstract :: shape                      !<Type - extendable shape class        
        type(vector) :: pt              !< Coordinates of a point
    contains
    procedure :: getnp                  !<Gets the number of points that define that geometry (based on GLOBALS::dp)
    procedure :: getpointdistribution   !<Gets the actual list of points always referant to the origin (based on GLOBALS::dp)
    end type
    
    type, extends(shape) :: point   !<Type - point class
    end type

    type, extends(shape) :: line    !<Type - line class        
        type(vector) :: last            !< Coordinates of the end point
    end type

    type, extends(shape) :: sphere  !<Type - sphere class        
        real(prec) :: radius            !< Sphere radius
    end type
    
    type, extends(shape) :: box     !<Type - point class        
        type(vector) :: size            !< Box size
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
        np = max((int(self%size%x/dp)+1)*(int(self%size%y/dp)+1)*(int(self%size%z/dp)+1),1)
    class is (point)
        np = 1
    class is (line)
        temp = self%pt-self%last !simply the difference
        np = max(int(temp%normL2()/dp),1)
    class is (sphere)
        call sphere_np_count(dp, self%radius, np)
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
        call box_grid(dp, self%size, np, ptlist)
    class is (point)
        ptlist(1)=self%pt
    class is (line)
        call line_grid(dp, self%last-self%pt, np, ptlist)
    class is (sphere)
        call sphere_grid(dp, self%radius, np, ptlist)
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
    !> private routine that returns the number of points distributed on a grid
    !> with spacing dp inside a sphere
    !
    !> @param[in] dp, r, np
    !---------------------------------------------------------------------------
    subroutine sphere_np_count(dp, r, np)
    implicit none
    real(prec), intent(in) :: dp
    real(prec), intent(in) :: r
    integer, intent(out) :: np
    integer :: i, j, k, n
    type(vector) :: pts
    n=int(3*r/dp)
    do i=1, n
        do j=1, n
            do k=1, n
                pts = dp*(ex*(i-1)+ey*(j-1)+ez*(k-1)) - r*(ex+ey+ez)
                if (pts%normL2() .le. r) then
                    np=np+1
                end if
            end do
        end do
    end do
    if (np == 0) then !Just the center point
        np=1
    end if
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp inside a sphere
    !
    !> @param[in] dp, r, np, ptlist
    !---------------------------------------------------------------------------
    subroutine sphere_grid(dp, r, np, ptlist)
    implicit none
    real(prec), intent(in) :: dp
    real(prec), intent(in) :: r
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    integer :: i, j, k, p, n
    type(vector) :: pts
    n=int(3*r/dp)
    p=0
    do i=1, n
        do j=1, n
            do k=1, n
                pts = dp*(ex*(i-1)+ey*(j-1)+ez*(k-1)) - r*(ex+ey+ez)
                if (pts%normL2() .le. r) then
                    p=p+1
                    ptlist(p)=pts
                end if
            end do
        end do
    end do
    if (np == 1) then !Just the center point
        ptlist(1)= 0*ex + 0*ey +0*ez
    end if
    return
    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp inside a box
    !
    !> @param[in] dp, size, np, ptlist
    !---------------------------------------------------------------------------
    subroutine box_grid(dp, size, np, ptlist)
    implicit none
    real(prec), intent(in) :: dp
    type(vector), intent(in) :: size
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    integer :: i, j, k, p
    p=0
    do i=1, int(size%x/dp)+1
        do j=1, int(size%y/dp)+1
            do k=1, int(size%z/dp)+1
                p=p+1
                ptlist(p) = dp*(ex*(i-1)+ey*(j-1)+ez*(k-1))
            end do
        end do
    end do
    if (np == 1) then !Just the origin
        ptlist(1)= 0*ex + 0*ey +0*ez
    end if
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the points distributed on a grid
    !> with spacing dp along a line
    !
    !> @param[in] dp, size, np, ptlist
    !---------------------------------------------------------------------------
    subroutine line_grid(dp, dist, np, ptlist)
    implicit none
    real(prec), intent(in) :: dp
    type(vector), intent(in) :: dist
    integer, intent(in)::  np
    type(vector), intent(out) :: ptlist(np)
    integer :: i, j, k, p

    do p=1, np
        ptlist(p) = dp/np*(dist*(p-1))
    end do
    if (np == 1) then !Just the origin
        ptlist(1)= 0*ex + 0*ey +0*ez
    end if
    return
    end subroutine

    end module geometry
