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

    module geometry_mod

    use vecfor
    use stringifor
    use simulation_precision_mod
    use simulation_logger_mod
    use simulation_globals_mod

    implicit none
    private

    !Update with any new additions as they are added
    type :: geometry_class
      type(string), allocatable, dimension(:) :: list !< String list (array) with the name of possible geometry types.
    contains
      procedure :: initialize => allocatelist
      procedure :: inlist
      procedure :: fillsize !<Gets the number of points that fill a geometry (based on GLOBALS::dp)
      procedure :: fill !<Gets the list of points that fill a geometry (based on GLOBALS::dp)
    end type

    type :: shape                      !<Type - extendable shape class
        type(vector) :: pt              !< Coordinates of a point
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

    type(geometry_class) :: Geometry

    public :: shape, point, line, sphere, box
    public :: Geometry

    contains


      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Public routine to allocate the possible geometry name list
      !---------------------------------------------------------------------------
      subroutine allocatelist(self)
      implicit none
      class(geometry_class), intent(inout) :: self
      allocate(self%list(4))
      self%list(1) ='point'
      self%list(2) ='line'
      self%list(3) ='box'
      self%list(4) ='sphere'
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
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> method to get the number of points that fill a given geometry
      !
      !> @param[in] shapetype
      !---------------------------------------------------------------------------
      function fillsize(self, shapetype)
      implicit none
      class(geometry_class), intent(in) :: self
      class(shape), intent(in) :: shapetype
      real(prec) :: dp
      integer :: fillsize
      type(vector) :: temp
      type(string) :: outext

      dp = Globals%SimDefs%Dp
      select type (shapetype)
      type is (shape)
      class is (box)
          fillsize = max((int(shapetype%size%x/dp)+1)*(int(shapetype%size%y/dp)+1)*(int(shapetype%size%z/dp)+1),1)
      class is (point)
          fillsize = 1
      class is (line)
          temp = shapetype%pt-shapetype%last
          fillsize = max(int(temp%normL2()/dp),1)
      class is (sphere)
          fillsize = sphere_np_count(dp, shapetype%radius)
          class default
          outext='[geometry::np] : unexpected type for geometry object, stoping'
          call Log%put(outext)
          stop
      end select

    end function fillsize

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to get the list of points that fill a given geometry
    !
    !> @param[in] shapetype,fillsize,ptlist
    !---------------------------------------------------------------------------
    subroutine fill(self,shapetype,fillsize,ptlist)
    implicit none
    class(geometry_class), intent(in) :: self
    class(shape) :: shapetype
    integer, intent(in) :: fillsize
    type(vector), intent(out) :: ptlist(fillsize)
    type(vector) :: temp
    type(string) :: outext

    select type (shapetype)
    type is (shape)
    class is (box)
        call box_grid(Globals%SimDefs%Dp, shapetype%size, fillsize, ptlist)
    class is (point)
        ptlist(1)=shapetype%pt
    class is (line)
        call line_grid(Globals%SimDefs%Dp, shapetype%last-shapetype%pt, fillsize, ptlist)
    class is (sphere)
        call sphere_grid(Globals%SimDefs%Dp, shapetype%radius, fillsize, ptlist)
        class default
        outext='[geometry::fill] : unexpected type for geometry object, stoping'
        call Log%put(outext)
        stop
    end select

  end subroutine fill


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private function that returns the number of points distributed on a grid
    !> with spacing dp inside a sphere
    !
    !> @param[in] dp, r
    !---------------------------------------------------------------------------
    function sphere_np_count(dp, r) result(np)
    implicit none
    real(prec), intent(in) :: dp
    real(prec), intent(in) :: r
    integer :: np
    integer :: i, j, k, n
    type(vector) :: pts
    np=0
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
    end function

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

  end module geometry_mod
