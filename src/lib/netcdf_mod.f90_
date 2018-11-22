    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : common_modules
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Module read netcdf in a easy way.
    !------------------------------------------------------------------------------

    module netcdf_mod

    use common_modules
    use netcdf

    implicit none
    private
    
    type :: grid
        real, dimension(:), allocatable :: vector
    end type grid

    type :: field
        real, dimension(:,:,:,:), allocatable :: data
    end type field

    type :: nc_header
        type(string) :: file_name
        type(string), dimension(:), allocatable :: dim_name
        type(string), dimension(:), allocatable :: var_name
    end type

    type :: nc_class
        type(nc_header) :: header
        integer :: dim_size, var_size
        integer, dimension(:), allocatable :: dim_id, dim_len, var_id
        type(grid), dimension(:), allocatable :: dim_data
        type(field), dimension(:), allocatable :: var_data
        integer :: nc_id
        integer :: status
    contains
    procedure :: initNc
    procedure :: getNcid
    procedure :: getData
    procedure :: closeNcid
    procedure :: check
    end type nc_class
    
    !Exposed public class
    public :: nc_class
    
    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Get the ID for the Netcdf file readed.
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNcid(self)
    class(nc_class),intent(inout) :: self
    self%status = NF90_open(self%header%file_name%chars(), NF90_NOWRITE, self%nc_id)
    call self%check()
    end subroutine getNcid
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Allocate the grid structure and the ammout of data.
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine initNc(self)
    class(nc_class),intent(inout) :: self
    integer :: var_size, dim_size, i

    self%dim_size = size(self%header%dim_name)
    self%var_size = size(self%header%var_name)
    allocate(self%dim_id(self%dim_size), self%dim_len(self%dim_size))
    allocate(self%dim_data(self%dim_size))
    allocate(self%var_id(self%var_size))
    allocate(self%var_data(self%var_size))


    do i=1,self%dim_size
        self%status = NF90_inq_dimid(self%nc_id, self%header%dim_name(i)%chars(), self%dim_id(i))
        call self%check()
        self%status = nf90_inquire_dimension(self%nc_id, self%dim_id(i), len=self%dim_len(i))
        call self%check()
        self%status = nf90_inq_varid(self%nc_id, self%header%dim_name(i)%chars(), self%dim_id(i))
        call self%check()
        allocate(self%dim_data(i)%vector(self%dim_len(i)))

        print*,'Dimension: ', self%header%dim_name(i)%chars()
        print*,'Id dimension: ', self%dim_id(i)
        print*,'Dimension size:' ,self%dim_len(i)
    end do


    do i=1,self%var_size
        self%status = NF90_inq_varid(self%nc_id, self%header%var_name(i)%chars(), self%var_id(i))
        call self%check()
        allocate(self%var_data(i)%data(self%dim_len(1),self%dim_len(2),self%dim_len(3),self%dim_len(4)))

        print*,'Field: ', self%header%var_name(i)%chars()
        !print*,'Field shape:' ,shape(self%var_data(i)%data)
        print*,'Field id: ', self%var_id(i)
    end do

    end subroutine initNc
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Read the data dimensions anda data variables
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getData(self)
    class(nc_class), intent(inout) :: self
    integer :: i

    do i=1,self%dim_size
        self%status = NF90_GET_VAR(self%nc_id, self%dim_id(i), self%dim_data(i)%vector)
        call self%check()
    end do
    do i=1,self%var_size
        self%status = NF90_GET_VAR(self%nc_id, self%var_id(i), self%var_data(i)%data)
        !print*, shape(self%var_data(i)%data)
        call self%check()
    end do

    end subroutine getData
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> close the netcdf file
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine closeNcid(self)
    class(nc_class),intent(inout) :: self
    self%status = NF90_close(self%nc_id)
    call self%check()
    end subroutine closeNcid


    subroutine check(self)
    class(nc_class), intent(inout) :: self
    if(self%status /= NF90_noerr) then
        print *, trim(NF90_strerror(self%status))
        stop "Stopped"
    end if
    end subroutine check


    end module netcdf_mod
