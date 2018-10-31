module netcdf_mod

use netcdf

implicit none

type :: string
	character(len=:),allocatable :: name
end type string

type,extends(string) :: nc_class
	type(string) :: file_name
	integer :: nc_id
	contains
end type nc_class

type,extends(nc_class) :: nc_dims_class
	real,dimension(:),allocatable :: dims_data
	real,dimension(:),allocatable :: dims_ids, dims_lens
	integer :: status
	contains

end type nc_dims_class

type,extends(nc_dims) :: nc_var_class
	real :: varid
	real :: var_data

contains
	procedure :: getNcId
	procedure :: getDims
	procedure :: check
	procedure :: closeNcId

end type nc_class

contains


subroutine getNcid(self)
	class(nc_var_class),intent(inout) :: self
	self%status = NF90_open(trim(self%file_name), NF90_NOWRITE, self%nc_id)
	call self%check()
end subroutine getNcid


subroutine getDims(self)
	class(nc_var_class), intent(inout) :: self

	self%status = NF90_inq_dimid(self%nc_id, self%name, self%dims_ids)
	self%status = NF90_inq_dimlen(self%nc_id, self%dimid,self%dims_lens)
	allocate(self%dim%dims_data(dimlens))
	self%status = NF90_get_var(self%ncid, self%dims_id,self%dim_data)

	call self%check()
end subroutine getDims


! subroutine getVarid(self)
! 	class(nc_class),intent(inout) :: self
! 	self%status = NF90_inq_varid(self%ncid, trim(self%varname), self%varid)
! 	self%status = NF90_get_var(self%ncid, self%varid, self%variable)
! 	call self%check()
! end subroutine



subroutine closeNcid(self)
	class(nc_class),intent(inout) :: self
	self%status = NF90_close(self%ncid)
	call self%check()
end subroutine closeNcid


subroutine check(self)
    class(nc_class), intent(inout) :: self
    
    if(self%status /= NF90_noerr) then
      print *, trim(NF90_strerror(self%status))
      stop "Stopped"
    end if
end subroutine check


! subroutine readNc(self)

! 	class(nc_class) :: self

! 	call self%getNcid()
! 	print*,self%ncid,'readed'
! 	call self%getVarid()
! 	print*,self%varid
! 	call self%closeNcid()

! end subroutine readNc

end module netcdf_mod


