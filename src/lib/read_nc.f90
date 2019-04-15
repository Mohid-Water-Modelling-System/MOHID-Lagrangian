program test


use netcdf_mod

type(nc_class) :: input

allocate(input%header%dim_name(4))
allocate(input%header%var_name(2))


input%header%file_name%name = 'MOHID_Vigo_20180904_0000.nc4'
input%header%dim_name(1)%name = 'lon'
input%header%dim_name(2)%name = 'lat'
input%header%dim_name(3)%name = 'depth'
input%header%dim_name(4)%name = 'time'

input%header%var_name(1)%name = 'u'
input%header%var_name(2)%name = 'v'



!call input%getDims()
!call input%getVars()
call input%getncid()
call input%initNc()
call input%getData()


end program test