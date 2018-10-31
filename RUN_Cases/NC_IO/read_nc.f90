program test


use netcdf_mod

	type(nc_class) :: input
	
	allocate(input%dimensions(4))
	allocate(input%dims(4))
	allocate(input)
	allocate(input%file(1))


	input%file%name = 'CMEMS_201008_sfc.nc'
	input%dimensions(1)%name= 'longitude'
	input%dimensions(2)%name= 'latitude'
	input%dimensions(3)%name= 'time'
	input%dimensions(4)%name= 'depth'

	allocate(input%dimsids(4))
	call input_flow_field%getNcid()
	call input_flow_field%getDims()


end program test
