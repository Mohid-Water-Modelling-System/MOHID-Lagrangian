    program test
    
    use netcdf_mod

    type(nc_class) :: input

    allocate(input%header%dim_name(4))
    allocate(input%header%var_name(2))
    
    input%header%file_name = 'MOHID_Vigo_20180904_0000.nc4'
    input%header%dim_name(1) = 'lon'
    input%header%dim_name(2) = 'lat'
    input%header%dim_name(3) = 'depth'
    input%header%dim_name(4) = 'time'
    
    input%header%var_name(1) = 'u'
    input%header%var_name(2) = 'v'
    
    call input%getNcid()
    !call input%initNc()
    !call input%getData()


    end program test