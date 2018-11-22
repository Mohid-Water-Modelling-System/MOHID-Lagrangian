To build the files you must follow these procedure:

Before anything, and because the cmake chain in hdf5 and netcdf-c is broken beyond salvation, you need to install netCDF4 in your system. In unix you just install it from a repo (apt-get install netcdf or yum install netcdf, or whatever you use), in windows you need to download an installer and make sure everythings gets properly put on your PATH.

- Linux -> you are in luck. Just run ./MakeLibraries.sh in this folder and everything should go smooth. 

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
	-run RunCmake.bat in the Libs folder
	-if there are any mistakes, you are on your own, maybe try the hardcore mode
	-go to the folder named 'Libs_build_win' in this directory
	-open the project (or run the make if you use nmake) and build the libs by this order
		-penf
		-vecfor
		-befor64
		-stringifor
		-foxy
		-vtkfortran
		-face
		-flap
		-datetime
		-fortran-csv-module
	-now go to the Libs/fox/build directory
	-open the project and compile it
	-copy the outputs in the directories modules/ and lib/ to the same folders in Libs/Libs_build_win
	-Go to the Libs/netcdf-fortran/build directory
	-open the project and compile it

	Hardcore mode
	-make a folder named 'Libs_build_win' in this directory
	-run cmake pointed at that folder, with the source on Libs/ to produce the build project
	-go to the folder, open the project (or run the make if you use nmake) and build the libs by this order
		-penf
		-vecfor
		-befor64
		-stringifor
		-foxy
		-vtkfortran
		-face
		-flap
		-datetime
		-fortran-csv-module
	-now go to the Libs/fox/ directory
	-create a build folder and point cmake at it, using Libs/fox/ as the source
	-open the build folder, open the project and compile it
	-copy the outputs in the directories modules/ and lib/ to the same folders in Libs/Libs_build_win
	-Go to the Libs/netcdf-fortran directory
	-create a build folder and point cmake at it, using Libs/netcdf-fortran/ as the source
	-open the build folder, open the project and compile it

