To build the files you must follow these procedure:

Before anything, and because the cmake chain in hdf5 and netcdf-c is broken beyond salvation, you need to install netCDF4 in your system. In unix you just install it from a repo (apt-get install netcdf or yum install netcdf, or whatever you use), in windows you need to download an installer and make sure everythings gets properly put on your PATH.

- Linux -> you are in luck. Just run ./MakeLibraries.sh in this directory and everything should go smooth. 

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
	-run RunCmake.bat in this directory
	-open the generated projects and compile them individually

