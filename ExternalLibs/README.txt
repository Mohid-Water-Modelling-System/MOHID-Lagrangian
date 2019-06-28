To build the files you must follow these procedure:

- Linux -> you are in luck. Just run ./MakeLibraries.sh in this directory and everything should go smooth. 

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
	-run RunCmake.bat in this directory
	-open the generated projects and compile them individually
	-heavy libs such as HDF5 and NETCDF are already precompiled as shared, lucky you