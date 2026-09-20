To build the files you must follow these procedure:

- Linux -> you are in luck. Just run ./MOHIDLagrangian_install.sh in this directory and everything should go smooth. 
	- open the script and change the configuration parameters if necessary

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
	-run RunCmake.bat in this directory
	-open the generated projects and compile them individually
	-heavy libs such as HDF5 and NETCDF are already precompiled as shared, lucky you