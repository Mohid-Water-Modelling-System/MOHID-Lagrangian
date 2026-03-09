To build the MOHID-Lagrangian. Please check the installation guide in UserGuides for more details

- Linux:
    -> Compile libraries: /MOHID-Lagrangian/ExternalLibs/MakeLibraries.sh
    -> Compile MOHID-Lagrangian: ./MakeMOHIDLagrangian.sh.
    -> Create MOHID-Lagrangian python environment: ./MakeEnvMOHIDLagrangian.sh ('conda' must be in path)

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
	-> Run RunCmakeMOHIDLagrangian.bat in this directory (requires Cmake in the path)
	-> Open the generated projects and compile them individually (requires VisualStudio with Ifort)
	-> Heavy libs such as HDF5 and NETCDF are already precompiled as shared, lucky you.
    -> Create MOHID-Lagrangian python environment: ./MakeEnvMOHIDLagrangian.bat ('conda' must be in path)
    
