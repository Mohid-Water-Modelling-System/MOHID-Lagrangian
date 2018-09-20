To build the files you must follow these procedure:

- Linux -> you are in luck. Just run ./MakeLibraries.sh in this folder and everything should go smooth. 

- Windows -> Not so lucky. Automated building in windows is tricky, so we need to go by parts
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
	-now go to the Libs/fox/ directory
	-create a build folder and point cmake at it, using Libs/fox/ as the source
	-open the build folder, open the project and compile it
	-copy the outputs in the directories modules/ and lib/ to the same folders in Libs/Libs_build_win

