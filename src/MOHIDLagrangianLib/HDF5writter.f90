    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_vtkwritter
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : July 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines a vtk writer class with an object exposable to the Output streamer.
    !> Writes files in .xml vtk, both in serial and parallel model. Uses an
    !> unstructured mesh format specifier to store any type of data, both meshes and
    !> Tracers. Supports scalar and vectorial data.
    !------------------------------------------------------------------------------

    module hdf5Writter_mod

    use common_modules

    use ModuleHDF5
    
    use ModuleTime
    
    use ModuleGlobalData
    
    implicit none
    private
    
    public :: writeTestmatrix
    
    contains
    
    subroutine writeTestmatrix(mat)
    integer, dimension(:, :), pointer, intent(in) :: mat
    integer :: ID, HDF5_CREATE, STAT_CALL
    
    ID = 0
    
        !Gets File Access Code
        call GetHDF5FileAccess(HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5(ID, "filehdf", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'writeTestmatrix - module hdf5Writter_mod - ERR01'
        
        call HDF5WriteData(ID, "/Data", "TestMatrix", "-",         &
                              Array2D = mat,            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'writeTestmatrix - module hdf5Writter_mod - ERR07'
        
    end subroutine writeTestmatrix
    

    end module hdf5Writter_mod
