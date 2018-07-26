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
    !> Defines a vtk writer class for the exposable to the Simulation 
    !------------------------------------------------------------------------------

    module simulation_vtkwritter_mod

    use common_modules
    use vtk_fortran

    implicit none
    private
    
    type :: vtkwritter_class
        private
        integer :: vtk_unit = 10
    contains
    procedure :: DomainVTK
    end type vtkwritter_class
    
    type(vtkwritter_class) :: vtkWritter

    !Public access vars
    public :: vtkWritter

    !Public access procedures
    !public :: getTimeStamp

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public Tracer writting routine. Writes Tracer data in binary XML VTK
    !> format using an unstructured grid. Serial writer for serial files.
    !---------------------------------------------------------------------------
    subroutine TracerSerialVTK(self, filename)
    implicit none
    class(vtkwritter_class), intent(inout) :: self
    type(string), intent(in) :: filename
    
    integer :: error
    type(vtk_file) :: vtk_file
    
    error = vtk_file%initialize(format='binary', filename=filename%chars()//'.vtu', mesh_topology='UnstructuredGrid')
    
    error = vtk_file%finalize()
    
    end subroutine TracerSerialVTK
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation domain writting routine. Writes binary XML VTK
    !> format using an unstructured grid.
    !---------------------------------------------------------------------------
    subroutine DomainVTK(self, filename)
    implicit none
    class(vtkwritter_class), intent(inout) :: self
    type(string), intent(in) :: filename
    
    type(vtk_file) :: vtk_file
    type(string) :: fullfilename
    type(string) :: outext
    integer :: error
    
    
    
    fullfilename = filename%chars()//'_Domain.vtu'
    outext = '->Writting Domain file '//fullfilename
    call Log%put(outext)
    fullfilename = Globals%Names%outpath//'/'//fullfilename
    
    
    
    error = vtk_file%initialize(format='binary', filename=fullfilename%chars(), mesh_topology='UnstructuredGrid')
    
    error = vtk_file%finalize()
    
    end subroutine DomainVTK
    

  end module simulation_vtkwritter_mod
