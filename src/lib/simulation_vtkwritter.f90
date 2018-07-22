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
        integer :: vtk_unit = 2
    contains
    !procedure :: initialize => initLog
    !procedure :: finalize   => closeLog
    !procedure :: put        => put_inLog
    end type vtkwritter_class
    
    type(vtkwritter_class) :: vtkWritter

    !Public access vars
    public :: vtkWritter

    !Public access procedures
    !public :: getTimeStamp

    contains

    

  end module simulation_vtkwritter_mod
