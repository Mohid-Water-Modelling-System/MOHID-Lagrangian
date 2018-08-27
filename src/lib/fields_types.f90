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

    module field_types_mod

    use common_modules

    implicit none
    private
    
    type :: scalar_field_class      !< a scalar field class
        type(string) :: name        !< name of the field
        type(string) :: units       !< units of the field, preferably SI please
        real(prec), allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    !procedure :: initialize => initScalarField
    end type scalar_field_class

    !Public access vars
    public :: scalar_field_class

    contains
    
    
  end module field_types_mod
