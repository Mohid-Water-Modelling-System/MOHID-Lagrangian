    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_input_streamer
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an inpup file reader class with an object exposable to the Simulation
    !> This class is in charge of selectig the correct reader for the selected input
    !> file format.
    !------------------------------------------------------------------------------

    module simulation_input_streamer_mod

    use common_modules

    implicit none
    private

    type :: input_streamer_class
        integer :: InputFormat = -1
    contains
    procedure :: initialize => initInputStreamer
    
    end type input_streamer_class

    type(input_streamer_class) :: InputStreamer

    !Public access vars
    public :: InputStreamer

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the Input writer object
    !---------------------------------------------------------------------------
    subroutine initInputStreamer(self)
    class(output_streamer_class), intent(inout) :: self
    self%InputFormat = Globals%Parameters%InputFormat
    if (self%InputFormat == 1) then
        !initialize netcdf reader class
    
    end if    
    end subroutine initInputStreamer

    

    end module simulation_input_streamer_mod
