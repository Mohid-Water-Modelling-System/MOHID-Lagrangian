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
    !> 
    !------------------------------------------------------------------------------

    module simulation_input_streamer_mod

    use common_modules

    implicit none
    private

    type :: input_streamer_class
        integer :: InputFormat = -1
    contains
    procedure :: initialize => initOutputStreamer
    procedure :: WriteDomain
    procedure :: WriteStepSerial
    end type input_streamer_class

    type(input_streamer_class) :: InputStreamer

    !Public access vars
    public :: InputStreamer

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the Output writer object
    !---------------------------------------------------------------------------
    subroutine initOutputStreamer(self)
    class(output_streamer_class), intent(inout) :: self
    
    
    end subroutine initOutputStreamer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Streamer method to call a simulation step writer. Writes binary XML VTK
    !> format using an unstructured grid.
    !> @param[in] self, blocks
    !---------------------------------------------------------------------------
    subroutine WriteStepSerial(self)
    class(output_streamer_class), intent(inout) :: self

    end subroutine WriteStepSerial

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation domain writting routine. Writes binary XML VTK
    !> format using an unstructured grid.
    !> @param[in] self, filename, bbox, npbbox, blocks
    !---------------------------------------------------------------------------
    subroutine WriteDomain(self)
    class(output_streamer_class), intent(inout) :: self

    end subroutine WriteDomain

    end module simulation_input_streamer_mod
