    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_output_streamer
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : July 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines a output file writer class with an object exposable to the Simulation
    !> This class is in charge of selectig the correct writter for the selected output
    !> file format.
    !------------------------------------------------------------------------------

    module simulationOutputStreamer_mod

    use common_modules
    use vtkWritter_mod
    use boundingbox_mod
    use blocks_mod

    implicit none
    private

    type :: output_streamer_class       !< Output Streamer class
        real(prec) :: OutputFrequency = MV      !< Output frequency to write simulation outputs
        real(prec) :: LastWriteTime = MV        !< Time stamp of the last output write
        integer :: OutputFormat = -1            !< Switch for output format
        type(vtkwritter_class) :: vtkWritter    !< The vtk writter object
    contains
    procedure :: initialize => initOutputStreamer
    procedure :: WriteDomain
    procedure :: WriteStepSerial
    procedure :: finalize => closeOutputStreamer
    procedure, private :: CheckWriteTime
    end type output_streamer_class

    !Public access vars
    public :: output_streamer_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Streamer method to call a simulation step writer. Writes binary XML VTK
    !> format using an unstructured grid.
    !> @param[in] self, blocks
    !---------------------------------------------------------------------------
    subroutine WriteStepSerial(self, blocks)
    class(output_streamer_class), intent(inout) :: self
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    type(string) :: filename                                !< name of the case to add
    if (self%CheckWriteTime()) then
        filename = Globals%Names%casename//'_'//Utils%int2str('(i5.5)',Globals%Sim%getnumoutfile())
        if (self%OutputFormat == 2) then !VTK file selected
            call self%vtkWritter%TracerSerial(filename, blocks)
        end if
        call Globals%Sim%increment_numoutfile()
        self%LastWriteTime = Globals%SimTime%CurrTime
    end if
    end subroutine WriteStepSerial

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation domain writting routine. Writes binary XML VTK
    !> format using an unstructured grid.
    !> @param[in] self, filename, bbox, npbbox, blocks
    !---------------------------------------------------------------------------
    subroutine WriteDomain(self, filename, bbox, npbbox, blocks)
    class(output_streamer_class), intent(inout) :: self
    type(string), intent(in) :: filename                    !< name of the case to add
    class(boundingbox_class), intent(in) :: bbox            !< Case bounding box
    integer, intent(in) :: npbbox                           !< number of points of the bbox geometry
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    if (self%OutputFormat == 2) then !VTK file selected
        call self%vtkWritter%Domain(filename, bbox, npbbox, blocks)
    end if
    end subroutine WriteDomain

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Streamer method to check if this timestep is appropriate to write an
    !> output file
    !> @param[in] self
    !---------------------------------------------------------------------------
    logical function CheckWriteTime(self)
    class(output_streamer_class), intent(inout) :: self
    CheckWriteTime = .false.
    if ((Globals%SimTime%CurrTime - self%LastWriteTime) >= 1.0/self%OutputFrequency) CheckWriteTime = .true.
    if (Globals%SimTime%CurrTime == 0.0) CheckWriteTime = .true.
    end function CheckWriteTime

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the Output writer object
    !---------------------------------------------------------------------------
    subroutine initOutputStreamer(self)
    class(output_streamer_class), intent(inout) :: self
    self%OutputFormat = Globals%Parameters%OutputFormat
    self%OutputFrequency = Globals%Parameters%TimeOut
    self%LastWriteTime = Globals%SimTime%CurrTime
    if (self%OutputFormat == 2) then !VTK file selected
        call self%vtkWritter%initialize()
    end if
    end subroutine initOutputStreamer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Closes the Output writer object
    !---------------------------------------------------------------------------
    subroutine closeOutputStreamer(self)
    class(output_streamer_class), intent(inout) :: self
    if (self%OutputFormat == 2) then !VTK file selected
        call self%vtkWritter%finalize()
    end if
    end subroutine closeOutputStreamer


    end module simulationOutputStreamer_mod