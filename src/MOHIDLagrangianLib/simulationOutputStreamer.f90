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
    use hdf5Writter_mod
    use boundingbox_mod
    use blocks_mod

    implicit none
    private

    type :: output_streamer_class       !< Output Streamer class
        real(prec) :: OutputIntervalTime = MV      !< Output interval to write simulation outputs
        real(prec) :: LastWriteTime = MV        !< Time stamp of the last output write
        integer :: OutputFormat = -1            !< Switch for output format
        type(string), allocatable, dimension(:) :: outputVariables  !< list of optional variables to output
        type(vtkwritter_class) :: vtkWritter    !< The vtk writter object
        type(hdf5writter_class) :: hdf5Writter    !< The vtk writter object
    contains
    procedure :: initialize => initOutputStreamer
    procedure :: writeOutputHeader
    procedure, private :: writeOutputSummary
    procedure :: WriteDomain
    procedure :: WriteStep
    procedure, private :: WriteStepSerial
    procedure :: finalize => closeOutputStreamer
    procedure, private :: CheckWriteTime
    end type output_streamer_class

    !Public access vars
    public :: output_streamer_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> output streamer method to check if it is writ time, and call a 
    !> step writer. Assembles output file name and updates streamer data.
    !> @param[in] self, blocks, numTracers, simTimer
    !---------------------------------------------------------------------------
    subroutine WriteStep(self, blocks, numTracers, simTimer)
    class(output_streamer_class), intent(inout) :: self
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    integer, intent(in) :: numTracers
    type(timer_class), intent(in) :: simTimer
    type(string) :: fileName
    
    if (self%CheckWriteTime()) then
        fileName = Globals%Names%casename//'_'//Utils%int2str('(i5.5)',Globals%Output%getnumOutFile())
        call self%WriteStepSerial(fileName, numTracers, blocks, self%outputVariables)
        call self%writeOutputSummary(numTracers, simTimer, fileName)
        call Globals%Output%setlastOutNumDt(Globals%Sim%getnumdt())
        call Globals%Output%increment_numOutFile()
        self%LastWriteTime = Globals%SimTime%CurrTime
    end if

    end subroutine WriteStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Streamer method to call an appropriate writer.
    !> @param[in] self, filename, numTracers, blocks, outputVars
    !---------------------------------------------------------------------------
    subroutine WriteStepSerial(self, filename, numTracers, blocks, outputVars)
    class(output_streamer_class), intent(inout) :: self
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    integer, intent(in) :: numTracers
    type(string), intent(in) :: filename                    !< name of the case to add
    type(string), dimension(:), intent(in) :: outputVars    !< names of the output variables to print
    

    if (self%OutputFormat == 2) then !VTK file selected
        call self%vtkWritter%TracerSerial(filename, numTracers, blocks, outputVars)
        !call self%hdf5Writter%TracerSerial(filename, blocks)
    end if

    end subroutine WriteStepSerial

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation domain writting routine.
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
    if ((Globals%SimTime%CurrTime - self%LastWriteTime) >= self%OutputIntervalTime) CheckWriteTime = .true.
    if (Globals%SimTime%CurrTime == 0.0) CheckWriteTime = .true.
    end function CheckWriteTime
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Writes simulation log header
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine writeOutputHeader(self)
    class(output_streamer_class), intent(in) :: self
    type(string) :: outext
    outext =         '==================================================================================================================='//new_line('a')
    outext = outext//'                                             Simulation starting'//new_line('a')
    outext = outext//' ==================================================================================================================='
    call Log%put(outext,.false.)
    outext = '    Output time    |   Simulation time   |      Finish time    |     % | Tracer # |  Steps | sim/time | output file'
    call Log%put(outext, .false.)
    end subroutine writeOutputHeader
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> writes log entry with data regarding current output
    !> @param[in] self, numTracers, simTimer, fileName
    !---------------------------------------------------------------------------
    subroutine writeOutputSummary(self, numTracers, simTimer, fileName)
    class(output_streamer_class), intent(in) :: self
    integer, intent(in) :: numTracers
    type(timer_class), intent(in) :: simTimer
    type(string), intent(in) :: fileName
    type(string) :: outext, temp
    type(datetime) :: finishDateTime
    type(timedelta) :: estimTimeDelta
    integer :: totalSecsToFinish
    temp = Globals%SimTime%CurrDate%isoformat(' ')
    temp = temp%basename(strip_last_extension=.true.)
    outext = '| '//temp
    finishDateTime = finishDateTime%now()
    estimTimeDelta = Globals%SimTime%EndDate - Globals%SimTime%StartDate
    if (simTimer%getElapsedLast() == 0.0) then
        totalSecsToFinish = 0.0
    else
        totalSecsToFinish = estimTimeDelta%total_seconds()/(Globals%SimDefs%dt/simTimer%getElapsedLast())
    endif
    estimTimeDelta = timedelta(0,0,0,totalSecsToFinish,0)
    finishDateTime = finishDateTime + estimTimeDelta
    temp = finishDateTime%isoformat(' ')
    temp = temp%basename(strip_last_extension=.true.)
    outext = outext//' | '//temp
    outext = outext//' | '//Utils%real2str('(f5.1)', (Globals%SimTime%CurrTime/Globals%Parameters%TimeMax*100.0))
    outext = outext//' | '//Utils%int2str('(i8.1)', numTracers)
    outext = outext//' | '//Utils%int2str('(i6.1)', Globals%Sim%getnumdt())
    if (simTimer%getElapsedLast() == 0.0) then
        outext = outext//' | '// Utils%real2str('(f8.2)',Globals%SimDefs%dt*0.0)
    else
        outext = outext//' | '//Utils%real2str('(f8.2)',(Globals%SimDefs%dt/simTimer%getElapsedLast()))
    endif
    outext = outext//' | '//fileName
    call Log%put(outext)
    end subroutine writeOutputSummary

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the Output writer object
    !---------------------------------------------------------------------------
    subroutine initOutputStreamer(self)
    class(output_streamer_class), intent(inout) :: self
    self%OutputFormat = Globals%Parameters%OutputFormat
    self%OutputIntervalTime = Globals%Parameters%OutputWriteTime
    self%LastWriteTime = Globals%SimTime%CurrTime
    call Globals%Output%getOutputPoolArray(self%outputVariables)
    if (self%OutputFormat == 2) then !VTK file selected
        call self%vtkWritter%initialize()
        !call self%hdf5Writter%initialize()
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