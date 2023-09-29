    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_logger
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold all the simulation logger related definitions and methods
    !------------------------------------------------------------------------------

    module simulationLogger_mod

    use penf
    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    
    use datetime_module

    implicit none
    private

    type :: logger_class
        private
        integer :: log_unit = -1
    contains
    procedure :: initialize => initLog
    procedure :: finalize   => closeLog
    procedure :: put        => put_inLog
    end type logger_class

    type(logger_class) :: Log

    !Public access vars
    public :: Log

    !Public access procedures
    public :: getTimeStamp

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Log file initizalization routine.
    !> @param[in] self, outpath, filename
    !---------------------------------------------------------------------------
    subroutine initLog(self,outpath, filename)
    implicit none
    class(logger_class), intent(inout) :: self
    type(string), intent(in) :: outpath !< output path were to point the logger
    type(string), intent(in) :: filename
    type(string) :: logfile
    logfile = outpath//filename//'_run.log'
    self%log_unit = 0

    open (unit=self%log_unit,file=logfile%chars(),action="write",status="replace")

    end subroutine initLog

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Log file closure routine.
    !---------------------------------------------------------------------------
    subroutine closeLog(self)
    implicit none
    class(logger_class), intent(inout) :: self
    close(self%log_unit)
    end subroutine closeLog

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Log serialization routine
    !> @param[in] self,tologstr,timeoption
    !---------------------------------------------------------------------------
    subroutine put_inLog(self,tologstr,timeoption)
    implicit none
    class(logger_class), intent(in) :: self
    type(string), intent(inout) :: tologstr
    logical, intent(in), optional :: timeoption
    type(string) :: timestamp

    call getTimeStamp(timestamp)
    if (present(timeoption)) then
        if (.not.timeoption) then
            timestamp=''
        endif
    endif
    tologstr=timestamp//' '//tologstr
    write(self%log_unit,"(A)") tologstr%chars()
    print'(A)', tologstr%chars()

    end subroutine put_inLog

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public timestamp builder
    !> @param[in] timestamp
    !---------------------------------------------------------------------------
    subroutine getTimeStamp(timestamp)
    type(string), intent(out) :: timestamp
    type(datetime) :: date   
    date = date%now()
    timestamp = date%isoformat(' ')
    timestamp = timestamp%basename(strip_last_extension=.true.) !getting rid of miliseconds
    end subroutine getTimeStamp


    end module simulationLogger_mod
