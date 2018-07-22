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

    module simulation_vtkwritter

    use common_modules

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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Log file initizalization routine.
    !
    !> @param[in] outpath
    !---------------------------------------------------------------------------
    subroutine initLog(self,outpath)
    implicit none
    class(logger_class), intent(inout) :: self
    type(string), intent(in) :: outpath !< output path were to point the logger
    type(string) :: logfile
    
    logfile = outpath//'MOHIDLagrangianRun.out'
    self%log_unit = 0
    open (unit=self%log_unit,file=logfile%chars(),action="write",status="replace")
    
    end subroutine initLog
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Log serialization routine
    !
    !> @param[in] tologstr,timeoption
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public timestamp builder
    !
    !> @param[in] timestamp
    !---------------------------------------------------------------------------
    subroutine getTimeStamp(timestamp)
    implicit none
    type(string), intent(out) :: timestamp
    character(80) :: temp(8)
    integer :: values(8),i

    call date_and_time(values=values)
    do i=1,8
        write(temp(i),*) values(i)
    enddo
    timestamp=trim(adjustl(temp(1)))//'-'//trim(adjustl(temp(2)))//'-'//trim(adjustl(temp(3)))//' @'//trim(adjustl(temp(5)))//':'//trim(adjustl(temp(6)))//':'//trim(adjustl(temp(7)))
    end subroutine getTimeStamp

  end module simulation_vtkwritter
