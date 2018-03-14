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
    
    module simulation_logger
    
    use penf
    use vecfor
    use stringifor
    
    implicit none
    private
    
    !Public access vars
    public :: Log_unit
    
    !File handling
    integer, parameter :: Log_unit=0    !> 'Number' of log file
    
     !Public access procedures
    public :: initMohidLagrangianLog, ToLog, getTimeStamp
    
    contains
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
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
    
    end subroutine
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public log serialization routine
    !
    !> @param[in] outpath
    !---------------------------------------------------------------------------
    subroutine ToLog(tologstr)
    implicit none
    type(string), intent(in) :: tologstr
    type(string) :: timestamp
    call getTimeStamp(timestamp)
    timestamp=timestamp//' '//tologstr
    write(Log_unit,"(A)") timestamp%chars()
    print'(A)', timestamp%chars()
    end subroutine
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public log file initizalization routine.
    !
    !> @param[in] outpath
    !---------------------------------------------------------------------------
    subroutine initMohidLagrangianLog(outpath)
    implicit none
    type(string), intent(in) :: outpath    
    type(string) :: logfile        
    logfile = outpath//'MOHIDLagrangianRun.out'
    open (unit=Log_unit,file=logfile%chars(),action="write",status="replace")    
    call getTimeStamp(logfile)        
    end subroutine
    
    end module simulation_logger
