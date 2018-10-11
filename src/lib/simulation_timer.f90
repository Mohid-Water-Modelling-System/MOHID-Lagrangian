    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : October 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a timer class using the OMP timer function to get the
    !> wall time. API supports accumulation, several start-stop cycles and printing.
    !------------------------------------------------------------------------------

    module simulation_timer_mod

    use simulation_precision_mod
    use simulation_logger_mod
    use stringifor
    use omp_lib

    implicit none
    private

    !Public access vars
    public :: timer_class

    type :: timer_class
        private !full information hidding - contents are only accessible through the methods
        logical :: initialized = .false.
        type(string) :: name               !< name of the timer
        real(prec_time) :: elapsed = 0.0   !< full elapsed time of all cycles
        real(prec_time) :: last = 0.0      !< elapsed time of the last cycle
        real(prec_time) :: start           !< time of the tic at the current cycle
        real(prec_time) :: stop            !< time of the tac at the current cycle
    contains
    procedure :: initialize => initTimer
    procedure :: getElapsed             !< returns the elasped time on this timer
    procedure :: Tic                    !< starts timming cycle
    procedure :: Toc                    !< ends timming cycle and accumulates the total
    procedure :: print => printElapsed  !< prints elapsed time of this timmer
    procedure :: printLast              !< prints the elapsed time of the last ran cycle
    procedure :: printCurrent           !< prints time of last cycle
    procedure :: printNotInitialized    
    end type timer_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a timer object. Optionaly receives a time to
    !> start the timer at a specified elapsed time
    !> @param[in] this, name, restartime
    !---------------------------------------------------------------------------
    subroutine initTimer(this, name, restartime)
    class(timer_class), intent(inout) :: this
    type(string), intent(in) :: name
    real, optional, intent(in) :: restartime
    this%initialized = .true.
    this%name = name
    this%elapsed = 0.0
    this%last = 0.0
    this%start = 0.0
    this%stop = 0.0
    if (present(restartime)) then
        this%elapsed = restartime
    end if
    end subroutine initTimer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the total elapsed time on this timer
    !---------------------------------------------------------------------------
    real(prec_time) function getElapsed(this)
    class(timer_class), intent(in) :: this
    getElapsed = this%elapsed
    end function getElapsed

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that starts the timmer in the current cycle
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine Tic(this)
    class (timer_class), intent(inout) :: this
    this%start = omp_get_wtime()
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that stops the timmer in the current cycle and stores the timmings
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine Toc(this)
    class (timer_class), intent(inout) :: this
    this%stop = omp_get_wtime()
    this%last = this%stop - this%start
    this%elapsed = this%elapsed + this%last
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the total elapsed time. Rudimentary at best.
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine printElapsed(this)
    class(timer_class), intent(in) :: this
    type(string) :: outext, temp
    if (this%initialized) then
        temp = this%elapsed
        outext = 'Total elapsed time for '// this%name //' is '// temp//' s'
        call Log%put(outext,.false.)
    else
        call this%printNotInitialized()
    end if    
    end subroutine printElapsed

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the elapsed time of the last cycle. Rudimentary at best.
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine printLast(this)
    class(timer_class), intent(in) :: this
    print*, "[timer_class::printLast]: elapsed time from last cycle is ", this%last
    end subroutine printLast

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the current elapsed time. Rudimentary at best.
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine printCurrent(this)
    class(timer_class), intent(in) :: this
    print*, "[timer_class::printCurrent]: current elapsed time is ", omp_get_wtime() - this%start
    end subroutine printCurrent
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the total elapsed time. Rudimentary at best.
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine printNotInitialized(this)
    class(timer_class), intent(in) :: this
    type(string) :: outext
    outext = 'Requested timer not initialized, you might have a bug...'
    call Log%put(outext,.false.)
    end subroutine printNotInitialized

    end module simulation_timer_mod
