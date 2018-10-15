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

    module simulation_parallel_omp_mod

    use penf
    use stringifor
    use omp_lib

    implicit none
    private

    type :: parallel_omp_class
        private !full information hidding - contents are only accessible through the methods
        integer :: numThreads       !< Number of OMP threads to use
    contains
    procedure :: initialize => initOMPManager
    procedure :: getMaxThreads
    procedure :: getThreads
    procedure :: setThreads
    end type parallel_omp_class

    type(parallel_omp_class) :: OMPManager

    !Public access vars
    public :: OMPManager

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes an omp data manager object. Sets the number of 
    !> threads as the maximum given by the system.
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine initOMPManager(this)
    class(parallel_omp_class), intent(inout) :: this
    this%numThreads = this%getMaxThreads()
    end subroutine initOMPManager

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the number of threads to the given value. Checks for
    !> the maximum admissible number of threads
    !> @param[in] this, tvalue
    !---------------------------------------------------------------------------
    subroutine setThreads(this, tvalue)
    class(parallel_omp_class), intent(inout) :: this
    integer, intent(in) :: tvalue
    this%numThreads = min(tvalue, this%getMaxThreads())
    CALL OMP_SET_NUM_THREADS(this%numThreads)
    end subroutine setThreads

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the maximum number of OMP threads as given by the system.
    !> @param[in] this
    !---------------------------------------------------------------------------
    integer function getMaxThreads(this)
    class(parallel_omp_class), intent(in) :: this
    getMaxThreads = OMP_GET_NUM_PROCS()
    end function getMaxThreads
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the set number of OMP threads for the simulation.
    !> @param[in] this
    !---------------------------------------------------------------------------
    integer function getThreads(this)
    class(parallel_omp_class), intent(in) :: this
    getThreads = this%numThreads
    end function getThreads


    end module simulation_parallel_omp_mod