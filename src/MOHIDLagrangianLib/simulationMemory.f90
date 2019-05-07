    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_memory
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the simulation memory managment class and its methods
    !------------------------------------------------------------------------------

    module simulationMemory_mod

    use stringifor
    use simulationLogger_mod
    use simulationPrecision_mod

    implicit none
    private

    type :: memory_t       !< Case memory occupation logger
        private
        integer :: size_of_sources   !< Size of the sources in memory (bytes)
        integer :: size_of_tracers   !< Size of the tracers in memory (bytes)
        integer :: size_of_defs      !< Size of the parameters and definitions in memory (bytes)
        integer :: size_of_blocks    !< Size of the Blocks in memory (bytes)
        integer :: ntrc              !< Expected number of Tracers for the simulation (by Source emission at least)
        integer :: sizeTrc           !< Size of a dummy Tracer, in bytes
    contains
    procedure :: initialize => initializeMemory
    procedure :: addblock
    procedure :: addsource
    procedure :: setracer
    procedure :: adddef
    procedure :: getotal
    procedure :: setNtrc
    procedure :: setsizeTrc
    procedure :: print => printmemory
    procedure :: detailedprint => printmemorydetailed
    end type memory_t

    !Simulation variables
    type(memory_t) :: SimMemory

    !Public access vars
    public :: SimMemory

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Memory logger initialization method.
    !---------------------------------------------------------------------------
    subroutine initializeMemory(self)
    implicit none
    class(memory_t), intent(inout) :: self
    self%size_of_sources = 0
    self%size_of_tracers = 0
    self%size_of_defs = 0
    self%size_of_blocks = 0
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to retreive the total size of the allocated memory.
    !---------------------------------------------------------------------------
    subroutine getotal(self,size)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(out) :: size
    size = self%size_of_sources + self%size_of_tracers + self%size_of_defs + self%size_of_blocks
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to set the total expected number of Tracers.
    !---------------------------------------------------------------------------
    subroutine setNtrc(self,ntrc)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(in) :: ntrc
    self%ntrc = ntrc
    end subroutine setNtrc

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to set the size of a typical Tracer.
    !---------------------------------------------------------------------------
    subroutine setsizeTrc(self,sizeTrc)
    implicit none
    class(memory_t), intent(inout) :: self
    integer*8, intent(in) :: sizeTrc
    self%sizeTrc = sizeTrc
    end subroutine setsizeTrc

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to add the size of a Block to the memory log.
    !---------------------------------------------------------------------------
    subroutine addblock(self,size)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(in) :: size
    self%size_of_blocks = self%size_of_blocks + size
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to add the size of a Source to the memory log.
    !---------------------------------------------------------------------------
    subroutine addsource(self,size)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(in) :: size
    self%size_of_sources = self%size_of_sources + size
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to add the size of a Tracer to the memory log.
    !---------------------------------------------------------------------------
    subroutine setracer(self,size)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(in) :: size
    self%size_of_tracers = size
    end subroutine setracer


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to add the size of a definition to the memory log.
    !---------------------------------------------------------------------------
    subroutine adddef(self,size)
    implicit none
    class(memory_t), intent(inout) :: self
    integer, intent(in) :: size
    self%size_of_defs = self%size_of_defs + size
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the total allocated memory.
    !---------------------------------------------------------------------------
    subroutine printmemory(self)
    implicit none
    class(memory_t), intent(inout) :: self
    integer :: size
    real(prec) :: sizemb
    type(string) :: outext,temp
    call self%getotal(size)
    sizemb = size*1E-6
    temp= sizemb
    outext='->Total allocated memory: '//temp//' mb'
    call Log%put(outext)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the allocated memory.
    !---------------------------------------------------------------------------
    subroutine printmemorydetailed(self)
    implicit none
    class(memory_t), intent(inout) :: self
    integer :: size
    real(prec) :: sizemb
    type(string) :: outext,temp(6)
    call self%getotal(size)
    sizemb = size*1E-6
    temp(6) = 2.25*sizemb
    temp(1)= sizemb
    sizemb = self%size_of_sources*1E-6
    temp(2)= sizemb
    sizemb = self%size_of_tracers*1E-6
    temp(3)= sizemb
    sizemb = self%size_of_defs*1E-6
    temp(4)= sizemb
    sizemb = self%size_of_blocks*1E-6
    temp(5)= sizemb
    sizemb = self%ntrc*self%sizeTrc*1E-6
    outext='->Total allocated memory: '//temp(1)//' mb'//new_line('a')//&
        '       Allocated memory for Blocks  = '//temp(5)//' mb'//new_line('a')//&
        '       Allocated memory for Sources = '//temp(2)//' mb'//new_line('a')//&
        '       Allocated memory for Tracers = '//temp(3)//' mb'//new_line('a')//&
        '       Allocated memory for Consts  = '//temp(4)//' mb'//new_line('a')//&
        '       Expected memory requirements exceed '//temp(6)//' mb'
    call Log%put(outext)
    end subroutine

    end module simulationMemory_mod
