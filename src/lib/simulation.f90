    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the simulation class and its methods
    !------------------------------------------------------------------------------

    module simulation_mod

    use commom_modules

    use initialize_mod
    use finalize_mod

    implicit none
    private

    type :: simulation_class   !< Parameters class

    contains
    procedure :: initialize => initSimulation
    end type


    !Simulation variables
    public :: simulation_class
    !Public access vars


    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Globals default setting routine.
    !---------------------------------------------------------------------------
    subroutine initSimulation(self, casefilename, outpath)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string), intent(in) :: casefilename         !< case file name
    type(string), intent(in) :: outpath              !< Output path
    
    ! Initialize logger
    call initMohidLagrangianLog(outpath)
    ! Initialization routines to build the simulation from the input case file
    call initMohidLagrangian(casefilename)
    
    end subroutine

    end module simulation_mod
