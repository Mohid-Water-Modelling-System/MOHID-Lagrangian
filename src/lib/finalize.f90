    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : finalize
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module with the simulation closing related definitions and methods
    !------------------------------------------------------------------------------
    
    module finalize

    use tracer_base
    use simulation_globals
    use source_identity

    use commom_modules

    implicit none
    private

    !Public access procedures
    public :: finalizeMohidLagrangian

    contains
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private closing statement routine.
    !---------------------------------------------------------------------------
    subroutine simulation_end
    implicit none
    type(string) :: outext
    outext='Simulation ended, freeing resources. See you next time'
    call ToLog(outext)
    
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private globals deallocation routine.
    !---------------------------------------------------------------------------
    subroutine deallocate_simulation
    implicit none    
    !deallocating Sources
    deallocate(Source)
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private logfile closing routine.
    !---------------------------------------------------------------------------
    subroutine closelog
    implicit none
    close(Log_unit)    
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public finalization routine. Destroys, deallocates and closes the simulation space
    !---------------------------------------------------------------------------
    subroutine finalizeMohidLagrangian
    implicit none    
    
    call simulation_end
    call deallocate_simulation
    call closelog
    
    end subroutine

    end module finalize