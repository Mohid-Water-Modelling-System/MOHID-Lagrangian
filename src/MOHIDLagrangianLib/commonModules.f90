    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : common_modules
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold all of the commonly used base modules.
    !------------------------------------------------------------------------------

    module common_modules

    use penf
    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    use datetime_module

    use omp_lib

    use geometry_mod
    use simulationPrecision_mod
    use simulationLogger_mod
    use simulationTimer_mod
    use simulationMemory_mod
    use simulationGlobals_mod
    use simulationParallel_omp_mod
    use utilities_mod

    end module common_modules
