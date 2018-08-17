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
    use vecfor_r4p
    use stringifor
    use datetime_module

    use geometry_mod
    use simulation_precision_mod
    use simulation_logger_mod
    use simulation_memory_mod
    use simulation_globals_mod
    use utilities_mod

end module common_modules
