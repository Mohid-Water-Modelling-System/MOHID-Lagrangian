!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracers
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION:
!> Module to hold and 'wrap' all the tracer respective modules.
!------------------------------------------------------------------------------

module tracers_mod

    use tracer_base_mod     ! The base type - pure Lagrangian Tracer
    use tracer_plastic_mod  ! Plastic particles
    use tracer_paper_mod    ! Paper particles

end module tracers_mod
