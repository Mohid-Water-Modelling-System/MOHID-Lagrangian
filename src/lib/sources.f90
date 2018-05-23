!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : source
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION:
!> Module to hold and wrap all the tracer sources respective modules. Defines a source class and related methods.
!------------------------------------------------------------------------------

module sources_mod

    use tracers_mod
    use initialize_mod
    use finalize_mod

    use source_identity_mod
    use source_emitter_mod

end module sources_mod
