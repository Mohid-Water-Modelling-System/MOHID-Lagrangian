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
    
module sources 
    
    use tracers    
    use initialize
    use finalize

    use source_identity
    use source_emitter
    
end module sources 
