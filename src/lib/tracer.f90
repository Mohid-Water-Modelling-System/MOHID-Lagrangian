!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracer
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module to hold and wrap all the tracer respective modules. Defines a pure Lagrangian tracer class. 
!> This is intended to serve as the base class for every type of tracer class needed, that should be     
!> built as derived of this class, with the necessary modifiers to model the desired behaviour.
!> Basic tracer data (parameters, variables) are implemented.
!> Tracer methods such as I/O, integration and interpolation routines are implemented.
!------------------------------------------------------------------------------
    
module tracer 
    
    use simulation_precision
    use tracer_interp    
    use tracer_base

end module tracer 
