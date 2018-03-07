!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : source_identity
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module that defines a source class and related methods.
!------------------------------------------------------------------------------
    
module source_identity
    
    use tracer    
    use initialize
    use commom_modules
    
    implicit none
    private
    
    type source_par_class               !>Type - parameters of a source object
        integer :: id                       !> unique source identification (integer)
        real(prec) :: emitting_rate         !> Emitting rate of the source (Hz)
        type(string) :: name                !> source name
        type(string) :: property_name       !> source property name
        type(string) :: source_type         !> Source type : 'point', 'line', 'sphere', 'box'
    end type
    
    type source_state_class             !>Type - state variables of a source object
        real(prec_time) :: age              ! time variables
        logical :: active                   !> active switch
        type(vector) :: pos                 !> Position of the source baricenter (m)
        type(vector) :: vel                 !> Velocity of the source (m s-1)        
        real(prec) :: depth                 !> Depth of the source baricenter (m)
        real(prec) :: T                     !> Temperature of the source (Celcius)
    end type 

    type source_stats_class             !>Type - statistical variables of a source object        
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        integer ::particles_emitted         !> Number of emitted particles by this source
        real(prec_wrt) :: acc_T             !> Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !> Number of sampling steps
    end type    

    type source_class                   !>Type - The source class
        type(source_par_class)   :: par     !>To access parameters
        type(source_state_class) :: now     !>To access state variables
        type(source_stats_class) :: stats   !>To access statistics
    end type 
    
    ! General public 
    public :: source_class
    
    
    
    
    
    
end module source_identity 
