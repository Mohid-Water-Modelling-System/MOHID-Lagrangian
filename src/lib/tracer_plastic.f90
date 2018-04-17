    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_plastic
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for plastic modelling and related methods. 
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracer_plastic

    use tracer_base
    use commom_modules

    implicit none
    private

    type :: plastic_par_class               !<Type - parameters of a Lagrangian tracer object representing a plastic material
        real(prec) :: density                       !< density of the material        
        real(prec) :: degradation_rate              !< degradation rate of the material
        logical    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type
    
    type :: plastic_state_class             !<Type - State variables of a tracer object representing a plastic material
        real(prec) :: radius                        !< Tracer radius (m)
        real(prec) :: condition                     !< Material condition (1-0)
        real(prec) :: Concentration                 !< Particle concentration
    end type

    type, extends(tracer_class) :: plastic_class    !<Type - The plastic material Lagrangian tracer class
        type(plastic_par_class)   :: mpar     !<To access material parameters
        type(plastic_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: initialize => plastic_initialize
    end type

    !Public access vars
    public :: plastic_class

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Tracer initialization method
    !
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine plastic_initialize(trc,id,id_source,time,pt)
    implicit none
    class(plastic_class) :: trc
    integer, intent(in) :: id
    integer, intent(in) :: id_source
    type(vector), intent(in) :: pt
    real(prec_time), intent(in) :: time

    ! initialize parameters
    trc%par%id = id
    trc%par%idsource = id_source
    trc%par%velmax = 15.0 !(m/s, just a placeholder)
    ! interp_method - TODO
    ! initialize tracer state
    trc%now%age=0.0
    trc%now%active = .false.
    trc%now%pos = pt
    trc%now%vel = 0.0
    trc%now%acc = 0.0
    trc%now%depth = 0.0
    ! Initialize statistical accumulator variables
    trc%stats%acc_pos = 0.0
    trc%stats%acc_vel = 0.0
    trc%stats%acc_depth = 0.0
    trc%stats%ns = 0

    return
    end subroutine

    end module tracer_plastic


