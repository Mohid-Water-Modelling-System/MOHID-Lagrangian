    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_user
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a generic Lagrangian tracer class for material modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods. It holds an arbitrarily sized array for properties and their names,
    !> that a user can then explore in the methods.
    !------------------------------------------------------------------------------

    module tracerUser_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: user_par_class               !<Type - parameters of a Lagrangian tracer object representing a user defined type
        type(string) :: name                                !< Name of the Tracer
        type(string), dimension(:), allocatable :: propName !< Name of the property on the Tracer property array
        integer      :: particulate                         !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec)   :: size                                !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type user_par_class

    type :: user_state_class             !<Type - State variables of a tracer object representing a user defined type
        real(prec) :: radius                            !< Tracer radius (m)
        real(prec), dimension(:), allocatable :: prop   !< Tracer properties array
    end type user_state_class

    type, extends(tracer_class) :: userTracer_class   !<Type - The user-defined Lagrangian tracer class
        type(user_par_class)   :: mpar                  !<To access material parameters
        type(user_state_class) :: mnow                  !<To access material state variables
    contains
    end type userTracer_class

    !Public access vars
    public :: userTracer_class

    !Public access routines
    public :: userTracer

    interface userTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> User-defined Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(userTracer_class) :: constructor
    integer, intent(in) :: id
    class(source_class), intent(in) :: src
    real(prec), intent(in) :: time
    integer, intent(in) :: p

    !use the base class constructor to build the base of our new derived type
    constructor%tracer_class = Tracer(id, src, time, p)
    !VERY NICE IFORT BUG (I think) - only some of the variables get used using the base constructor...
    constructor%par%id = id !forcing
    constructor%par%idsource = src%par%id !forcing
    !now initialize the specific components of this derived type

    end function constructor


    end module tracerUser_mod
