    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : blocks
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a block class and related methods. A block is a fundamental
    !> type of the model. It maps to a sub-domain of the simulation bounding box,
    !> holding all entities inside that sub-domain. It maps to a domain decomposition
    !> parallelization strategy, mainly in distributed memory architectures.
    !------------------------------------------------------------------------------

    module blocks

    use commom_modules
    use simulation_globals

    use pointer_linked_list
    !use container_linked_list

    implicit none
    private

    type block_class
      integer :: id
      !type(container_list) :: Source_stack
      !type(container_list) :: Tracer_stack
      type(pointer_list)   :: Tracer_list

      real(prec) :: xmin = 0._R8P  !> x limit min
      real(prec) :: xmax = 0._R8P  !> x limit max
      real(prec) :: ymin = 0._R8P  !> y limit min
      real(prec) :: ymax = 0._R8P  !> y limit max
    end type block_class

    !Public access vars
    public :: block_class
    !Public access procedures

    contains



  end module blocks
