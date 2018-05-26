    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : blocks
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a block class and related methods. A block is a fundamental
    !> type of the model. It contains a sub-domain of the simulation bounding box,
    !> holding all entities inside that sub-domain. It maps to a domain decomposition
    !> parallelization strategy, mainly in distributed memory architectures.
    !------------------------------------------------------------------------------

    module blocks_mod

    use commom_modules
    use simulation_globals_mod
    use tracer_array_mod
    use sources_array_mod

    implicit none
    private

    type block_class
      integer :: id
      type(SourceArray) :: Source_stack
      type(TracerArray) :: Tracer_stack

      type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the block
      type(vector)    ::  Pointmax        !< Point that defines the upper corner of the block
    end type block_class

    !Simulation variables
    type(block_class), allocatable, dimension(:) :: Block

    !Public access vars
    !Public access procedures

    contains



  end module blocks_mod