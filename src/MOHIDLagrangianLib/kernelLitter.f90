    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelLitter
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for litter associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelLitter_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod


    type :: kernelLitter_class        !< Litter kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelLitter
    procedure :: DegradationLinear
    !procedure :: Buoyancy
    end type kernelLitter_class

    public :: kernelLitter_class
    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Linear degradation kernel.
    !> @param[in] self, sv
    !---------------------------------------------------------------------------
    function DegradationLinear(self, sv)
    class(kernelLitter_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DegradationLinear
    integer :: nf, idx
    type(string) :: tag

    DegradationLinear = 0.0
    tag = 'condition'
    nf = Utils%find_str(sv%varName, tag, .true.)
    tag = 'degradation_rate'
    idx = Utils%find_str(sv%varName, tag, .true.)

    DegradationLinear(:,nf) = -sv%state(:,idx)
    where(sv%state(:,nf) < 0.0) sv%active = .false.

    end function DegradationLinear

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelLitter(self)
    class(kernelLitter_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernelLitter

    end module kernelLitter_mod