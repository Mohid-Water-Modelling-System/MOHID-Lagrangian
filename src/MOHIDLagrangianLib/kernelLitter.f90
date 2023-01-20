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
    procedure :: BioFouling
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
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Linear growth of biofouling, which increases density of floating object
    !> @param[in] self, sv
    !---------------------------------------------------------------------------
    function BioFouling(self, sv, dt)
    class(kernelLitter_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Biofouling
    real(prec), intent(in) :: dt
    real(prec), dimension(size(sv%state,1)) :: mass_litter, mass_biofouling, volume_biofouling
    real(prec), dimension(size(sv%state,1)) :: total_mass, total_volume
    real(prec) :: density_biofouling, max_age
    integer :: col_density, col_volume, col_age
    type(string) :: tag
    real(prec):: compute_rate = 86400
    !begin---------------------------------------------------------------------------
    Biofouling = 0.0
    tag = 'density'
    col_density = Utils%find_str(sv%varName, tag, .true.)
    tag = 'volume'
    col_volume = Utils%find_str(sv%varName, tag, .true.)
    tag = 'age'
    col_age = Utils%find_str(sv%varName, tag, .true.)
    
    !Compute only when time of the oldest particle is a multiple of 86400s, to reduce simulation time
    max_age = maxval(sv%state(:,col_age))
    if (mod(max_age,compute_rate) == 0) then    
        !kg/m3
        density_biofouling = 2500.0
        !kg         =        kg/m3            *         m3
        mass_litter = sv%state(:,col_density) * sv%state(:,col_volume)
        
        !This is the main part.
        !If new methods are to be included, they should be added here as new function or routine calls
        where (sv%state(:, col_age) == 0)
            !initial biofouling mass in kg is a function of size (or weight) of litter
            mass_biofouling = 0.01 * mass_litter
        elsewhere
            ! Kg            =        Kg       *               1/s                             *     s
            mass_biofouling = mass_biofouling * Globals%Sources%biofouling_rate(sv%source(:)) * compute_rate
        endwhere
    
        !m3               =       Kg        / Kg/m3
        volume_biofouling = mass_biofouling / density_biofouling
        !increase mass. biofilm includes shellfish
        total_mass = mass_litter + mass_biofouling
        !total volume is the sum of litter volume and biofouling volume
        total_volume = sv%state(:,col_volume) + volume_biofouling
        
        !increase in density
        !Kg/(m3.dt)               = Kg/m3 new                  -     Kg/m3 old           / s
        Biofouling(:,col_density) = ((total_mass/total_volume) - sv%state(:,col_density))/dt
        where(sv%state(:,col_density) > 2000.0) sv%active = .false.
        
    end if

    end function BioFouling

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