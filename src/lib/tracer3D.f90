!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracer3D
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module that defines a pure Lagrangian tracer class and related methods.
!
! REVISION HISTORY:
! 26 02 2018 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
    
module tracer3D 

    use tracer_precision
    use tracer_interp   

    implicit none 

    type tracer_par_trans_class !>Type - transient parameters of a pure Lagrangian tracer object        
        
    end type 

    type tracer_par_class !>Type - parameters of a pure Lagrangian tracer object
        integer :: n, n_active, n_max_dep, id_max 
        logical :: is_sigma                     ! Is the defined z-axis in sigma coords
        real(prec_time) :: dt, dt_dep, dt_write 
        real(prec) :: thk_min                   ! Minimum thickness of tracer (m)
        real(prec) :: H_min                     ! Minimum ice thickness to track (m)
        real(prec) :: depth_max                 ! Maximum depth of tracer (fraction)
        real(prec) :: U_max                     ! Maximum horizontal velocity of tracer to track (m/a)
        real(prec) :: U_max_dep                 ! Maximum horizontal velocity allowed for tracer deposition (m/a)
        real(prec) :: H_min_dep                 ! Minimum ice thickness for tracer deposition (m)
        real(prec) :: alpha                     ! Slope of probability function
        character(len=56) :: weight             ! Weighting function for generating prob. distribution
        logical    :: noise                     ! Add noise to gridded deposition location
        real(prec) :: dens_z_lim                ! Distance from surface to count density
        integer    :: dens_max                  ! Max allowed density of particles at surface
        character(len=56) :: interp_method  

        ! Transient parameters 
        character(len=512) :: par_trans_file 
        logical            :: use_par_trans
        type(tracer_par_trans_class) :: tpar 

    end type 

    type tracer_state_class !>Type - state variables of a pure Lagrangian tracer object
        real(prec_time) :: time, time_old
        real(prec_time) :: time_dep, time_write 
        real(prec_time) :: dt  
        integer, allocatable :: active(:), id(:)
        real(prec), allocatable :: x(:), y(:), z(:), sigma(:)
        real(prec), allocatable :: ux(:), uy(:), uz(:)
        real(prec), allocatable :: ax(:), ay(:), az(:)
        real(prec), allocatable :: dpth(:), z_srf(:)
        real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
        real(prec), allocatable :: H(:)

    end type 

    type tracer_stats_class !>Type - statistical variables of a pure Lagrangian tracer object
        
        ! All stats variables at writing precision (prec_wrt), default is single precision
        real(prec_wrt), allocatable :: x(:), y(:)
        real(prec_wrt), allocatable :: depth_norm(:)
        real(prec_wrt), allocatable :: age_iso(:) 

        real(prec_wrt), allocatable :: depth_iso(:,:,:)
        real(prec_wrt), allocatable :: depth_iso_err(:,:,:)
        real(prec_wrt), allocatable :: dep_z_iso(:,:,:)
        integer,    allocatable :: density_iso(:,:,:)
        
        real(prec_wrt), allocatable :: ice_age(:,:,:)
        real(prec_wrt), allocatable :: ice_age_err(:,:,:)
        integer,    allocatable :: density(:,:,:)
           
    end type

    type tracer_dep_class 
        ! Standard deposition information (time and place)
        real(prec), allocatable :: time(:) 
        real(prec), allocatable :: H(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: lon(:), lat(:) 

        ! Additional tracer deposition information (climate, isotopes, etc)
        real(prec), allocatable :: t2m_ann(:), t2m_sum(:)      
        real(prec), allocatable :: pr_ann(:), pr_sum(:)     
        real(prec), allocatable :: t2m_prann(:) ! Precip-weighted temp
        real(prec), allocatable :: d18O_ann(:)
!         real(prec), allocatable :: dD(:)

    end type 

    type tracer_class !>Type - a pure Lagrangian tracer object
        type(tracer_par_class)   :: par 
        type(tracer_state_class) :: now 
        type(tracer_dep_class)   :: dep
        type(tracer_stats_class) :: stats 

    end type 

    private 

    ! For other tracer modules
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_dep_class
    public :: tracer_stats_class

    ! General public 
    public :: tracer_class 
    public :: tracer_init

contains 
    
   !---------------------------------------------------------------------------  
   !> @Ricardo Birjukovs Canelas - MARETEC
   !> Routine Author Name and Affiliation.
   !
   ! DESCRIPTION: 
   !> Brief description of routine. 
   !> @brief
   !> Tracer inititialization routine - Generates a tracer collection and initializes their variables
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[out] trc      
   !> @param[in] filename
   !---------------------------------------------------------------------------
    subroutine tracer_init(trc,filename,time,x,y,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:), y(:)
        logical,    intent(IN) :: is_sigma  
        real(prec_time), intent(IN) :: time 

        ! Local variables 
        integer :: i         

        return 

    end subroutine tracer_init   

end module tracer3D 


