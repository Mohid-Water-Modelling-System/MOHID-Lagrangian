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
!> Module that defines an emitter class and related methods. This module is 
!> responsible for building a potential tracer list based on the availble 
!> sources and calling their initializers.
!------------------------------------------------------------------------------
    
module source_emitter
    
    use commom_modules
    use source_identity
    
    implicit none
    private
    
    type :: emitter_t
        integer :: emitted
    contains
    procedure :: initialize
    end type
    
    type(emitter_t) ::  Emitter
    
    !Public access vars
    public :: Emitter
    
    contains
    
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that returns the total number of tracers an input 
    !> source can potentially create
    !
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine initialize(self, srcs, nsrcs)
    implicit none
    class(emitter_t), intent(in) :: self
    class(source_class), intent(inout) :: srcs(nsrcs)
    integer, intent(in) :: nsrcs
    integer :: i
    
    do i=1, nsrcs
        call setotalnp(srcs(i))
        !print*, srcs(i)%stencil%total_np
    end do
    
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the total number of tracers an input 
    !> source will potentially create
    !
    !> @param[in] src
    !---------------------------------------------------------------------------
    subroutine setotalnp(src)
    implicit none
    class(source_class), intent(inout) :: src
    src%stencil%total_np=(src%par%stoptime-src%par%startime)*src%par%emitting_rate*src%stencil%np
    end subroutine
    
    
    
        
    
end module source_emitter 
