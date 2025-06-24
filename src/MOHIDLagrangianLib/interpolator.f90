    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : interpolator
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an Interpolator class.
    !------------------------------------------------------------------------------

    module interpolator_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use fieldTypes_mod
    USE ieee_arithmetic

    implicit none
    private
    
    type :: interpolator_class        !< Interpolator class
        integer :: interpType = 1     !< Interpolation Algorithm 1: linear
        type(string) :: name                !< Name of the Interpolation algorithm
    contains
    procedure :: run
    procedure, private :: getArrayCoord
    procedure, private :: getArrayCoordRegular
    procedure, private :: getPointCoordRegular
    procedure, private :: getPointCoordNonRegular
    procedure, private :: getArrayCoordNonRegular
    procedure, private :: getArrayCoord_curv
    procedure, private :: Trc2Grid_4D
    procedure, private :: isPointInsidePolygon
    procedure :: initialize => initInterpolator
    procedure :: print => printInterpolator
    procedure, private :: test4D
    procedure, private :: interp4D
    procedure, private :: interp4D_Hor
    procedure, private :: interp3D
    end type interpolator_class

    !Public access vars
    public :: interpolator_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> V2 in Dec 2022 by Joao Sobrinho - Colab Atlantic
    !> Method that runs the chosen interpolator method on the given data.
    !> @param[in] self, state, bdata, time, var_dt, var_name, toInterp
    !---------------------------------------------------------------------------
    subroutine run(self, state, bdata, time, var_dt, var_name, toInterp, reqVertInt)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:,:), intent(in) :: state
    type(background_class), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(:,:), intent(out) :: var_dt
    type(string), dimension(:), intent(out) :: var_name
    type(string), dimension(:), intent(in), optional :: toInterp
    logical, intent(in), optional :: reqVertInt
    real(prec), dimension(:,:,:,:), pointer :: bathymetry
    real(prec), dimension(:), allocatable :: var_dt_aux
    logical :: interp, requireVertInt
    real(prec) :: newtime
    class(*), pointer :: aField
    integer :: i, count
    logical :: gridIsCurvilinear
    type(string) :: outext
    real(prec), dimension(size(state,1)) :: xx, yy, zz
    logical, dimension(size(state,1)) :: outOfBounds
    real(prec) :: tt
    !begin----------------------------------------------------------------------------
    !Check field extents and what particles will be interpolated
    !interpolate each field to the correspoing slice in var_dt
    !write(*,*)"Entrei interpolador"
    i = 1
    call bdata%fields%reset()                   ! reset list iterator
    do while(bdata%fields%moreValues())         ! loop while there are values
        interp = .true.
        requireVertInt = .true.
        aField => bdata%fields%currentValue()   ! get current value
        select type(aField)
        class is(scalar4d_field_class)          !4D interpolation is possible
            !write(*,*)"Variavel que entro no while = ", trim(aField%name)
            if (self%interpType == 1) then !linear interpolation in space and time
                if (present(toInterp)) then
                    !write(*,*)"toInterp present"
                    !do count=1,size(toInterp)
                    !    write(*,*)"variavel para interpolar = ", trim(toInterp(count))
                    !    write(*,*)"variavel no field = ", trim(aField%name)
                    !enddo
                    if (.not.(any(toInterp == aField%name))) then
                        !write(*,*)"toInterp false"
                        interp = .false.
                    end if
                end if
                if (interp) then
                    outOfBounds = .false.
                    var_name(i) = aField%name
                    gridIsCurvilinear = bdata%getGridIsCurvilinear()
                    !write(*,*)"Sai do gridIsCurvilinear", gridIsCurvilinear
                    
                    !produce xx, yy and zz vectors and tt value according to grid type
                    call self%Trc2Grid_4D(state,bdata,time,xx,yy,zz,tt,outOfBounds,gridIsCurvilinear)
                    !write(*,*)"Sai do Trc2Grid_4D"
                    
                    if (var_name(i) == Globals%Var%landIntMask) then
                        !adjust interpolation to bathymetry rather than vertical layers. important for ressuspension processes
                        !write(*,*)"get bathymetry"
                        call bdata%getVarByName4D(varName = Globals%Var%bathymetry, outField_4D = bathymetry, origVar = aField%name)
                        allocate(var_dt_aux(size(state,1)))
                        !write(*,*)"Entrar em interp4D_Hor"
                        var_dt_aux = self%interp4D_Hor(xx, yy, zz, tt, outOfBounds, bathymetry, size(bathymetry,1), size(bathymetry,2), size(bathymetry,3), size(bathymetry,4), size(state,1))
                        !write(*,*)"Saida de interp4D_Hor"
                        zz = self%getArrayCoord(state(:,3), bdata, Globals%Var%level, outOfBounds, bat = var_dt_aux)
                        !write(*,*)"Saida de getArrayCoord"
                        deallocate(var_dt_aux)
                        nullify(bathymetry)
                    end if
                    !write(*,*)"entrada para if de reqVertInt"
                    if (present(reqVertInt)) then
                        !Interpolate on 2D even if the field is 3D (usefull for Bottom stress)
                        requireVertInt = reqVertInt
                    end if
                    if (requireVertInt) then
                        !write(*,*)"entrada requireVertInt...........Variavel = ", trim(aField%name)
                        var_dt(:,i) = self%interp4D(xx, yy, zz, tt, outOfBounds, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(aField%field,4), size(state,1))
                        !write(*,*)"Saida requireVertInt. var_dt(1,i) = ", var_dt(1,i)
                    else
                        !write(*,*)"Entrada requireVertInt false"
                        var_dt(:,i) = self%interp4D_Hor(xx, yy, zz, tt, outOfBounds, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(aField%field,4), size(state,1))
                        !write(*,*)"Saia requireVertInt false. var_dt(1,i) = ", var_dt(1,i)
                    end if
                end if
            end if !add more interpolation types here
        class is(scalar3d_field_class)          !3D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                if (present(toInterp)) then
                    if (.not.(any(toInterp == aField%name))) then
                        interp = .false.
                    end if
                end if
                if (interp) then
                    var_name(i) = aField%name
                    outOfBounds = .false.
                    xx = self%getArrayCoord(state(:,1), bdata, Globals%Var%lon, outOfBounds)
                    yy = self%getArrayCoord(state(:,2), bdata, Globals%Var%lat, outOfBounds)
                    tt = self%getPointCoordNonRegular(time, bdata, Globals%Var%time)
                    var_dt(:,i) = self%interp3D(xx, yy, tt, outOfBounds, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(state,1))
                end if
            end if !add more interpolation types here
            !add more field types here
        class default
            outext = '[Interpolator::Run] Unexepected type of field, not correct or supported at the time'
            call Log%put(outext)
            stop
        end select
        call bdata%fields%next()                ! increment the list iterator
        i = i+1 !to select the correct slice of var_dt for the corresponding field
    end do
    call bdata%fields%reset()                   ! reset list iterator
    !write(*,*)"Sai interpolador"
    end subroutine run

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 4d interpolation is a weighted average of 16
    !> neighbors. Consider the 4D domain between the 16 neighbors. The hypercube is
    !> divided into 16 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube.
    !> @param[in] self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp4D(self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y, z                       !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                           !< time to interpolate to in array coordinates
    logical, dimension(:), intent(in) :: out
    real(prec), dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: field    !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv                         !< field dimensions
    integer, intent(in) :: n_e                                            !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    real(prec), dimension(n_e) :: xd, yd, zd, c000, c100, c010, c110, c001
    real(prec), dimension(n_e) :: c101, c011, c111, c00, c10, c01, c11, c0, c1
    real(prec) :: td
    integer :: i, t0, t1
    real(prec), dimension(n_e) :: interp4D                                !< Field evaluated at x,y,z,t
    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
    !write(*,*)"Entrada interp4D"
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            x0(i) = floor(x(i))
            x1(i) = ceiling(x(i))
            y0(i) = floor(y(i))
            y1(i) = ceiling(y(i))
            z0(i) = floor(z(i))
            z1(i) = ceiling(z(i))
        endif
    end do
    
    t0 = floor(t)
    t1 = ceiling(t)

    ! If depth layer has one layer
    if (n_pv == 1) then
        z0 = 1
        z1 = 1
    end if
    xd = 0.
    yd = 0.
    zd = 0.
    td = 0.
    ! Compute the "normalized coordinates" of the particle inside the data field box
    where (x1 /= x0) xd = (x-x0)/(x1-x0)
    where (y1 /= y0) yd = (y-y0)/(y1-y0)
    where (z1 /= z0) zd = (z-z0)/(z1-z0)
    if (t1 /= t0) td = (t-t0)/(t1-t0)
    ! Interpolation on the first dimension and collapse it to a three dimension problem
    interp4D = 0.0
    !write(*,*)"tamanho field : ", size(field,1), size(field,2), size(field,3), size(field,4)
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            !write(*,*)"y0(i),x0(i),z0(i),t0 = ", y0(i),x0(i),z0(i),t0
            !write(*,*)"y1(i),x1(i),z1(i),t1 = ", y1(i),x1(i),z1(i),t0
            !write(*,*)"xd(i),yd(i),zd(i),td = ", xd(i),yd(i),zd(i),td
            !write(*,*)"field(y0(i),x0(i),z0(i),t0) = ", field(y0(i),x0(i),z0(i),t0)
            !write(*,*)"field(y0(i),x1(i),z0(i),t0) = ", field(y0(i),x1(i),z0(i),t0)
            !write(*,*)"field(y0(i),x0(i),z1(i),t0) = ", field(y0(i),x0(i),z1(i),t0)
            !write(*,*)"field(y0(i),x1(i),z1(i),t0) = ", field(y0(i),x1(i),z1(i),t0)
        
            c000(i) = field(y0(i),x0(i),z0(i),t0)*(1.-xd(i)) + field(y0(i),x1(i),z0(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
            !write(*,*)"c000 = ", c000(i)
            c100(i) = field(y1(i),x0(i),z0(i),t0)*(1.-xd(i)) + field(y1(i),x1(i),z0(i),t0)*xd(i)
            !write(*,*)"c100 = ", c000(i)
            c010(i) = field(y0(i),x0(i),z1(i),t0)*(1.-xd(i)) + field(y0(i),x1(i),z1(i),t0)*xd(i)
            !write(*,*)"c010 = ", c000(i)
            c110(i) = field(y1(i),x0(i),z1(i),t0)*(1.-xd(i)) + field(y1(i),x1(i),z1(i),t0)*xd(i)
            !write(*,*)"c110 = ", c000(i)
        
            c001(i) = field(y0(i),x0(i),z0(i),t1)*(1.-xd(i)) + field(y0(i),x1(i),z0(i),t1)*xd(i) !y0x0z0t0!  y0x1z0t0
            c101(i) = field(y1(i),x0(i),z0(i),t1)*(1.-xd(i)) + field(y1(i),x1(i),z0(i),t1)*xd(i)
            c011(i) = field(y0(i),x0(i),z1(i),t1)*(1.-xd(i)) + field(y0(i),x1(i),z1(i),t1)*xd(i)
            c111(i) = field(y0(i),x1(i),z1(i),t1)*(1.-xd(i)) + field(y1(i),x1(i),z1(i),t1)*xd(i)
            !write(*,*)"c001 = ", c001(i)
            !write(*,*)"c101 = ", c101(i)
            !write(*,*)"c011 = ", c011(i)
            !write(*,*)"c111 = ", c111(i)
        ! Interpolation on the second dimension and collapse it to a two dimension problem
            c00(i) = c000(i)*(1.-yd(i))+c100(i)*yd(i)
            c10(i) = c010(i)*(1.-yd(i))+c110(i)*yd(i)
            c01(i) = c001(i)*(1.-yd(i))+c101(i)*yd(i)
            c11(i) = c011(i)*(1.-yd(i))+c111(i)*yd(i)
            !write(*,*)"c00 = ", c00(i)
            !write(*,*)"c10 = ", c10(i)
            !write(*,*)"c01 = ", c01(i)
            !write(*,*)"c11 = ", c11(i)
        ! Interpolation on the third dimension and collapse it to a one dimension problem
            c0(i) = c00(i)*(1.-zd(i))+c10(i)*zd(i)
            !write(*,*)"c0 = ", c0(i)
            c1(i) = c01(i)*(1.-zd(i))+c11(i)*zd(i)
            !write(*,*)"c1 = ", c1(i)
        ! Interpolation on the time dimension and get the final result.
            interp4D(i) = c0(i)*(1.-td)+c1(i)*td
            !write(*,*)"interp4D(i) = ", interp4D(i)
        endif
    end do
    end function interp4D
    
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 3d interpolation is a weighted average of 8
    !> neighbors. Consider the 3D domain between the 8 neighbors. The hypercube is
    !> divided into 8 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube. vertical interpolation is excluded
    !> @param[in] self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp4D_Hor(self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y, z                       !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                           !< time to interpolate to in array coordinates
    logical, dimension(:), intent(in) :: out
    real(prec), dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: field    !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv                         !< field dimensions
    integer, intent(in) :: n_e                                            !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    real(prec), dimension(n_e) :: xd, yd
    real(prec), dimension(n_e) :: c00, c10, c01, c11, c0, c1
    real(prec) :: td
    real(prec) :: bottom_lower_left_t0, bottom_lower_right_t0
    real(prec) :: bottom_upper_left_t0, bottom_upper_right_t0
    real(prec) :: ceiling_lower_left_t0, ceiling_lower_right_t0
    real(prec) :: ceiling_upper_left_t0, ceiling_upper_right_t0
    real(prec) :: bottom_lower_left_t1, bottom_lower_right_t1
    real(prec) :: bottom_upper_left_t1, bottom_upper_right_t1
    real(prec) :: ceiling_lower_left_t1, ceiling_lower_right_t1
    real(prec) :: ceiling_upper_left_t1, ceiling_upper_right_t1
    integer :: i, t0, t1
    real(prec), dimension(n_e) :: interp4D_Hor                                !< Field evaluated at x,y,z,t
    !Begin---------------------------------------------------------------------------
    !write(*,*)"Entrada interp4D_Hor"
    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            x0(i) = floor(x(i))
            x1(i) = ceiling(x(i))
            y0(i) = floor(y(i))
            y1(i) = ceiling(y(i))
            z0(i) = floor(z(i))
            z1(i) = ceiling(z(i))
        endif
    end do
    !write(*,*)"interp4D_Hor 1"
    t0 = floor(t)
    t1 = ceiling(t)

    ! If depth layer has one layer
    if (n_pv == 1) then
        z0 = 1
        z1 = 1
    end if

    xd = 0.
    yd = 0.
    td = 0.
    ! Compute the "normalized coordinates" of the particle inside the data field box
    where (x1 /= x0) xd = (x-x0)/(x1-x0)
    where (y1 /= y0) yd = (y-y0)/(y1-y0)
    if (t1 /= t0) td = (t-t0)/(t1-t0)
    ! Interpolation on the first dimension and collapse it to a three dimension problem
    interp4D_Hor = 0.0
    !write(*,*)"interp4D_Hor 2"
    !write(*,*)"tamanho field : ", size(field,1), size(field,2), size(field,3), size(field,4)
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            !Use the first available value in the water column that is not 0
            bottom_lower_left_t0   = field(y0(i),x0(i),z0(i),t0)
            bottom_lower_right_t0  = field(y0(i),x1(i),z0(i),t0)
            bottom_upper_left_t0   = field(y1(i),x0(i),z0(i),t0)
            bottom_upper_right_t0  = field(y1(i),x1(i),z0(i),t0)
            bottom_lower_left_t1   = field(y0(i),x0(i),z0(i),t1)
            bottom_lower_right_t1  = field(y0(i),x1(i),z0(i),t1)
            bottom_upper_left_t1   = field(y1(i),x0(i),z0(i),t1)
            bottom_upper_right_t1  = field(y1(i),x1(i),z0(i),t1)
            ceiling_lower_left_t0  = field(y0(i),x0(i),z1(i),t0)
            ceiling_lower_right_t0 = field(y0(i),x1(i),z1(i),t0)
            ceiling_upper_left_t0  = field(y1(i),x0(i),z1(i),t0)
            ceiling_upper_right_t0 = field(y1(i),x1(i),z1(i),t0)
            ceiling_lower_left_t1  = field(y0(i),x0(i),z1(i),t1)
            ceiling_lower_right_t1 = field(y0(i),x1(i),z1(i),t1)
            ceiling_upper_left_t1  = field(y1(i),x0(i),z1(i),t1)
            ceiling_upper_right_t1 = field(y1(i),x1(i),z1(i),t1)
        
            if (bottom_lower_left_t0 == 0) then
                bottom_lower_left_t0 = ceiling_lower_left_t0
                bottom_lower_left_t1 = ceiling_lower_left_t1
            end if
            if (bottom_lower_right_t0 == 0) then
                bottom_lower_right_t0 = ceiling_lower_right_t0
                bottom_lower_right_t1 = ceiling_lower_right_t1
            end if
            if (bottom_upper_left_t0 == 0) then
                bottom_upper_left_t0 = ceiling_upper_left_t0
                bottom_upper_left_t1 = ceiling_upper_left_t1
            end if
            if (bottom_upper_right_t0 == 0) then
                bottom_upper_right_t0 = ceiling_upper_right_t0
                bottom_upper_right_t1 = ceiling_upper_right_t1
            end if

            c00(i) = bottom_lower_left_t0*(1.-xd(i)) + bottom_lower_right_t0*xd(i) !y0x0z0t0!  y0x1z0t0
            c10(i) = bottom_upper_left_t0*(1.-xd(i)) + bottom_upper_right_t0*xd(i)
            c01(i) = bottom_lower_left_t1*(1.-xd(i)) + bottom_lower_right_t1*xd(i)
            c11(i) = bottom_upper_left_t1*(1.-xd(i)) + bottom_upper_right_t1*xd(i)
        
        ! Interpolation on the second dimension and collapse it to a two dimension problem
            c0(i) = c00(i)*(1.-yd(i))+c10(i)*yd(i)
            c1(i) = c01(i)*(1.-yd(i))+c11(i)*yd(i)
        ! Interpolation on the time dimension and get the final result.
            interp4D_Hor(i) = c0(i)*(1.-td)+c1(i)*td
        endif
    end do
    !write(*,*)"Saida interp4D_Hor"
    end function interp4D_Hor

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 3d interpolation is a weighted average of 8
    !> neighbors. Consider the 4D domain between the 8 neighbors. The hypercube is
    !> divided into 4 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube.
    !> @param[in] self, x, y, t, out, field, n_fv, n_cv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp3D(self, x, y, t, out, field, n_fv, n_cv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y                        !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                         !< time to interpolate to in array coordinates
    logical, dimension(:), intent(in) :: out
    real(prec), dimension(n_fv, n_cv, n_tv), intent(in) :: field        !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv,n_tv                              !< field dimensions
    integer, intent(in) :: n_e                                          !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, x1, y1
    real(prec), dimension(n_e) :: xd, yd, c00, c10, c01, c11
    real(prec), dimension(n_e) :: c0, c1
    real(prec) :: td
    integer :: i, t0, t1
    real(prec), dimension(n_e) :: interp3D                              !< Field evaluated at x,y,z,t

    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            x0(i) = floor(x(i))
            x1(i) = ceiling(x(i))
            y0(i) = floor(y(i))
            y1(i) = ceiling(y(i))
        endif
    end do

    t0 = floor(t)
    t1 = ceiling(t)

    ! Compute the "normalized coordinates" of the particle inside the data field box
    xd = 0.
    yd = 0.
    td = 0.
    
    ! Compute the "normalized coordinates" of the particle inside the data field box
    where (x1 /= x0) xd = (x-x0)/(x1-x0)
    where (y1 /= y0) yd = (y-y0)/(y1-y0)
    if (t1 /= t0) td = (t-t0)/(t1-t0)


    interp3D = 0.0

    ! Interpolation on the first dimension and collapse it to a three dimension problem
    !do concurrent(i=1:n_e, .not. out(i))
    do i = 1, n_e
        if (.not. out(i)) then
            c00(i) = field(y0(i),x0(i),t0)*(1.-xd(i)) + field(y0(i),x1(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
            c10(i) = field(y1(i),x0(i),t0)*(1.-xd(i)) + field(y1(i),x1(i),t0)*xd(i)
            c01(i) = field(y0(i),x0(i),t1)*(1.-xd(i)) + field(y0(i),x1(i),t1)*xd(i)
            c11(i) = field(y1(i),x0(i),t1)*(1.-xd(i)) + field(y1(i),x1(i),t1)*xd(i)
        ! Interpolation on the second dimension and collapse it to a two dimension problem
            c0(i) = c00(i)*(1.-yd(i))+c10(i)*yd(i)
            c1(i) = c01(i)*(1.-yd(i))+c11(i)*yd(i)
        ! Interpolation on the time dimension and get the final result.
            interp3D(i) = c0(i)*(1.-td)+c1(i)*td
        endif
    end do

    end function interp3D

    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Returns the grid coordinates arrays of a set of tracer coordinates
    !> @param[in] self, state, bdata, time, xx, yy, zz, tt, outOfBounds, gridIsCurvilinear
    !---------------------------------------------------------------------------
    subroutine Trc2Grid_4D(self, state, bdata, time, xx, yy, zz, tt, outOfBounds, gridIsCurvilinear)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:,:), intent(in) :: state
    type(background_class), intent(in) :: bdata                 !< Background to use
    real(prec), intent(in) :: time
    real(prec), dimension(size(state,1)), intent(out) :: xx, yy, zz
    real(prec), intent(out) :: tt
    logical, dimension(:), intent(inout) :: outOfBounds
    logical, intent(in) :: gridIsCurvilinear
    integer, dimension(size(state,1), 2) :: grdCoord
    integer                              :: dim_level
    !Begin---------------------------------------------------------------
    if (gridIsCurvilinear) then
        !write(*,*)"Entrei curvilinear"
        grdCoord = self%getArrayCoord_curv(state, bdata, outOfBounds)
        xx=grdCoord(:,2) !Alterado para xx ser coluna 2
        yy=grdCoord(:,1) !Alterado para yy ser coluna 1
    else
        !write(*,*)"Entrei normais"
        xx = self%getArrayCoord(state(:,1), bdata, Globals%Var%lon, outOfBounds) !State(:,1) is always Lon
        !write(*,*)"xx(1) = ", xx(1)
        yy = self%getArrayCoord(state(:,2), bdata, Globals%Var%lat, outOfBounds)
        !write(*,*)"yy(1) = ", yy(1)
    endif
    !Get time
    tt = self%getPointCoordNonRegular(time, bdata, Globals%Var%time)
    
    !Get vertical position
    dim_level = bdata%getDimIndex(Globals%Var%level)
    if (allocated(bdata%dim(dim_level)%field4D) .and. Globals%SimDefs%inputFromHDF5) then!Means it is Verticalz from MOHID hdf5.
        zz = self%getArrayCoord(state(:,3), bdata, Globals%Var%level, outOfBounds, xx = xx, yy = yy, tt = tt)
        !write(*,*)"zz(1) = ",zz(1)
    else
        zz = self%getArrayCoord(state(:,3), bdata, Globals%Var%level, outOfBounds)
        !write(*,*)"zz(1) = ", zz(1)
    endif
    end subroutine Trc2Grid_4D
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Returns the relative grid array coordinates of a point. For curvilinear grids
    !> (or grids where lat or lon are 2D)
    !> @param[in] self, xdata, bdata, out
    !---------------------------------------------------------------------------
    function getArrayCoord_curv(self, xdata, bdata, out)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:,:), intent(in):: xdata              !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    logical, dimension(:), intent(inout) :: out
    real(prec) :: minBound_lat, minBound_lon, maxBound_lat, maxBound_lon
    real(prec) :: sumAux, perc_x, perc_y, da, db, dc, dd, x, y
    integer :: dim_lat, dim_lon, i, j, id                        !< corresponding background dimension
    real(prec), dimension(size(xdata,1),2) :: getArrayCoord_curv   !< coordinates (x,y) in array index
    type(vector), dimension(5) :: cellPolygon !Vertices of the grid cell polygon
    real(prec) :: limLeft, limRight, limBottom, limTop !limits of the polygon
    integer :: pw = 2
    logical :: foundCell
    !Begin-----------------------------------------------------------------
    dim_lon = bdata%getDimIndex(Globals%Var%lon)
    dim_lat = bdata%getDimIndex(Globals%Var%lat)
    
    minBound_lon = bdata%dim(dim_lon)%getFieldMinBound()
    where (xdata(:,1) < minBound_lon) out = .true.
    maxBound_lon = bdata%dim(dim_lon)%getFieldMaxBound()
    where (xdata(:,1) > maxBound_lon) out = .true.
    minBound_lat = bdata%dim(dim_lat)%getFieldMinBound()
    where (xdata(:,2) < minBound_lat) out = .true.
    maxBound_lat = bdata%dim(dim_lat)%getFieldMaxBound()
    where (xdata(:,2) > maxBound_lat) out = .true.
    
    !write(*,*)"tamanho da dimensao 1 da matriz longitude = ", size(bdata%dim(dim_lon)%field2D,1)
    do id = 1, size(xdata,1)
        if (.not. out(id)) then
            x = xdata(id,1)
            y = xdata(id,2)
            dj: do j = 2, size(bdata%dim(dim_lon)%field2D,2)-1
                do i = 2, size(bdata%dim(dim_lat)%field2D,1)-1
                    !Define polygon of each grid cell
                    cellPolygon(1)%x = bdata%dim(dim_lon)%field2D(i,j) !lower left corner
                    cellPolygon(1)%y = bdata%dim(dim_lat)%field2D(i,j) 
                    cellPolygon(2)%x = bdata%dim(dim_lon)%field2D(i+1,j) !upper left corner
                    cellPolygon(2)%y = bdata%dim(dim_lat)%field2D(i+1,j)
                    cellPolygon(3)%x = bdata%dim(dim_lon)%field2D(i+1,j+1) !upper right corner
                    cellPolygon(3)%y = bdata%dim(dim_lat)%field2D(i+1,j+1)
                    cellPolygon(4)%x = bdata%dim(dim_lon)%field2D(i,j+1) !lower right corner
                    cellPolygon(4)%y = bdata%dim(dim_lat)%field2D(i,j+1)
                    cellPolygon(5)%x = bdata%dim(dim_lon)%field2D(i,j) !close polygon
                    cellPolygon(5)%y = bdata%dim(dim_lat)%field2D(i,j)
                
                    limLeft = min(cellPolygon(1)%x, cellPolygon(2)%x) !limit left
                    limBottom = min(cellPolygon(1)%y, cellPolygon(4)%y) !limit bottom
                    limRight = max(cellPolygon(3)%x, cellPolygon(4)%x) !limit right
                    limTop = max(cellPolygon(2)%y, cellPolygon(3)%y) !limit top
                    
                    if((x<limLeft) .or. (x>limRight) .or. (y>limTop) .or. (y<limBottom)) cycle !Skip point
                    
                    foundCell = self%isPointInsidePolygon(x, y, cellPolygon)
                    if (foundCell) then
                        !using a IDW - Inverse distance weighting method - reference is point c
                        !(Xa,Ya) = (perc_x = 0, perc_y =1)
                        !(Xb,Yb) = (perc_x = 1, perc_y =1)        
                        !(Xc,Yc) = (perc_x = 0, perc_y =0)                
                        !(Xd,Yd) = (perc_x = 1, perc_y =0)
                        !
                        !  a_______b
                        !  |       | 
                        !  |       |         
                        !  c_______d

                        da = sqrt((x-cellPolygon(2)%x)**2+(y-cellPolygon(2)%y)**2)
                        db = sqrt((x-cellPolygon(3)%x)**2+(y-cellPolygon(3)%y)**2)
                        dc = sqrt((x-cellPolygon(1)%x)**2+(y-cellPolygon(1)%y)**2)
                        dd = sqrt((x-cellPolygon(4)%x)**2+(y-cellPolygon(4)%y)**2)        
        
                        if     (da == 0) then
                            perc_x = 0
                            perc_y = 1
                        elseif (db == 0) then
                            perc_x = 1
                            perc_y = 1
                        elseif (dc == 0) then
                            perc_x = 0
                            perc_y = 0
                        elseif (dd == 0) then
                            perc_x = 1
                            perc_y = 0
                        else
                            sumAux = ((1./da)**pw+(1./db)**pw+(1./dc)**pw+(1./dd)**pw)
                            perc_x    = ((1./db)**pw+(1./dd)**pw) / sumAux
                            perc_y    = ((1./da)**pw+(1./db)**pw) / sumAux            
                        endif
                        getArrayCoord_curv(id,1) = j + j*perc_x !relative position to J index
                        getArrayCoord_curv(id,2) = i + i*perc_y !relative position to I index
                        exit dj
                    endif
                enddo
            enddo dj
        endif
    enddo

    end function getArrayCoord_curv
    
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Imported from MOHID horizontalgrid module
    !> The first action performed is an acceleration test. if 
    !>                   point is not inside maximum limits of the polygon
    !>                   then it is not inside the polygon. If point is inside
    !>                   the maximum bounds of the polygon then there is a 
    !>                   possibility of being inside the polygon. Thus, it is 
    !>                   drawn a "semi-recta" on the X axis to the right of the
    !>                   point. If the "semi-recta" intersects the polygon an 
    !>                   odd ("ímpar") number of times then the point is 
    !>                   inside the polygon. The intersection test is performed 
    !>                   for every segment of the polygon. If the point belongs 
    !>                   is a vertix of the polygon or belongs to one of the 
    !>                   segments that defines the polygon then it is considered 
    !>                   to be inside the polygon. There are a few cases in which
    !>                   special care needs to be taken, regarding the intersection
    !>                   counting. If the intersection point coincides with one of 
    !>                   the vertices then one intersection is counted only if it is
    !>                   the vertice with the higher Y value. This is the way to 
    !>                   count only one intersection once one polygon vertix belongs
    !>                   to two segments and otherwise it would be tested twice.
    !>                   Another case is when the segment belongs to the "semi-recta".
    !>                   In this case no intersection is counted. This in the next 
    !>                   segment to be tested the "semi-recta" will intersect the first
    !>                   vertice of the segment. If the segments slope is negative then
    !>                   the vertice has the higher Y value and one intersection is 
    !>                   counted.
    !> @param[in] self, x, y, polygon
    logical function isPointInsidePolygon(self, x, y, polygon)
    class(interpolator_class), intent(in) :: self
    real(prec), intent(in) :: x, y
    type(vector), dimension(5), intent(in) :: polygon !Vertices of the grid cell polygon
    integer :: i
    real(prec) :: segStart_x, segEnd_x, segStart_y, segEnd_y
    integer :: numberOfIntersections
    real(prec) :: slope, intersectionPointX, higherY
    !Begin-----------------------------------------------------------------

    isPointInsidePolygon  = .false.
    NumberOfIntersections = 0

    do i = 1, 4 !Go through all vertices
        !construct segment
        segStart_x = polygon(i)%x
        segStart_y = polygon(i)%y
        segEnd_x = polygon(i+1)%x
        segEnd_y = polygon(i+1)%y
            
        !if point coincides with one of the segments vertices then
        !it in inside the polygon
        if(((x == segStart_x) .and. (y == segStart_y)) .or. ((x == segEnd_x) .and. (y == segEnd_y))) then
            isPointInsidePolygon = .true.
            return
        !if point is placed between the segment vertices in the Y axis
        elseif(((y <= segStart_y) .and. (y >= segEnd_y)) .or. ((y >= segStart_y) .and. (y <= segEnd_y)))then
            !if segment has a slope
            if(segStart_y /= segEnd_y)then 
                !compute slope
                slope = (segStart_x - segEnd_x) / (segStart_y - segEnd_y)
                !compute intersection point X coordinate of the "semi-recta" with the segment 
                intersectionPointX = segStart_x - slope * (segStart_y - y)
                !if point belongs to the segment then it is inside the polygon
                if(intersectionPointX == x)then
                    isPointInsidePolygon = .true.
                    return
                elseif(intersectionPointX > x)then
                    !if the intersection point coincides with one of the vertices then one intersection
                    !is counted only if it is the vertice with the higher Y value 
                    if((y == segStart_y) .or. (y == segEnd_y)) then
                        !find higher segment ends higher Y coordinate
                        higherY = max(segStart_y, segEnd_y)
                        if(y == higherY) numberOfIntersections = numberOfIntersections + 1
                    else
                        !if the intersection point is placed to the right of the point or
                        !if the point belongs to the segment then one intersection is counted
                        numberOfIntersections = numberOfIntersections + 1
                    end if
                end if
            elseif(segStart_y == segEnd_y)then
                if(y == segStart_y)then
                    if(segStart_x > segEnd_x)then
                        if((x <= segStart_x) .and. (x >= segEnd_x))then
                            !if point belongs to the segment then it is inside the polygon
                            isPointInsidePolygon = .true.
                            return
                        end if
                    elseif(segStart_x < segEnd_x)then
                        if((x >= segStart_x) .and. (x <= segEnd_x))then
                            !if point belongs to the segment then it is inside the polygon
                            isPointInsidePolygon = .true.
                            return
                        end if
                    elseif(segStart_x == segEnd_x)then
                        if(x == segStart_x)then
                            !if point belongs to the segment then it is inside the polygon
                            isPointInsidePolygon = .true.
                            return
                        end if
                    end if                       
                end if
            end if
        end if
    enddo

    !if number of intersections is odd (odd = ímpar) then
    !point is inside the polygon
    if(numberOfIntersections/2. > int(numberOfIntersections/2.)) isPointInsidePolygon = .true.

    end function isPointInsidePolygon
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinates of a set of points, given a coordinate
    !> array.
    !> @param[in] self, xdata, bdata, dimName, out, bat, xx, yy, tt
    !---------------------------------------------------------------------------
    function getArrayCoord(self, xdata, bdata, dimName, out, bat, xx, yy, tt)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    type(string), intent(in) :: dimName
    logical, dimension(:), intent(inout) :: out
    real(prec), dimension(:), optional, intent(in) :: bat
    real(prec), dimension(:), optional, intent(in):: xx, yy
    real(prec), optional, intent(in):: tt
    integer :: dim                                              !< corresponding background dimension
    real(prec), dimension(size(xdata)) :: getArrayCoord         !< coordinates in array index
    dim = bdata%getDimIndex(dimName)
    if (bdata%regularDim(dim)) then
        if (present(xx)) then
            getArrayCoord = self%getArrayCoordRegular(xdata, bdata, dim, out, xx = xx, yy = yy, tt = tt)
        else
            getArrayCoord = self%getArrayCoordRegular(xdata, bdata, dim, out)
            !write(*,*)"getArrayCoord(1) = ", getArrayCoord(1)
        endif
        
    endif
    
    if (.not.bdata%regularDim(dim)) then
        if (present(bat)) then
            getArrayCoord = self%getArrayCoordNonRegular(xdata, bdata, dim, out, bat=bat)
        else
            if (present(xx)) then
                !VerticalZ
                !write(*,*)"Entrei .not.bdata%regularDim(dim) sem bat e com xx", trim(dimName)
                getArrayCoord = self%getArrayCoordNonRegular(xdata, bdata, dim, out, xx = xx, yy = yy, tt = tt)
            else
                getArrayCoord = self%getArrayCoordNonRegular(xdata, bdata, dim, out)
            endif
        end if
    end if 

    end function getArrayCoord

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> V2 in Dec 2022 by Joao Sobrinho
    !> Returns the array coordinates of a set of points, given a coordinate
    !> array. Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dim, out
    !---------------------------------------------------------------------------
    function getArrayCoordRegular(self, xdata, bdata, dim, out, xx, yy, tt)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    integer, intent(in) :: dim
    logical, dimension(:), intent(inout) :: out
    real(prec), dimension(:), optional, intent(in):: xx, yy
    real(prec), optional, intent(in):: tt
    real(prec), dimension(size(xdata)) :: getArrayCoordRegular  !< coordinates in array index
    real(prec) :: minBound, maxBound, res
    integer :: fieldLength, dimSize, indx, i
    !Begin----------------------------------------------------------
    !write(*,*)"Entrei getArrayCoordRegular"
    if (allocated(bdata%dim(dim)%field1D)) then !dimension variable is 1D
        fieldLength = size(bdata%dim(dim)%field1D)
        dimSize = 1
        if (fieldLength == 1) then
            getArrayCoordRegular = 1
            return
        endif
    elseif (allocated(bdata%dim(dim)%field2D)) then
        dimSize = 2
    endif
    !write(*,*)"dimsize =  ", dimSize
    if (dimSize == 1) then
        !write(*,*)"Entrei no dimsize1 ", trim(bdata%dim(dim)%name)
        minBound = bdata%dim(dim)%getFieldMinBound()
        maxBound = bdata%dim(dim)%getFieldMaxBound()
        res = abs(maxBound - minBound)/(size(bdata%dim(dim)%field1D)-1.0)
        getArrayCoordRegular = (xdata - minBound)/res + 1.0
    else
        !write(*,*)"Entrei no dimsize2 ", trim(bdata%dim(dim)%name)
        minBound = bdata%dim(dim)%getFieldMinBound()
        !write(*,*)"minBound = ", minBound
        maxBound = bdata%dim(dim)%getFieldMaxBound()
        !write(*,*)"maxBound = ", maxBound
        !write(*,*)"tamanho = ", size(bdata%dim(dim)%field2D,indx)
        if (bdata%dim(dim)%name == Globals%Var%lat) indx = 1 !lat varies with index 2 of field2D
        if (bdata%dim(dim)%name == Globals%Var%lon) indx = 2 !lon varies with index 2 of field2D
        res = abs(maxBound - minBound)/(size(bdata%dim(dim)%field2D,indx)-1.0)
        !write(*,*)"res = ", res
        !write(*,*)"xdata(1) = ", xdata(1)
        !write(*,*)"xdata - minBound = ", xdata(1) - minBound
        !write(*,*)"getArrayCoordRegular 1 = ", (xdata(1) - minBound)/res
        !write(*,*)"getArrayCoordRegular 2 = ", (xdata(1) - minBound)/res + 1.0
        getArrayCoordRegular = (xdata - minBound)/res + 1.0 !this will give the index along the dimension indx which intersects with the xdata
        !do i = 1, size(xdata)
        !    write(*,*)"getArrayCoordRegular (i) = ", getArrayCoordRegular(i)
        !enddo
        
    endif
    !write(*,*)"out 1 = ", out(1)
    where (xdata < minBound) out = .true.
    where (xdata > maxBound) out = .true.
    !write(*,*)"out 1 saida = ", out(1)
    end function getArrayCoordRegular
    
    
    !!---------------------------------------------------------------------------
    !!> @author Ricardo Birjukovs Canelas - MARETEC
    !!> @brief
    !!> Returns the array coordinates of a set of points, given a coordinate
    !!> array. Works only for regularly spaced data.
    !!> @param[in] self, xdata, bdata, dim, out
    !!---------------------------------------------------------------------------
    !function getArrayCoordRegular(self, xdata, bdata, dim, out)
    !class(interpolator_class), intent(in) :: self
    !real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    !type(background_class), intent(in) :: bdata                 !< Background to use
    !integer, intent(in) :: dim
    !logical, dimension(:), intent(inout) :: out
    !real(prec), dimension(size(xdata)) :: getArrayCoordRegular  !< coordinates in array index
    !real(prec) :: minBound, maxBound, res
    ! 
    !if(size(bdata%dim(dim)%field) == 1) then
    !    getArrayCoordRegular = 1
    !    return
    !end if
    !minBound = bdata%dim(dim)%getFieldMinBound()
    !maxBound = bdata%dim(dim)%getFieldMaxBound()
    !res = abs(maxBound - minBound)/(size(bdata%dim(dim)%field)-1.0)
    !getArrayCoordRegular = (xdata - minBound)/res + 1.0
    !where (xdata < minBound) out = .true.
    !where (xdata > maxBound) out = .true.
    !end function getArrayCoordRegular

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> @param[in] self, xdata, bdata, dim, out, bat, xx, yy, tt
    ! !---------------------------------------------------------------------------
    function getArrayCoordNonRegular(self, xdata, bdata, dim, out, bat, xx, yy, tt)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                    !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                     !< Background to use
    integer, intent(in) :: dim
    logical, dimension(:), intent(inout) :: out
    real(prec), dimension(:), optional, intent(in) :: bat
    real(prec), dimension(:), optional, intent(in):: xx, yy
    real(prec), optional, intent(in):: tt
    integer :: i, j, k                              
    integer :: id, idx_1, idx_2, indx, fieldLength, dimSize, t1, t2                             
    real(prec), dimension(size(xdata)) :: getArrayCoordNonRegular   !< coordinates in array index
    real(prec) :: minBound, maxBound, positionT1, positionT2, distance
    type(string) :: outext
    !Begin--------------------------------------------------------------------------    
    if (allocated(bdata%dim(dim)%field1D)) then
        fieldLength = size(bdata%dim(dim)%field1D)
        dimSize = 1
        if (fieldLength == 1) then
            getArrayCoordNonRegular = 1
            return
        endif
    elseif (allocated(bdata%dim(dim)%field2D)) then
        dimSize = 2
    elseif (allocated(bdata%dim(dim)%field4D)) then
        dimSize = 4 !For VerticalZ
    endif
    
    getArrayCoordNonRegular = 1
    
    if (dimSize == 1) then
        minBound = bdata%dim(dim)%getFieldMinBound()
        maxBound = bdata%dim(dim)%getFieldMaxBound()
        where (xdata < minBound) out = .true.
        where (xdata > maxBound) out = .true.
        if (present(bat)) then
            do id = 1,size(xdata)
                if (.not. out(id)) then
                    do i = 2, size(bdata%dim(dim)%field1D)
                        if (bdata%dim(dim)%field1D(i) >= xdata(id)) then
                            idx_1 = i-1
                            idx_2 = i
                            exit
                        end if
                    end do
                    !Bathymetric value is deeper than bottom face of layer where tracer is located
                    if (bat(id) <= bdata%dim(dim)%field1D(idx_1)) then
                       getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field1D(idx_1))/(bdata%dim(dim)%field1D(idx_2)-bdata%dim(dim)%field1D(idx_1))) 
                    else
                        !Bathymetry is inside bottom open cell. using max to avoid error when tracer is below the bathymetry (its position is corrected afterwards)
                        getArrayCoordNonRegular(id) = idx_1 + abs((max(xdata(id),bat(id))-bat(id))/(bdata%dim(dim)%field1D(idx_2)-bat(id))) 
                    end if
                endif
            enddo
        else
            do id = 1,size(xdata)
                if (.not. out(id)) then
                    do i = 2, size(bdata%dim(dim)%field1D)
                        if (bdata%dim(dim)%field1D(i) >= xdata(id)) then
                            idx_1 = i-1
                            idx_2 = i
                            exit
                        end if
                    end do
                    getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field1D(idx_1))/(bdata%dim(dim)%field1D(idx_2)-bdata%dim(dim)%field1D(idx_1)))
                endif
            enddo
        end if
    elseif (dimSize == 2) then !2D lat and lon grid (not yet ready for a 3D vertical layer dimension
        if (present(bat)) then
            outext = '[Interpolator::getArrayCoordNonRegular] correction with bathymetry cannot yet be done when level dimension is not 1D, stoping'
            call Log%put(outext)
            stop
        endif
        
        minBound = bdata%dim(dim)%getFieldMinBound()
        maxBound = bdata%dim(dim)%getFieldMaxBound()
        where (xdata < minBound) out = .true.
        where (xdata > maxBound) out = .true.
        
        if (bdata%dim(dim)%name == Globals%Var%lat)indx = 1
        if (bdata%dim(dim)%name == Globals%Var%lon)indx = 2
        
        if (indx == 1) then!Lat
            do id = 1,size(xdata)
                if (.not. out(id)) then
                    do i = 2, size(bdata%dim(dim)%field2D, indx)
                        if (bdata%dim(dim)%field2D(i,1) >= xdata(id)) then
                            idx_1 = i-1
                            idx_2 = i
                            exit
                        end if
                    end do
                    getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field2D(idx_1,1))/(bdata%dim(dim)%field2D(idx_2,1)-bdata%dim(dim)%field2D(idx_1,1)))
                end if
            enddo
        else
            do id = 1,size(xdata)
                if (.not. out(id)) then
                    do i = 2, size(bdata%dim(dim)%field2D, indx)
                        if (bdata%dim(dim)%field2D(1,i) >= xdata(id)) then
                            idx_1 = i-1
                            idx_2 = i
                            exit
                        end if
                    end do
                    getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field2D(1,idx_1))/(bdata%dim(dim)%field2D(1,idx_2)-bdata%dim(dim)%field2D(1,idx_1)))
                endif
            enddo
        endif
    elseif (dimSize == 4 .and. present(xx)) then
        !Get position for instant 1
        t1 = floor(tt)
        t2 = ceiling(tt)
        
        !do concurrent(id = 1:size(xdata), .not. out(id))
        do id = 1,size(xdata)
            if (.not. out(id)) then
                i = floor(yy(id))
                j = floor(xx(id))
            
                !write(*,*)"Depth k =1 = ", bdata%dim(dim)%field4D(i, j, 1, t1)
                !write(*,*)"Depth k =2 = ", bdata%dim(dim)%field4D(i, j, 2, t1)
                !write(*,*)"Depth k =3 = ", bdata%dim(dim)%field4D(i, j, 3, t1)
                idx_1 = size(bdata%dim(dim)%field4D,3)
                idx_2 = idx_1
                do k = 2, size(bdata%dim(dim)%field4D, 3)
                    if (bdata%dim(dim)%field4D(i, j, k, t1) >= xdata(id)) then
                        idx_1 = k-1
                        idx_2 = k
                        exit
                    end if
                end do
            
                distance = (bdata%dim(dim)%field4D(i,j,idx_2,t1)-bdata%dim(dim)%field4D(i,j,idx_1,t1))
                if (distance == 0 .or. xdata(id) < bdata%dim(dim)%field4D(i, j, idx_1, t1)) then
                    !Also use idx_1 if xdata is outside the vertical bounds
                    positionT1 = idx_1
                else
                    positionT1 = idx_1 + abs((xdata(id)-bdata%dim(dim)%field4D(i,j,idx_1,t1))/distance)
                endif
            
                !Get position for instant 2
                if (t1 /= t2) then
                    positionT2 = positionT1
                    idx_1 = size(bdata%dim(dim)%field4D,3)
                    idx_2 = idx_1
                    do k = 2, size(bdata%dim(dim)%field4D, 3)
                        if (bdata%dim(dim)%field4D(i, j, k, t2) >= xdata(id)) then
                            idx_1 = k-1
                            idx_2 = k
                            exit
                        end if
                    end do
                
                    distance = (bdata%dim(dim)%field4D(i,j,idx_2,t2)-bdata%dim(dim)%field4D(i,j,idx_1,t2))
                    if (distance == 0 .or. xdata(id) < bdata%dim(dim)%field4D(i, j, idx_1, t2)) then
                        !Also use idx_1 if xdata is outside the vertical bounds
                        positionT2 = idx_1
                    else
                        positionT2 = idx_1 + abs((xdata(id)-bdata%dim(dim)%field4D(i,j,idx_1,t2))/distance)
                    endif
                    getArrayCoordNonRegular(id) = (positionT1 + positionT2) / 2
                else
                    getArrayCoordNonRegular(id) = positionT1
                endif
            endif
        end do
    endif

    end function getArrayCoordNonRegular
    
    !!> @author Daniel Garaboa Paz - USC
    !!> @brief
    !!> Returns the array coordinate of a point, along a given dimension.
    !!> @param[in] self, xdata, bdata, dim, out
    !! !---------------------------------------------------------------------------
    !function getArrayCoordNonRegular(self, xdata, bdata, dim, out, bat)
    !class(interpolator_class), intent(in) :: self
    !real(prec), dimension(:), intent(in):: xdata                    !< Tracer coordinate component
    !type(background_class), intent(in) :: bdata                     !< Background to use
    !integer, intent(in) :: dim
    !logical, dimension(:), intent(inout) :: out
    !real(prec), dimension(:), optional, intent(in) :: bat
    !integer :: i                                                
    !integer :: id, idx_1, idx_2                                 
    !real(prec), dimension(size(xdata)) :: getArrayCoordNonRegular   !< coordinates in array index
    !real(prec) :: minBound, maxBound
    !
    !if(size(bdata%dim(dim)%field) == 1) then
    !    getArrayCoordNonRegular = 1
    !    return
    !end if
    !getArrayCoordNonRegular = 1
    !minBound = bdata%dim(dim)%getFieldMinBound()
    !maxBound = bdata%dim(dim)%getFieldMaxBound()
    !where (xdata < minBound) out = .true.
    !where (xdata > maxBound) out = .true.
    !if (present(bat)) then
    !    do concurrent(id = 1:size(xdata), .not. out(id))
    !        do i = 2, size(bdata%dim(dim)%field)
    !            if (bdata%dim(dim)%field(i) >= xdata(id)) then
    !                idx_1 = i-1
    !                idx_2 = i
    !                exit
    !            end if
    !        end do
    !        !Bathymetric value is deeper than bottom face of layer where tracer is located
    !        if (bat(id) <= bdata%dim(dim)%field(idx_1)) then
    !           getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field(idx_1))/(bdata%dim(dim)%field(idx_2)-bdata%dim(dim)%field(idx_1))) 
    !        else
    !            !Bathymetry is inside bottom open cell. using max to avoid error when tracer is below the bathymetry (its position is corrected afterwards)
    !            getArrayCoordNonRegular(id) = idx_1 + abs((max(xdata(id),bat(id))-bat(id))/(bdata%dim(dim)%field(idx_2)-bat(id))) 
    !        end if
    !        
    !    end do
    !else
    !    do concurrent(id = 1:size(xdata), .not. out(id))
    !        do i = 2, size(bdata%dim(dim)%field)
    !            if (bdata%dim(dim)%field(i) >= xdata(id)) then
    !                idx_1 = i-1
    !                idx_2 = i
    !                exit
    !            end if
    !        end do
    !        getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field(idx_1))/(bdata%dim(dim)%field(idx_2)-bdata%dim(dim)%field(idx_1)))
    !    end do
    !end if
    !end function getArrayCoordNonRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName
    !---------------------------------------------------------------------------
    function getPointCoordNonRegular(self, xdata, bdata, dimName)
    class(interpolator_class), intent(in) :: self
    real(prec), intent(in):: xdata              !< Tracer coordinate component
    type(background_class), intent(in) :: bdata !< Background to use
    type(string), intent(in) :: dimName
    integer :: dim                              !< corresponding background dimension
    real(prec) :: getPointCoordNonRegular       !< coordinates in array index
    type(string) :: outext
    integer :: i
    integer :: idx_1, idx_2, n_idx
    logical :: found

    found = .false.
    dim = bdata%getDimIndex(dimName)
    if (.not. allocated(bdata%dim(dim)%field1D)) then
        outext = '[Interpolator::getPointCoordNonRegular] variable "'//dimName//'" needs to be 1d. stoping'
        call Log%put(outext)
        stop
    endif
    
    if(size(bdata%dim(dim)%field1D) == 1) then
        getPointCoordNonRegular = 1
        return
    end if
    
    n_idx = size(bdata%dim(dim)%field1D)
    do i = 2, n_idx
        if (bdata%dim(dim)%field1D(i) >= xdata) then
            idx_1 = i-1
            idx_2 = i
            found = .true.
            exit
        end if
    end do
    if (.not.found) then
        outext = '[Interpolator::getPointCoordNonRegular] Point not contained in "'//dimName//'" dimension, stoping'
        call Log%put(outext)
        stop
    end if
    getPointCoordNonRegular = idx_1 + abs((xdata-bdata%dim(dim)%field1D(idx_1))/(bdata%dim(dim)%field1D(idx_2)-bdata%dim(dim)%field1D(idx_1)))
    end function getPointCoordNonRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName, eta
    !---------------------------------------------------------------------------
    function getPointCoordRegular(self, xdata, bdata, dimName, eta)
    class(interpolator_class), intent(in) :: self
    real(prec), intent(in):: xdata              !< Tracer coordinate component
    type(background_class), intent(in) :: bdata !< Background to use
    type(string), intent(in) :: dimName
    real(prec), intent(in), optional :: eta
    integer :: dim                              !< corresponding background dimension
    real(prec) :: getPointCoordRegular          !< coordinates in array index
    real(prec) :: minBound, maxBound, res, ieta
    type(string) :: outext
    integer :: id                                 !< corresponding background dimension
    dim = bdata%getDimIndex(dimName)
    if (.not. allocated(bdata%dim(dim)%field1D)) then
        outext = '[Interpolator::getPointCoordRegular] variable "'//dimName//'" needs to be 1d. stoping'
        call Log%put(outext)
        stop
    endif
    res = size(bdata%dim(dim)%field1D)-1
    minBound = bdata%dim(dim)%getFieldMinBound()
    maxBound = bdata%dim(dim)%getFieldMaxBound()
    res = abs(maxBound - minBound)/res
    getPointCoordRegular = (xdata - minBound)/res+1
    ieta = -res/10.0
    if (present(eta)) ieta = eta
    if (.not.Utils%isBounded(xdata, minBound, maxBound, ieta)) then
        outext = '[Interpolator::getPointCoordRegular] Point not contained in "'//dimName//'" dimension, stoping'
        print*, xdata
        print*, minBound, maxBound+ieta
        call Log%put(outext)
        stop
    end if
    end function getPointCoordRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializer method for the Interpolator class. Sets the type of interpolator
    !> and name of the algorithm this Interpolator will call
    !> @param[in] self, flag, name
    !---------------------------------------------------------------------------
    subroutine initInterpolator(self, flag, name)
    class(interpolator_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    self%interpType = flag
    self%name = name
    end subroutine initInterpolator

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Test for interp 4D
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine test4D(self)
    class(interpolator_class), intent(inout) :: self
    real(prec), dimension(:,:,:,:), allocatable :: field
    real(prec), dimension(:), allocatable :: xx, yy, zz
    logical, dimension(:), allocatable :: out
    real(prec) :: time
    integer :: npts, fieldDims
    real(prec) :: fieldVal
    fieldDims = 10
    fieldVal = 1.9
    allocate(field(fieldDims,fieldDims,fieldDims,fieldDims))
    field = fieldVal
    npts = 1
    allocate(xx(npts), yy(npts), zz(npts))
    allocate(out(npts))
    xx = 13.45
    yy = xx
    zz = xx
    time = 1
    print*, 'testing 4D interpolation, expected result is ', fieldVal
    out = .false.
    xx = self%interp4D(xx, yy, zz, time, out, field, fieldDims, fieldDims, fieldDims, fieldDims, npts)
    print*, 'result = ', xx
    if (xx(1) == fieldVal) then
        print*, 'Test: SUCCESS'
    else
        print*, 'Test: FAILED'
    end if
    read(*,*)
    end subroutine test4D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the Interpolator information
    !---------------------------------------------------------------------------
    subroutine printInterpolator(self)
    class(interpolator_class), intent(inout) :: self
    type(string) :: outext
    outext = 'Interpolation algorithm is '//self%name
    call Log%put(outext,.false.)
    end subroutine printInterpolator

    end module interpolator_mod