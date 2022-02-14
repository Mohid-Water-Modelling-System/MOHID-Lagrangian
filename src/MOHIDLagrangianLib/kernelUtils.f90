    module kernelUtils_mod
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernel
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Defines a kernel utilities class.
    !> 
    !------------------------------------------------------------------------------
    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod

    type :: kernelUtils_class        !< Kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelUtils
    procedure :: getInterpolatedFields => runInterpolatorOnVars
    procedure, private :: getBackgroundDictIntersection
    end type kernelUtils_class
     
    public :: kernelUtils_class
    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> calls interpolator for each required var. Seaches in every background
    !---------------------------------------------------------------------------
    subroutine runInterpolatorOnVars(self, sv, bdata, time, requiredVars, var_dt, var_name, justRequired)
    class(kernelUtils_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), dimension(:,:), allocatable, intent(inout) :: var_dt
    type(string), dimension(:), allocatable, intent(inout) :: var_name
    type(string), dimension(:), intent(in) :: requiredVars
    logical, optional, intent(in) :: justRequired
    logical :: localjustRequired
    real(prec), intent(in) :: time
    integer :: np, nf, bkg, i, j
    type(varBackground_t), dimension(:), allocatable :: bkgToInterpolate
    type(string), allocatable, dimension(:) :: currList
    
    call self%getBackgroundDictIntersection(requiredVars, bkgToInterpolate)
    
    localjustRequired = .false.
    if (present(justRequired)) then
        localjustRequired = justRequired
    end if
    
    np = size(sv%active) !number of Tracers
    
    if (localjustRequired) then
        allocate(var_dt(np,size(requiredVars)))
        allocate(var_name(size(requiredVars)))    
        j=1
        do i=1, size(bkgToInterpolate)
            call bkgToInterpolate(i)%bkgVars%toArray(currList)
            call self%Interpolator%run(sv%state, bdata(bkgToInterpolate(i)%bkgIndex), time, var_dt(:,j:j+bkgToInterpolate(i)%bkgVars%getSize()-1), var_name(j:j+bkgToInterpolate(i)%bkgVars%getSize()-1), currList)
            j = j + bkgToInterpolate(i)%bkgVars%getSize()
            deallocate (currList)
        end do
    else
        nf = 0
        do i=1, size(bkgToInterpolate)
            nf = nf + bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize() !number of fields to interpolate
        end do
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        j=1
        do i=1, size(bkgToInterpolate)
            call self%Interpolator%run(sv%state, bdata(bkgToInterpolate(i)%bkgIndex), time, var_dt(:,j:j+bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()-1), var_name(j:j+bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()-1))
            j = j + bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()
        end do
    end if
    
    end subroutine runInterpolatorOnVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Gets intersection between required var and background
    !---------------------------------------------------------------------------
    subroutine getBackgroundDictIntersection(self, reqVars, bkgVarsIntersect)
    class(kernelUtils_class), intent(inout) :: self
    type(string), dimension(:), intent(in) :: reqVars
    type(varBackground_t), dimension(:), allocatable, intent(out) :: bkgVarsIntersect
    integer, dimension(:), allocatable :: bkgVarsIntersectNumber
    integer :: i, j, k
    !To compensate for openMP bug :)
    type(varBackground_t), dimension(:), allocatable :: localBackgroundVarDict
    
    allocate(localBackgroundVarDict(size(Globals%BackgroundVarDict%bkgDict)), source=Globals%BackgroundVarDict%bkgDict)
     
    allocate(bkgVarsIntersectNumber(size(Globals%BackgroundVarDict%bkgDict)))
    bkgVarsIntersectNumber = 0
    
    do i=1, size(localBackgroundVarDict)
        do j=1, size(reqVars)
            if (.not.localBackgroundVarDict(i)%bkgVars%notRepeated(reqVars(j))) then
                bkgVarsIntersectNumber(i) = bkgVarsIntersectNumber(i) + 1
            end if
        end do
    end do
    
    i = 0
    do j =1, size(bkgVarsIntersectNumber)
        if(bkgVarsIntersectNumber(j)>0) then
            i = i+1
        end if
    end do

    allocate(bkgVarsIntersect(i))
    k = 1
    do i=1, size(localBackgroundVarDict)
        do j=1, size(reqVars)
            if (.not.localBackgroundVarDict(i)%bkgVars%notRepeated(reqVars(j))) then
                bkgVarsIntersect(k)%bkgIndex = i
                if (bkgVarsIntersect(k)%bkgVars%notRepeated(reqVars(j))) then
                    call bkgVarsIntersect(k)%bkgVars%add(reqVars(j))
                    if (bkgVarsIntersect(k)%bkgVars%getSize() == bkgVarsIntersectNumber(i)) then
                        k = k + 1
                    end if
                end if
            end if
        end do
    end do    
    end subroutine getBackgroundDictIntersection

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelUtils(self)
    class(kernelUtils_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernelUtils

    end module kernelUtils_mod