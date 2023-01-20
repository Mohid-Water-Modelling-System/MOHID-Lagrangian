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
    subroutine runInterpolatorOnVars(self, sv, bdata, time, requiredVars, var_dt, var_name, justRequired, reqVertInt)
    class(kernelUtils_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), dimension(:,:), allocatable, intent(inout) :: var_dt
    type(string), dimension(:), allocatable, intent(inout) :: var_name
    type(string), dimension(:), intent(in) :: requiredVars
    logical, optional, intent(in) :: justRequired, reqVertInt
    logical :: localjustRequired, localreqVertInt
    real(prec), intent(in) :: time
    integer :: np, nf, bkg, i, j, aux
    type(varBackground_t), dimension(:), allocatable :: bkgToInterpolate
    type(string), allocatable, dimension(:) :: currList
    !begin--------------------------------------------------------------------------
    !write(*,*)"Entrada runInterpolatorOnVars"
    call self%getBackgroundDictIntersection(requiredVars, bkgToInterpolate)
    !write(*,*)"Tamanho bkgToInterpolate = ", size(bkgToInterpolate)
    !do i=1,size(bkgToInterpolate)
    !    write(*,*)"index = ", bkgToInterpolate(i)%bkgIndex
    !    aux = bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()
    !    write(*,*)"tamanho bdata(bkgToInterpolate(i)%bkgIndex) = ", aux
    !enddo
    
    localjustRequired = .false.
    if (present(justRequired)) then
        localjustRequired = justRequired
    end if
    localreqVertInt = .true.
    if (present(reqVertInt)) then
        localreqVertInt = reqVertInt
    end if
    
    np = size(sv%active) !number of Tracers
    
    if (localjustRequired) then
        !write(*,*)"localjustRequired = true"
        allocate(var_dt(np,size(requiredVars)))
        allocate(var_name(size(requiredVars)))    
        j=1
        do i=1, size(bkgToInterpolate)
            !write(*,*)"i do bkgToInterpolate = ", i
            !write(*,*)"j do bkgToInterpolate = ", j
            call bkgToInterpolate(i)%bkgVars%toArray(currList)
            !write(*,*)"runInterpolatorOnVars 3 : indice do background = ", bkgToInterpolate(i)%bkgIndex
            call self%Interpolator%run(sv%state, bdata(bkgToInterpolate(i)%bkgIndex), time, var_dt(:,j:j+bkgToInterpolate(i)%bkgVars%getSize()-1), var_name(j:j+bkgToInterpolate(i)%bkgVars%getSize()-1), currList, reqVertInt=localreqVertInt)
            !write(*,*)"runInterpolatorOnVars 4 : saida :", bkgToInterpolate(i)%bkgVars%getSize()
            j = j + bkgToInterpolate(i)%bkgVars%getSize()
            deallocate (currList)
        end do
    else
        !write(*,*)"localjustRequired = false"
        nf = 0
        do i=1, size(bkgToInterpolate)
            aux = bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()
            !write(*,*)"tamanho bdata(bkgToInterpolate(i)%bkgIndex) = ", aux
            nf = nf + bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize() !number of fields to interpolate
        end do
        !write(*,*)"tamanho np = ", np
        !write(*,*)"tamanho nf = ", nf
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        j=1
        !write(*,*)"tamanho bkgToInterpolate = ", size(bkgToInterpolate)
        do i=1, size(bkgToInterpolate)
            !write(*,*)"i do bkgToInterpolate = ", i
            !write(*,*)"j do bkgToInterpolate = ", j
            aux = bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()
            !write(*,*)"tamanho 2 bdata(bkgToInterpolate(i)%bkgIndex) = ", aux
            call self%Interpolator%run(sv%state, bdata(bkgToInterpolate(i)%bkgIndex), time, var_dt(:,j:j+bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()-1), var_name(j:j+bdata(bkgToInterpolate(i)%bkgIndex)%fields%getSize()-1), reqVertInt=localreqVertInt)
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
    !write(*,*)"tamanho localBackgroundVarDict = ", size(localBackgroundVarDict)
    !write(*,*)"tamanho reqVars = ", size(reqVars)
    do i=1, size(localBackgroundVarDict)
        !write(*,*)"localBackgroundVarDict i = ", i
        do j=1, size(reqVars)
            !write(*,*)"localBackgroundVarDict j = ", j
            !write(*,*)"variavel = ", trim(reqVars(j))
            !call localBackgroundVarDict(i)%bkgVars%print()
            if (.not.localBackgroundVarDict(i)%bkgVars%notRepeated(reqVars(j))) then
                !write(*,*)"encontrei = ", trim(reqVars(j))
                bkgVarsIntersectNumber(i) = bkgVarsIntersectNumber(i) + 1
            end if
        end do
        !write(*,*)"items encontrados no localBackgroundVarDict = ", i, bkgVarsIntersectNumber(i)
    end do
    
    i = 0
    do j =1, size(bkgVarsIntersectNumber)
        !write(*,*)"j entrada = ", j
        !write(*,*)"i entrada = ", i
        if(bkgVarsIntersectNumber(j)>0) then
            i = i+1
        end if
        !write(*,*)"j saida = ", j
        !write(*,*)"i saida = ", i
    end do

    allocate(bkgVarsIntersect(i))
    k = 1
    do i=1, size(localBackgroundVarDict)
        !write(*,*)"i entrada localBackgroundVarDict = ", i
        do j=1, size(reqVars)
            !write(*,*)"j entrada localBackgroundVarDict = ", j
            !write(*,*)"k entrada localBackgroundVarDict = ", k
            !write(*,*)"required var = ", reqVars(j)
            !call localBackgroundVarDict(i)%bkgVars%print()
            if (.not.localBackgroundVarDict(i)%bkgVars%notRepeated(reqVars(j))) then
                !write(*,*)"encontrei reqVars(j) no localBackgroundVarDict"
                bkgVarsIntersect(k)%bkgIndex = i
                if (bkgVarsIntersect(k)%bkgVars%notRepeated(reqVars(j))) then
                    !write(*,*)"encontrei reqVars(j) no bkgVarsIntersect"
                    call bkgVarsIntersect(k)%bkgVars%add(reqVars(j))
                    !if number of variables in bkgVarsIntersect is the same as the number of variables identified in a background
                    !add 1 to k
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