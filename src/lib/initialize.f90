
    module initialize

    use tracer3D
    use simulation_parameters
    use source_identity

    use FoX_dom
    use commom_modules

    implicit none
    private

    !Public access procedures
    public :: initMohidLagrangian

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private source definitions parser routine. Builds the tracer sources from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_sources(parsedxml)
    implicit none
    type(Node), intent(in), pointer :: parsedxml    !>.xml file handle

    type(NodeList), pointer :: sourcedefList    !> Node list for simulationdefs
    type(NodeList), pointer :: sourceList       !> Node list for sources
    type(NodeList), pointer :: sourcedetailList
    type(Node), pointer :: sourcedef
    type(Node), pointer :: source               !> Single source block to process
    type(Node), pointer :: source_detail
    integer :: i, j, k
    !source vars
    integer :: id
    type(string) :: name, source_geometry
    character(80) :: name_char, source_geometry_char
    real(prec) :: emitting_rate
    class(*), allocatable :: geometry

    nullify(sourcedefList)
    nullify(sourceList)
    sourcedefList => getElementsByTagname(parsedxml, "sourcedef")       !searching for tags with the 'simulationdefs' name
    sourcedef => item(sourcedefList,0)
    sourceList => getElementsByTagname(sourcedef, "source")
    print*, "There are ", getLength(sourceList), " tracer sources"

    call allocSources(getLength(sourceList))                         !allocating the source objects

    do j = 0, getLength(sourceList) - 1
        source => item(sourceList,j)
        nullify(sourcedetailList)
        sourcedetailList => getElementsByTagname(source, "setsource")   !searching for tags with the 'setsource' name
        if (associated(sourcedetailList)) then
            source_detail => item(sourcedetailList, 0)
            call extractDataAttribute(source_detail, "id", id)
            call extractDataAttribute(source_detail, "name", name_char)
            name=trim(name_char)
        else
            print*, "Could not find any 'setsource' tag for source", j, "in input file, stoping."
            stop
        endif
        nullify(sourcedetailList)
        sourcedetailList => getElementsByTagname(source, "set")         !searching for tags with the 'set' name
        if (associated(sourcedetailList)) then
            source_detail => item(sourcedetailList, 0)
            call extractDataAttribute(source_detail, "emitting_rate", emitting_rate)
        else
            print*, "Could not find any 'set' tag for source", j, "in input file, stoping."
            stop
        endif
        !now we need to find out the geometry of the source and read accordingly
        nullify(sourcedetailList)
        source_detail => item(getChildNodes(source),5)                       !this is really ugly, but there seems to be a bug here. This effectivelly hardcodes the geometry position in the xml...
        source_geometry = getNodeName(source_detail)
        print*, "Source ",name ," is a ", source_geometry
        select case (name%chars())
        case ('point')
            allocate(type(point)::geometry)            
        case ('sphere')
            allocate(type(sphere)::geometry)
        case ('box')
            allocate(type(box)::geometry)
        case ('line')
            allocate(type(line)::geometry)
        case default
            print*, "Source", j, " geometry type is not supported, stopping"
            stop
        end select

    enddo





    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private simulation definitions parser routine. Builds the simulation geometric space from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_simdefs(parsedxml)
    implicit none
    type(Node), intent(in), pointer :: parsedxml    !>.xml file handle

    type(NodeList), pointer :: parameterList        !> Node list for parameters
    type(Node), pointer :: parmt                    !> Single parameter block to process
    real(prec) :: points(3), i
    type(string) :: pts(2)

    nullify(parameterList)
    parameterList => getElementsByTagname(parsedxml, "resolution")       !searching for tags with the 'simulationdefs' name
    if (associated(parameterList)) then
        parmt => item(parameterList, 0)
        call extractDataAttribute(parmt, "dp", points(1))
        call setSimDp(points(1))
    else
        print*, "Could not find any 'simulationdefs' tag in input file, stoping."
        stop
    endif

    pts=(/ 'pointmin', 'pointmax'/) !strings to search for
    do i=1, size(pts)
        nullify(parameterList)
        parameterList => getElementsByTagname(parsedxml, pts(i)%chars())       !searching for tags with the 'pts(i)' name
        if (associated(parameterList)) then
            parmt => item(parameterList, 0)
            call extractDataAttribute(parmt, "x", points(1))
            call extractDataAttribute(parmt, "y", points(2))
            call extractDataAttribute(parmt, "z", points(3))
            call setSimBounds(pts(i), points)
        else
            print*, "Could not find any ", pts(i)%chars(), " tag in input file, stoping"
            stop
        endif
    enddo

    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private case constant parser routine. Builds the simulation parametric space from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_caseconstants(parsedxml)
    implicit none
    type(Node), intent(in), pointer :: parsedxml    !>.xml file handle

    type(NodeList), pointer :: parameterList        !> Node list for parameters
    type(Node), pointer :: parmt                    !> Single parameter block to process
    real(prec) :: points(3)

    nullify(parameterList)
    parameterList => getElementsByTagname(parsedxml, "Gravity")       !searching for tags with the 'Gravity' name
    if (associated(parameterList)) then
        parmt => item(parameterList, 0)
        call extractDataAttribute(parmt, "x", points(1))
        call extractDataAttribute(parmt, "y", points(2))
        call extractDataAttribute(parmt, "z", points(3))
        call setSimGravity(points)
    else
        print*, "Could not find any 'Gravity' tag in input file, assuming (0,0,-9.81) m s-2."
        call setSimGravity((/0.0,0.0,-9.81/))
    endif
    nullify(parameterList)

    parameterList => getElementsByTagname(parsedxml, "Rho_ref")       !searching for tags with the 'Rho_ref' name
    if (associated(parameterList)) then
        parmt => item(parameterList, 0)
        call extractDataAttribute(parmt, "value", points(1))
        call setSimRho(points(1))
    else
        print*, "Could not find any 'Rho_ref' tag in input file, will assume 1000 kg m-3."
    endif

    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private parameter parser routine. Builds the simulation parametric space from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_parameters(parsedxml)
    implicit none
    type(Node), intent(in), pointer :: parsedxml    !>.xml file handle

    type(NodeList), pointer :: parameterList        !> Node list for parameters
    type(Node), pointer :: parmt                    !> Single parameter block to process
    integer :: i
    type(string) :: parmkey, parmvalue
    character(80) :: parmkey_char, parmvalue_char

    nullify(parameterList)
    parameterList => getElementsByTagname(parsedxml, "parameter")       !searching for tags with the 'parameter' name
    if (associated(parameterList)) then                                !checking if the list is not empty
        do i = 0, getLength(parameterList) - 1                          !extracting parameter tags one by one
            parmt => item(parameterList, i)
            call extractDataAttribute(parmt, "key", parmkey_char)       !name of the parameter
            call extractDataAttribute(parmt, "value", parmvalue_char)   !value of the parameter
            parmkey=trim(parmkey_char)
            parmvalue=trim(parmvalue_char)
            call setSimParameter(parmkey,parmvalue)
        enddo
    else
        print*, "Could not find any 'parameter' tag in input file, stopping."
        stop
    endif

    return
    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public xml parser routine. Builds the simulation space from the input xml case file.
    !
    !> @param[in] xmlfilename
    !---------------------------------------------------------------------------
    subroutine initMohidLagrangian(xmlfilename)
    implicit none
    type(string), intent(in) :: xmlfilename         !> .xml file name

    !local vars
    type(Node), pointer :: xmldoc                   !> .xml file handle
    integer :: i

    xmldoc => parseFile(xmlfilename%chars(), iostat=i)
    if (i/=0) then
        print*, "Could not open", xmlfilename%chars(), " input file, give me at least that!"
        stop
    endif

    call init_parameters(xmldoc)     !Reading the parameters
    call init_caseconstants(xmldoc)  !Reading the case constants
    call init_simdefs(xmldoc)
    call init_sources(xmldoc)

    call destroy(xmldoc)

    return

    end subroutine







    end module initialize
