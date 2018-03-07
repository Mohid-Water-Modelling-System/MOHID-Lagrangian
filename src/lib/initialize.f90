
    module initialize

    use tracer3D
    use simulation_parameters

    use FoX_dom
    use commom_modules

    implicit none
    private

    !Public access procedures
    public :: init_mohidlagrangian

    contains
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private simulation definitions parser routine. Builds the simulation geometric space from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine initSimDefs(parsedxml)
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
    subroutine initCaseConstants(parsedxml)
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
    subroutine initParameters(parsedxml)
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
    subroutine init_mohidlagrangian(xmlfilename)
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

    call initParameters(xmldoc)     !Reading the parameters
    call initCaseConstants(xmldoc)  !Reading the case constants
    call initSimDefs(xmldoc)
    
    
    call destroy(xmldoc)

    return

    end subroutine







    end module initialize
