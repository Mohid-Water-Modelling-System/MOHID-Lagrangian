
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

    parameterList => getElementsByTagname(parsedxml, "parameter")
    do i = 0, getLength(parameterList) - 1
        parmt => item(parameterList, i)
        call extractDataAttribute(parmt, "key", parmkey_char)
        call extractDataAttribute(parmt, "value", parmvalue_char)
        parmkey=trim(parmkey_char)
        parmvalue=trim(parmvalue_char)
        call setSimParameter(parmkey,parmvalue) 
    enddo

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

    call initParameters(xmldoc) !Reading the parameters


    call destroy(xmldoc)

    return

    end subroutine







    end module initialize
