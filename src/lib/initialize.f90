
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : initialize
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module with the simulation initialization related definitions and methods. Has one public access routine that is incharge of building the simulation space from input files.
    !------------------------------------------------------------------------------

    module initialize

    use commom_modules
    use tracer_base
    use simulation_xmlparser
    use source_identity
    use source_emitter
    use about

    use FoX_dom

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
    !> Private property xml parser routine. Reads the properties tab from the xml
    !> file and links these to the corresponding source
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine linkPropertySources(linksNode)
    implicit none
    type(Node), intent(in), pointer :: linksNode

    type(NodeList), pointer :: linkList
    type(Node), pointer :: linknode
    integer :: i
    character(80) :: sourceid_char, sourcetype_char, sourceprop_char
    type(string) :: sourceid, sourcetype, sourceprop

    linkList => getElementsByTagname(linksNode, "link")
    do i = 0, getLength(linkList) - 1
        linknode => item(linkList,i)
        call extractDataAttribute(linknode, "source", sourceid_char)
        call extractDataAttribute(linknode, "type", sourcetype_char)
        call extractDataAttribute(linknode, "property", sourceprop_char)
        sourceid=trim(sourceid_char)
        sourcetype=trim(sourcetype_char)
        sourceprop=trim(sourceprop_char)
        call setSourceProperties(sourceid,sourcetype,sourceprop)
    enddo

    return
    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private property xml parser routine. Reads the properties tab from the xml
    !> file and links these to the corresponding source
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_properties(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: props_node          !< Single properties block to process
    type(string) :: outext
    type(string) :: tag, att_name

    tag="properties"    !the node we want
    call gotoChildNode(case_node,props_node,tag)
    if (associated(props_node)) then
        tag="propertyfile"
        att_name="name"
        call readxmlatt(props_node, tag, att_name, Globals%FileNames%propsxmlfilename)  !getting the file name from that tag
        outext='-->Properties to link to Sources found at '//Globals%FileNames%propsxmlfilename
        call ToLog(outext)
        tag="links"
        call gotoChildNode(props_node,props_node,tag,.true.) !getting the links node
        call linkPropertySources(props_node)                 !calling the property linker
    else
        outext='-->No properties to link to Sources, assuming pure Lagrangian tracers'
        call ToLog(outext)
    endif

    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private geometry xml parser routine. Reads a geometry from the xml depending on the geometry type of the node
    !
    !> @param[in] source, geometry
    !---------------------------------------------------------------------------
    subroutine read_xml_geometry(source,source_detail,geometry)
    implicit none
    type(Node), intent(in), pointer :: source           !<Working xml node
    type(Node), intent(in), pointer :: source_detail    !<Working xml node details
    class(shape), intent(inout) :: geometry             !<Geometrical object to fill
    type(string) :: outext
    type(string) :: tag

    select type (geometry)
    type is (shape)
        !nothing to do
    class is (box)
        tag='point'
        call readxmlvector(source_detail,tag,geometry%pt)
        tag='size'
        call readxmlvector(source_detail,tag,geometry%size)
    class is (point)
        tag='point'
        call readxmlvector(source,tag,geometry%pt)
    class is (line)
        tag='pointa'
        call readxmlvector(source_detail,tag,geometry%pt)
        tag='pointb'
        call readxmlvector(source_detail,tag,geometry%last)
    class is (sphere)
        tag='point'
        call readxmlvector(source_detail,tag,geometry%pt)
        call extractDataAttribute(source_detail, "radius", geometry%radius)
        class default
        outext='read_xml_geometry: unexpected type for geometry object!'
        call ToLog(outext)
        stop
    end select

    end subroutine


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private source definitions parser routine. Builds the tracer sources from the input xml case file.
    !
    !> @param[in] parsedxml
    !---------------------------------------------------------------------------
    subroutine init_sources(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(string) :: outext
    type(NodeList), pointer :: sourcedefList        !< Node list for simulationdefs
    type(NodeList), pointer :: sourceList           !< Node list for sources
    type(NodeList), pointer :: sourceChildren      !< Node list for source node children nodes
    type(Node), pointer :: sourcedef
    type(Node), pointer :: source_node                   !< Single source block to process
    type(Node), pointer :: source_detail
    integer :: i, j
    !source vars
    integer :: id
    type(string) :: name, source_geometry, tag, att_name, att_val
    character(80) :: val, name_char, source_geometry_char
    real(prec) :: emitting_rate, start, finish
    class(shape), allocatable :: geometry

    outext='-->Reading case Sources'
    call ToLog(outext,.false.)

    tag="sourcedef"    !the node we want
    call gotoChildNode(case_node,sourcedef,tag,.true.)
    sourceList => getElementsByTagname(sourcedef, "source")

    call allocSources(getLength(sourceList))                         !allocating the source objects

    do j = 0, getLength(sourceList) - 1
        source_node => item(sourceList,j)
        tag="setsource"
        att_name="id"
        call ReadXMLatt(source_node, tag, att_name, att_val)
        id=att_val%to_number(kind=1_I1P)
        att_name="name"
        call ReadXMLatt(source_node, tag, att_name, name)
        tag="set"
        att_name="emitting_rate"
        call ReadXMLatt(source_node, tag, att_name, att_val)
        emitting_rate = att_val%to_number(kind=1._R4P)
        tag="active"
        att_name='start'
        call ReadXMLatt(source_node, tag, att_name, att_val,.false.)
        if (att_val%is_number()) then
            start = att_val%to_number(kind=1._R4P)
        else
            start = 0.0
        endif
        att_name='end'
        call ReadXMLatt(source_node, tag, att_name, att_val,.false.)
        if (att_val%is_number()) then
            finish = att_val%to_number(kind=1._R4P)
        else
            finish = Globals%Parameters%TimeMax
        endif
        !now we need to find out the geometry of the source and read accordingly
        sourceChildren => getChildNodes(source_node) !getting all of the nodes bellow the main source node (all of it's private info)
        do i=0, getLength(sourceChildren)-1
            source_detail => item(sourceChildren,i) !grabing a node
            source_geometry = getLocalName(source_detail)  !finding its name
            if (IsValidGeom(source_geometry)) then  !if the node is a valid geometry name
                select case (source_geometry%chars())
                case ('point')
                    allocate(point::geometry)
                case ('sphere')
                    allocate(sphere::geometry)
                case ('box')
                    allocate(box::geometry)
                case ('line')
                    allocate(line::geometry)
                    case default
                    outext='init_sources: unexpected type for geometry object!'
                    call ToLog(outext)
                    stop
                end select
                call read_xml_geometry(source_node,source_detail,geometry)
                exit
            endif
        enddo

        !initializing Source j
        call Source(j+1)%initialize(id,name,emitting_rate,start,finish,source_geometry,geometry)

        deallocate(geometry)
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
    subroutine init_simdefs(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(NodeList), pointer :: defsList       !< Node list for simdefs
    type(Node), pointer :: simdefs_node       !< Single simdefs block to process
    type(string) :: outext
    integer :: i
    type(string) :: pts(2), tag, att_name, att_val
    type(vector) :: coords

    outext='-->Reading case simulation definitions'
    call ToLog(outext,.false.)

    tag="simulationdefs"    !the node we want
    call gotoChildNode(case_node,simdefs_node,tag,.true.)
    tag="resolution"
    att_name="dp"
    call readxmlatt(simdefs_node, tag, att_name, att_val)
    call Globals%SimDefs%setdp(att_val)
    tag="timestep"
    att_name="dt"
    call readxmlatt(simdefs_node, tag, att_name, att_val)
    call Globals%SimDefs%setdt(att_val)
    pts=(/ 'pointmin', 'pointmax'/) !strings to search for
    do i=1, size(pts)
        call readxmlvector(simdefs_node, pts(i), coords)
        call Globals%SimDefs%setboundingbox(pts(i), coords)
    enddo
    call Globals%SimDefs%printout()

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
    subroutine init_caseconstants(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: constants_node       !< Single constants block to process
    type(string) :: outext
    type(string) :: tag, att_name, att_val
    type(vector) :: coords

    outext='-->Reading case constants'
    call ToLog(outext,.false.)

    tag="constantsdef"    !the node we want
    call gotoChildNode(case_node,constants_node,tag,.true.)
    tag="Gravity"
    call readxmlvector(constants_node, tag, coords,.false.)
    call Globals%Constants%setgravity(coords)
    tag="Rho_ref"
    att_name="value"
    call readxmlatt(constants_node, tag, att_name, att_val)
    call Globals%Constants%setrho(att_val)

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
    subroutine init_parameters(execution_node)
    implicit none
    type(Node), intent(in), pointer :: execution_node

    type(string) :: outext
    type(NodeList), pointer :: parameterList        !< Node list for parameters
    type(Node), pointer :: parmt, parameters_node
    integer :: i
    type(string) :: parmkey, parmvalue, tag
    character(80) :: parmkey_char, parmvalue_char

    outext='-->Reading case parameters'
    call ToLog(outext,.false.)

    tag="parameters"    !the node we want
    call gotoChildNode(execution_node,parameters_node,tag,.true.)
    parameterList => getElementsByTagname(parameters_node, "parameter")       !searching for tags with the 'parameter' name
    do i = 0, getLength(parameterList) - 1                          !extracting parameter tags one by one
        parmt => item(parameterList, i)
        call extractDataAttribute(parmt, "key", parmkey_char)       !name of the parameter
        call extractDataAttribute(parmt, "value", parmvalue_char)   !value of the parameter
        parmkey=trim(parmkey_char)
        parmvalue=trim(parmvalue_char)
        call Globals%Parameters%setparameter(parmkey,parmvalue)
    enddo
    call Globals%Parameters%check()
    call Globals%Parameters%printout()

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
    type(string), intent(in) :: xmlfilename         !< .xml file name
    type(string) :: outext, tag
    type(Node), pointer :: xmldoc                   !< .xml file handle
    type(Node), pointer :: case_node
    type(Node), pointer :: execution_node
    integer :: i

    call PrintLicPreamble
    call AllocateGeomList
    Globals%SimTime = 0.

    !check if log file was opened - if not stop here
    if (Log_unit==0) then
        outext='->Logger initialized'
        call ToLog(outext)
    else
        stop 'Logger has not been initialized, stopping'
    end if

    xmldoc => parseFile(xmlfilename%chars(), iostat=i)
    if (i==0) then
        outext='->Reading case definition from '//xmlfilename
        call ToLog(outext)
        Globals%FileNames%mainxmlfilename = xmlfilename
    else
        outext='Could not open '//xmlfilename//' input file, give me at least that!'
        call ToLog(outext)
        stop
    endif

    tag="case"          !base document node
    call gotoChildNode(xmldoc,execution_node,tag,.true.)
    tag="execution"     !finding execution node
    call gotoChildNode(execution_node,execution_node,tag,.true.)
    tag="case"          !base document node
    call gotoChildNode(xmldoc,case_node,tag,.true.)
    tag="casedef"     !finding execution node
    call gotoChildNode(case_node,case_node,tag,.true.)

    !initializing memory log
    call SimMemory%initialize()

    ! building the simulation basic structures according to the case definition file
    ! every other structure in the simulation is built from these, i.e., not defined by the user directly
    call init_parameters(execution_node)
    call init_caseconstants(case_node)
    call init_simdefs(case_node)
    call init_sources(case_node)
    call init_properties(case_node)

    !With the Sources initialized, now we initialize the Emmiter class, that automatically
    !allocates and initializes all of the useable tracers
    call Emitter%initialize(Source)

    !printing memory occupation at the time
    call SimMemory%detailedprintout()

    call destroy(xmldoc)

    return
    end subroutine

    end module initialize
