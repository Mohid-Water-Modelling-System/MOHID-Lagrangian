
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
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module with the simulation initialization related definitions and methods.
    !> Has one public access routine that is incharge of building the simulation
    !> space from input files.
    !------------------------------------------------------------------------------

    module simulationInitialize_mod

    use common_modules
    use tracerBase_mod
    use xmlParser_mod
    use sources_mod

    use FoX_dom

    implicit none
    private

    !Public access procedures
    public :: InitFromXml

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private property xml parser routine. Reads the properties tab from the xml
    !> file and links these to the corresponding Source
    !> @param[in] linksNode
    !---------------------------------------------------------------------------
    subroutine linkPropertySources(linksNode)
    implicit none
    type(Node), intent(in), pointer :: linksNode

    type(NodeList), pointer :: linkList
    type(Node), pointer :: anode
    type(Node), pointer :: xmlProps                   !< .xml file handle
    type(string) :: xmlfilename, outext
    integer :: i, p
    type(string) :: att_name, att_val, tag
    type(string) :: sourceid, sourcetype, sourceprop

    linkList => getElementsByTagname(linksNode, "link")
    do i = 0, getLength(linkList) - 1
        anode => item(linkList,i)
        att_name="source"
        call XMLReader%getLeafAttribute(anode,att_name,sourceid)
        att_name="type"
        call XMLReader%getLeafAttribute(anode,att_name,sourcetype)
        att_name="property"
        call XMLReader%getLeafAttribute(anode,att_name,sourceprop)
        !find the source and save the type and property name
        call tempSources%setPropertyNames(sourceid,sourcetype,sourceprop)
    enddo

    !parse the properties file
    xmlfilename = Globals%Names%propsxmlfilename
    call XMLReader%getFile(xmlProps,xmlfilename)

    !Go to the materials node
    tag = "materials"
    call XMLReader%gotoNode(xmlProps,xmlProps,tag)

    !find and set the actual atributes of the properties
    att_name="value"
    do i = 1, size(tempSources%src)
        tag = tempSources%src(i)%prop%property_type
        if (tag .ne. 'base') then
            call XMLReader%gotoNode(xmlProps,anode,tag) !finding the material type node
            tag = tempSources%src(i)%prop%property_name
            call XMLReader%gotoNode(anode,anode,tag)     !finding the actual material node
            do p = 1, size(Globals%SrcProp%baselist)
                call XMLReader%getNodeAttribute(anode, Globals%SrcProp%baselist(p), att_name, att_val)
                call tempSources%src(i)%setPropertyAtributes(Globals%SrcProp%baselist(p), att_val)
            end do
            if (tempSources%src(i)%isParticulate()) then
                do p = 1, size(Globals%SrcProp%particulatelist)
                    call XMLReader%getNodeAttribute(anode, Globals%SrcProp%particulatelist(p), att_name, att_val)
                    call tempSources%src(i)%setPropertyAtributes(Globals%SrcProp%particulatelist(p), att_val)
                end do
            end if
            !Run integrety check on the properties to see if Source is well defined
            call tempSources%src(i)%check()
        end if
    end do
    outext='-->Sources properties are set'
    call Log%put(outext,.false.)

    call XMLReader%closeFile(xmlProps)

    end subroutine linkPropertySources


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private property xml parser routine. Reads the properties tab from the xml
    !> file and links these to the corresponding source
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_properties(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: props_node          !< Single properties block to process
    type(string) :: outext
    type(string) :: tag, att_name

    tag="properties"    !the node we want
    call XMLReader%gotoNode(case_node,props_node,tag,mandatory =.false.)
    if (associated(props_node)) then
        tag="propertyfile"
        att_name="name"
        call XMLReader%getNodeAttribute(props_node, tag, att_name, Globals%Names%propsxmlfilename) !getting the file name from that tag
        outext='-->Properties to link to Sources found at '//Globals%Names%propsxmlfilename
        call Log%put(outext,.false.)
        tag="links"
        call XMLReader%gotoNode(props_node,props_node,tag) !getting the links node
        call linkPropertySources(props_node)          !calling the property linker
    else
        outext='-->No properties to link to Sources, assuming pure Lagrangian tracers'
        call Log%put(outext,.false.)
    endif

    end subroutine init_properties
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> naming xml parser routine. Reads the naming file(s), opens the file(s) and
    !> stores the naming conventions for input files
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_naming(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: naming_node          !< Single naming block to process
    type(Node), pointer :: temp
    type(NodeList), pointer :: namingfileList
    type(string) :: outext
    type(string) :: tag, att_name
    type(string), allocatable, dimension(:) :: namingFilename
    integer :: i

    tag="naming"    !the node we want
    call XMLReader%gotoNode(case_node,naming_node,tag,mandatory =.false.)
    if (associated(naming_node)) then
        namingfileList => getElementsByTagname(naming_node, "namingfile")       !searching for tags with the 'namingfile' name
        allocate(namingFilename(getLength(namingfileList)))
        do i = 0, getLength(namingfileList) - 1
            temp => item(namingfileList, i)
            att_name="name"
            call XMLReader%getLeafAttribute(temp,att_name,namingFilename(i+1))
        end do
        call Globals%setNamingConventions(namingFilename)
    else
        outext='-->No naming files, assuming basic naming settings for variables from input files'
        call Log%put(outext,.false.)
    endif

    end subroutine init_naming


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private geometry xml parser routine. Reads a geometry from the xml depending on the geometry type of the node
    !> @param[in] source, source_detail, source_shape
    !---------------------------------------------------------------------------
    subroutine read_xml_geometry(source,source_detail,source_shape)
    implicit none
    type(Node), intent(in), pointer :: source           !<Working xml node
    type(Node), intent(in), pointer :: source_detail    !<Working xml node details
    class(shape), intent(inout) :: source_shape         !<Geometrical object to fill
    type(string) :: outext
    type(string) :: tag
    select type (source_shape)
    type is (shape)
    class is (box)
        tag='point'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%pt)
        tag='size'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%size)
    class is (point)
        tag='point'
        call XMLReader%getNodeVector(source,tag,source_shape%pt)
    class is (line)
        tag='pointa'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%pt)
        tag='pointb'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%last)
    class is (sphere)
        tag='point'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%pt)
        call extractDataAttribute(source_detail, "radius", source_shape%radius)
        class default
        outext='[read_xml_geometry]: unexpected type for geometry object!'
        call Log%put(outext)
        stop
    end select
    end subroutine read_xml_geometry

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private source definitions parser routine. Builds the tracer sources from the input xml case file.
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_sources(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node
    type(string) :: outext
    type(NodeList), pointer :: sourceList           !< Node list for sources
    type(NodeList), pointer :: sourceChildren       !< Node list for source node children nodes
    type(Node), pointer :: sourcedef
    type(Node), pointer :: source_node              !< Single source block to process
    type(Node), pointer :: source_detail
    integer :: i, j
    logical :: readflag
    integer :: id
    type(string) :: name, source_geometry, tag, att_name, att_val, rate_file
    real(prec) :: emitting_rate, start, finish
    logical :: emitting_fixed
    class(shape), allocatable :: source_shape

    rate_file = 'not_set'
    readflag = .false.
    outext='-->Reading case Sources'
    call Log%put(outext,.false.)
    tag="sourcedef"    !the node we want
    call XMLReader%gotoNode(case_node,sourcedef,tag)
    sourceList => getElementsByTagname(sourcedef, "source")
    !allocating the temporary source objects
    call tempSources%initialize(getLength(sourceList))

    do j = 0, getLength(sourceList) - 1
        source_node => item(sourceList,j)
        tag="setsource"
        att_name="id"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val)
        id=att_val%to_number(kind=1_I1P)
        att_name="name"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, name)
        !reading emission rate, need to check for options
        tag="rate_dt"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            emitting_rate = 1.0/(att_val%to_number(kind=1._R4P)*Globals%SimDefs%dt)
            emitting_fixed = .true.
        end if
        tag="rate"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            emitting_rate = att_val%to_number(kind=1._R4P)
            emitting_fixed = .true.
        end if
        tag="rate_file"
        att_name="name"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            rate_file = att_val
            emitting_fixed = .false.
        end if
        tag="active"
        att_name="start"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            start = att_val%to_number(kind=1._R4P)
        else
            start = 0.0
        end if
        att_name="end"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val,readflag,.false.)
        if (readflag.and.att_val%is_number()) then
            finish = att_val%to_number(kind=1._R4P)
        else
            finish = Globals%Parameters%TimeMax
        end if
        !now we need to find out the geometry of the source and read accordingly
        sourceChildren => getChildNodes(source_node) !getting all of the nodes bellow the main source node (all of it's private info)
        do i=0, getLength(sourceChildren)-1
            source_detail => item(sourceChildren,i) !grabing a node
            source_geometry = getLocalName(source_detail)  !finding its name
            if (Geometry%inlist(source_geometry)) then  !if the node is a valid geometry name
                call Geometry%allocateShape(source_geometry,source_shape)                
                call read_xml_geometry(source_node,source_detail,source_shape)
                exit
            end if
        end do
        !initializing Source j
        call tempSources%src(j+1)%initialize(id,name,emitting_rate,emitting_fixed,rate_file,start,finish,source_geometry,source_shape)

        deallocate(source_shape)
    enddo

    end subroutine init_sources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private simulation definitions parser routine. Builds the simulation geometric space from the input xml case file.
    !> @param[in] case_node
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
    call Log%put(outext,.false.)

    tag="simulationdefs"    !the node we want
    call XMLReader%gotoNode(case_node,simdefs_node,tag)
    tag="resolution"
    att_name="dp"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val)
    call Globals%SimDefs%setdp(att_val)
    tag="timestep"
    att_name="dt"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val)
    call Globals%SimDefs%setdt(att_val)
    pts=(/ 'pointmin', 'pointmax'/) !strings to search for
    do i=1, size(pts)
        call XMLReader%getNodeVector(simdefs_node, pts(i), coords)
        call Globals%SimDefs%setboundingbox(pts(i), coords)
    enddo
    call Globals%SimDefs%print()

    end subroutine init_simdefs

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private case constant parser routine. Builds the simulation parametric space from the input xml case file.
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_caseconstants(case_node)
    implicit none
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: constants_node       !< Single constants block to process
    type(string) :: outext
    type(string) :: tag, att_name, att_val
    type(vector) :: coords
    logical :: readflag

    outext='-->Reading case constants'
    call Log%put(outext,.false.)

    tag="constantsdef"    !the node we want
    call XMLReader%gotoNode(case_node,constants_node,tag,readflag,.false.)
    if (readflag) then !if the node exists, since his one is not mandatory
        tag="Gravity"
        call XMLReader%getNodeVector(constants_node,tag,coords,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setgravity(coords)
        endif
        tag="Z0"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setz0(att_val)
        endif
        tag="Rho_ref"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setrho(att_val)
        endif
    endif
    call Globals%Constants%print()

    end subroutine init_caseconstants

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private parameter parser routine. Builds the simulation parametric space from the input xml case file.
    !> @param[in] execution_node
    !---------------------------------------------------------------------------
    subroutine init_parameters(execution_node)
    implicit none
    type(Node), intent(in), pointer :: execution_node

    type(string) :: outext
    type(NodeList), pointer :: parameterList        !< Node list for parameters
    type(Node), pointer :: parmt, parameters_node
    integer :: i
    type(string) :: parmkey, parmvalue, tag, att_name

    outext='-->Reading case parameters'
    call Log%put(outext,.false.)

    tag="parameters"    !the node we want
    call XMLReader%gotoNode(execution_node,parameters_node,tag)
    parameterList => getElementsByTagname(parameters_node, "parameter")       !searching for tags with the 'parameter' name
    do i = 0, getLength(parameterList) - 1                          !extracting parameter tags one by one
        parmt => item(parameterList, i)
        att_name="key"
        call XMLReader%getLeafAttribute(parmt,att_name,parmkey)
        att_name="value"
        call XMLReader%getLeafAttribute(parmt,att_name,parmvalue)
        call Globals%Parameters%setParam(parmkey,parmvalue)
    end do
    call Globals%Parameters%check()
    call Globals%Parameters%print()
    
    call Globals%setTimeDate()
    call Globals%SimTime%print()

    end subroutine init_parameters

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public xml parser routine. Builds the simulation space from the input xml case file.
    !> @param[in] xmlfilename
    !---------------------------------------------------------------------------
    subroutine InitFromXml(xmlfilename)
    implicit none
    type(string), intent(in) :: xmlfilename         !< .xml file name
    type(string) :: outext, tag
    type(Node), pointer :: xmldoc                   !< .xml file handle
    type(Node), pointer :: case_node
    type(Node), pointer :: execution_node

    call XMLReader%getFile(xmldoc,xmlfilename)
    Globals%Names%mainxmlfilename = xmlfilename
    Globals%Names%casename = xmlfilename%basename(extension='.xml')
    outext='->Case name is '//Globals%Names%casename
    call Log%put(outext)

    tag="case"          !base document node
    call XMLReader%gotoNode(xmldoc,execution_node,tag)
    tag="execution"     !finding execution node
    call XMLReader%gotoNode(execution_node,execution_node,tag)
    tag="case"          !base document node
    call XMLReader%gotoNode(xmldoc,case_node,tag)
    tag="casedef"     !finding execution node
    call XMLReader%gotoNode(case_node,case_node,tag)

    ! building the simulation basic structures according to the case definition file
    ! every other structure in the simulation is built from these, i.e., not defined by the user directly
    call init_parameters(execution_node)
    call init_caseconstants(case_node)
    call init_simdefs(case_node)
    call init_sources(case_node)
    call init_properties(case_node)
    call init_naming(case_node)

    !setting the number of blocks to the correct ammount of selected threads
    Globals%SimDefs%numblocks = Globals%Parameters%numOPMthreads

    call XMLReader%closeFile(xmldoc)

    end subroutine InitFromXml

    end module simulationInitialize_mod
