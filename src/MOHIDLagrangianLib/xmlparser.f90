
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : xmlparser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module with xml parsing class and methods, encapsulates the
    !> FOX_dom library.
    !------------------------------------------------------------------------------

    module xmlParser_mod

    use FoX_common, only : rts, countrts
    use FoX_dom
    use penf
    use vecfor_r8p
    use stringifor

    use simulationPrecision_mod
    use simulationLogger_mod

    implicit none
    private

    type :: xmlparser_class  !< The .xml parser class
    contains
    procedure :: getFile
    procedure :: closeFile
    procedure :: getLeafAttribute
    procedure :: getNodeAttribute
    procedure :: getNodeVector
    procedure :: gotoNode
    procedure :: getPolygonFromKMZFile
    end type xmlparser_class

    !Simulation variables
    type(xmlparser_class) :: XMLReader

    !Public access vars
    public :: XMLReader, xmlparser_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a xml file and returns a pointer to the master node.
    !> @param[in] self, xmldoc, xmlfilename, mandatory
    !---------------------------------------------------------------------------
    subroutine getFile(self, xmldoc, xmlfilename, mandatory)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(out), pointer :: xmldoc   !< Node that contains the parsed file
    type(string), intent(in) :: xmlfilename      !< File name
    logical, intent(in), optional :: mandatory   !<Swich for optional or mandatory tags
    logical :: mand
    integer :: err
    type(string) :: outext

    if (present(mandatory)) then
        mand = mandatory
    else
        mand = .true.
    end if

    xmldoc => parseFile(xmlfilename%chars(), iostat=err) !using FOX function
    if (err==0) then
        outext='->Reading .xml file '//xmlfilename
        call Log%put(outext)
    else
        nullify(xmldoc)
        if (.not.mand) then
            outext='[XMLReader::getFile]: no '//xmlfilename//' file, or file is invalid. Ignoring'
            call Log%put(outext)
        else
            outext='[XMLReader::getFile]: no '//xmlfilename//' file, or file is invalid. Stoping'
            call Log%put(outext)
            stop
        end if
    end if
    end subroutine getFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that closes a parsed xml file or node.
    !> @param[in] self, xmldoc
    !---------------------------------------------------------------------------
    subroutine closeFile(self, xmldoc)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(out), pointer :: xmldoc   !< Node that conatins the parsed file
    call destroy(xmldoc) !using FOX function
    end subroutine closeFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses an xml attribute. Reads the requested attribute
    !> from a given leaf node,
    !> @param[in] self, xmlnode, att_name, att_value
    !---------------------------------------------------------------------------
    subroutine getLeafAttribute(self, xmlnode, att_name, att_value)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: att_name        !<Atribute name to collect from tag
    type(string), intent(inout) :: att_value      !<Attribute value
    character(CHAR_LEN) :: att_value_chars
    call extractDataAttribute(xmlnode, att_name%chars(), att_value_chars) !using FOX function
    att_value=trim(att_value_chars)
    end subroutine getLeafAttribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses an attribute from an xml node. In the format
    !> '<Tag att_name="att_value"/>'
    !> @param[in] self, xmlnode, tag, att_name, att_value, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine getNodeAttribute(self, xmlnode, tag, att_name, att_value, read_flag, mandatory)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: tag             !<Tag to search in xml node
    type(string), intent(in) :: att_name        !<Atribute name to collect from tag
    type(string), intent(inout) :: att_value      !<Attribute value
    logical, intent(out), optional :: read_flag !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory  !<Swich for optional or mandatory tags
    logical :: mand
    type(string) :: outext, nodename
    character(CHAR_LEN) :: att_value_chars
    type(NodeList), pointer :: target_node_list, nodeChildren
    type(Node), pointer :: nodedetail
    logical :: validtag
    integer :: i

    if (present(mandatory)) then
        mand = mandatory
    else
        mand = .true.
    end if

    validtag = .false.
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info) !using FOX function
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node !using FOX function
        nodename = getLocalName(nodedetail)  !finding its name !using FOX function
        if (nodename == tag) then
            validtag=.true.
            exit
        end if
    end do
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name !using FOX function
        nodedetail => item(target_node_list, 0) !using FOX function
        call extractDataAttribute(nodedetail, att_name%chars(), att_value_chars) !using FOX function
        att_value=trim(att_value_chars)
        if (present(read_flag)) then
            read_flag = .true.
            !if (att_value%to_number(kind=1._R8P) <= 1.0/100000.0) read_flag = .false.
        end if
    else
        if(.not.mand) then
            if (present(read_flag)) then
                read_flag =.false.
            end if
        else
            outext='Could not find any "'//tag//'" tag for xml node "'//getNodeName(xmlnode)//'", stoping'
            call Log%put(outext)
            stop
        end if
    end if
    end subroutine getNodeAttribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to parse xyz vectors in xml files.
    !> Vector must be in format '<Tag x="vec%x" y="vec%y" z="vec%z"/>'
    !> @param[in] self, xmlnode, tag, vec, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine getNodeVector(self, xmlnode, tag, vec, read_flag, mandatory)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: tag             !<Tag to search in xml node
    type(vector), intent(inout) :: vec            !<Vector to fill with read contents
    logical, intent(out), optional :: read_flag  !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory  !<Swich for optional or mandatory tags
    logical :: mand
    type(string) :: outext, nodename
    type(NodeList), pointer :: target_node_list, nodeChildren
    type(Node), pointer :: nodedetail
    logical :: validtag
    integer :: i

    if (present(mandatory)) then
        mand = mandatory
    else
        mand = .true.
    end if

    vec%x=MV !marking the array as not read
    validtag = .false.
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info) !using FOX function
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node !using FOX function
        nodename = getLocalName(nodedetail)  !finding its name !using FOX function
        if (nodename == tag) then
            validtag =.true.
            exit
        end if
    end do
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name !using FOX function
        nodedetail => item(target_node_list, 0) !using FOX function
        call extractDataAttribute(nodedetail, "x", vec%x) !using FOX function
        call extractDataAttribute(nodedetail, "y", vec%y)
        call extractDataAttribute(nodedetail, "z", vec%z)
        if (present(read_flag)) then
            read_flag =.true.
            if (vec%normL2() <= 1.0/100000.0) read_flag = .false.
        end if
    else
        if(.not.mand) then
            if (present(read_flag)) then
                read_flag =.false.
            end if
        else
            outext='Could not find any "'//tag//'" tag for xml node "'//getNodeName(xmlnode)//'", stoping'
            call Log%put(outext)
            stop
        end if
    end if
    end subroutine getNodeVector

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that retrieves a node from within a node.
    !> Returns a nullifyed pointer if not found, stops if mandatory.
    !> @param[in] self, currentNode, targetNode, targetNodeName, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine gotoNode(self, currentNode, targetNode, targetNodeName, read_flag, mandatory)
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: currentNode
    type(Node), intent(out), pointer :: targetNode
    type(string), intent(in) :: targetNodeName
    logical, intent(out), optional :: read_flag  !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory   !<Swich for optional or mandatory tags
    logical :: mand
    type(NodeList), pointer :: target_node_list
    type(string) :: outext, nodename
    integer :: i
    logical :: target_node_exists

    if (present(mandatory)) then
        mand = mandatory
    else
        mand = .true.
    end if

    target_node_exists = .false.
    target_node_list => getChildNodes(currentNode) !using FOX function
    do i=0, getLength(target_node_list)-1
        targetNode => item(target_node_list,i) !grabing a node !using FOX function
        nodename = getLocalName(targetNode)  !finding its name !using FOX function
        if (nodename == targetNodeName) then !found our target node
            target_node_exists = .true.
            if (present(read_flag)) then
                read_flag =.true.
            end if
            exit
        end if
    end do
    if (.not.target_node_exists) then
        nullify(targetNode)
        if(.not.mand) then
            !outext='Could not find any node called "'//targetNodeName//'" in the xml file, ignoring'
            !call Log%put(outext)
            if (present(read_flag)) then
                read_flag =.false.
            end if
        else
            outext='Could not find any node called "'//targetNodeName//'" in the xml file, stoping'
            call Log%put(outext)
            stop
        end if
    end if
    end subroutine gotoNode

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that retrieves a node from within a node.
    !> Returns a nullifyed pointer if not found, stops if mandatory.
    !> @param[in] self, filename, vertex
    !---------------------------------------------------------------------------
    subroutine getPolygonFromKMZFile(self, filename, vertex)
    class(xmlparser_class), intent(in) :: self
    type(string), intent(in) :: filename
    type(vector), dimension(:), allocatable, intent(out) :: vertex
    type(Node), pointer :: xmlDoc, xmlNode
    type(string) :: outext, tag
    real(prec), dimension(:), allocatable :: temp
    integer :: i

    call self%getFile(xmlDoc,fileName)
    xmlNode => item(getElementsByTagName(xmlDoc, "coordinates"), 0)
    allocate(temp(countrts(getTextContent(xmlNode), 0.0d0)))
    call extractDataContent(xmlNode,temp)
    if (mod(size(temp),3) /= 0) then
        outext='[xmlParser::getPolygonFromKMZFile]: 3D Polygon is not well defined, stoping'
        call Log%put(outext)
        stop
    else
        allocate(vertex(size(temp)/3))
        do i=1, size(vertex)
            vertex(i) = temp(1+(i-1)*3)*ex + temp(2+(i-1)*3)*ey + temp(3+(i-1)*3)*ez
        end do
    end if    

    end subroutine getPolygonFromKMZFile

    end module xmlParser_mod
