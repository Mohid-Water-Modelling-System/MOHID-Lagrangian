
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_xmlparser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module with the simulation xml parsing class and methods, Encapsulates the
    !> FOX_dom library.
    !------------------------------------------------------------------------------

    module simulation_xmlparser_mod

    use FoX_dom
    use common_modules

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
    end type xmlparser_class

    !the object tho expose to the rest of the code
    type(xmlparser_class) :: XMLReader

    !Public access vars
    public :: XMLReader

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses an xml file and returns a pointer to the master node.
    !> @param[in] xmldoc, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, xmldoc, xmlfilename)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(out), pointer :: xmldoc   !< Node that conatins the parsed file
    type(string), intent(in) :: xmlfilename      !< File name
    integer :: err
    type(string) :: outext
    xmldoc => parseFile(xmlfilename%chars(), iostat=err) !using FOX function
    if (err==0) then
        outext='->Reading .xml file '//xmlfilename
        call Log%put(outext)
    else
        outext='[XMLReader::getFile]: no '//xmlfilename//' file, or file is invalid. Stoping'
        call Log%put(outext)
        stop
    endif
    end subroutine getFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that closes a parsed xml file or node.
    !> @param[in] xmldoc
    !---------------------------------------------------------------------------
    subroutine closeFile(self, xmldoc)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(out), pointer :: xmldoc   !< Node that conatins the parsed file
    call destroy(xmldoc) !using FOX function
    end subroutine closeFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses an xml attribute. Reads the requested attribute
    !> from a given leaf node,
    !> @param[in] xmlnode, att_name, att_value, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine getLeafAttribute(self, xmlnode, att_name, att_value)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: att_name        !<Atribute name to collect from tag
    type(string), intent(out) :: att_value      !<Attribute value
    character(80) :: att_value_chars
    call extractDataAttribute(xmlnode, att_name%chars(), att_value_chars) !using FOX function
    att_value=trim(att_value_chars)
    end subroutine getLeafAttribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses an attribute from an xml node. In the format
    !> <Tag att_name="att_value"/>
    !> @param[in] xmlnode, tag, att_name, att_value, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine getNodeAttribute(self, xmlnode, tag, att_name, att_value, read_flag, mandatory)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: tag             !<Tag to search in xml node
    type(string), intent(in) :: att_name        !<Atribute name to collect from tag
    type(string), intent(out) :: att_value      !<Attribute value
    logical, intent(out), optional :: read_flag  !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory  !<Swich for optional or mandatory tags

    type(string) :: outext, nodename
    character(80) :: att_value_chars
    type(NodeList), pointer :: target_node_list, nodeChildren
    type(Node), pointer :: nodedetail
    logical :: validtag
    integer :: i

    validtag = .false.
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info) !using FOX function
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node !using FOX function
        nodename = getLocalName(nodedetail)  !finding its name !using FOX function
        if (nodename == tag) then
            validtag=.true.
            exit
        endif
    enddo
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name !using FOX function
        nodedetail => item(target_node_list, 0) !using FOX function
        call extractDataAttribute(nodedetail, att_name%chars(), att_value_chars) !using FOX function
        att_value=trim(att_value_chars)
        if (present(read_flag)) then
            read_flag =.true.
        endif
    else
        if(present(mandatory)) then
            if(mandatory.eqv..false.) then
                if (present(read_flag)) then
                    read_flag =.false.
                endif
            endif
        else
            outext='Could not find any "'//tag//'" tag for xml node "'//getNodeName(xmlnode)//'", stoping'
            call Log%put(outext)
            stop
        endif
    endif
    end subroutine getNodeAttribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to parse xyz vectors in xml files.
    !> Vector must be in format <Tag x="vec%x" y="vec%y" z="vec%z"/>
    !> @param[in] xmlnode, tag, vec, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine getNodeVector(self, xmlnode, tag, vec, read_flag, mandatory)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: xmlnode  !<Working xml node
    type(string), intent(in) :: tag             !<Tag to search in xml node
    type(vector), intent(out) :: vec            !<Vector to fill with read contents
    logical, intent(out), optional :: read_flag  !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory  !<Swich for optional or mandatory tags

    type(string) :: outext, nodename
    type(NodeList), pointer :: target_node_list, nodeChildren
    type(Node), pointer :: nodedetail
    logical :: validtag
    integer :: i

    vec%x=MV !marking the array as not read
    validtag = .false.
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info) !using FOX function
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node !using FOX function
        nodename = getLocalName(nodedetail)  !finding its name !using FOX function
        if (nodename == tag) then
            validtag =.true.
            exit
        endif
    enddo
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name !using FOX function
        nodedetail => item(target_node_list, 0) !using FOX function
        call extractDataAttribute(nodedetail, "x", vec%x) !using FOX function
        call extractDataAttribute(nodedetail, "y", vec%y)
        call extractDataAttribute(nodedetail, "z", vec%z)
        if (present(read_flag)) then
            read_flag =.true.
        endif
    else
        if(present(mandatory)) then
            if(mandatory.eqv..false.) then
                if (present(read_flag)) then
                    read_flag =.false.
                endif
            endif
        else
            outext='Could not find any "'//tag//'" tag for xml node "'//getNodeName(xmlnode)//'", stoping'
            call Log%put(outext)
            stop
        endif
    endif
    end subroutine getNodeVector

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that retrieves a node from within a node.
    !> Returns a nullifyed pointer if not found, stops if mandatory.
    !> @param[in] currentNode, targetNode, targetNodeName, read_flag, mandatory
    !---------------------------------------------------------------------------
    subroutine gotoNode(self, currentNode, targetNode, targetNodeName, read_flag, mandatory)
    implicit none
    class(xmlparser_class), intent(in) :: self
    type(Node), intent(in), pointer :: currentNode
    type(Node), intent(out), pointer :: targetNode
    type(string), intent(in) :: targetNodeName
    logical, intent(out), optional :: read_flag  !< Optional flag to capture read/non-read status
    logical, intent(in), optional :: mandatory   !<Swich for optional or mandatory tags

    type(NodeList), pointer :: target_node_list
    type(string) :: outext, nodename
    integer :: i
    logical :: target_node_exists

    target_node_exists = .false.
    target_node_list => getChildNodes(currentNode) !using FOX function
    do i=0, getLength(target_node_list)-1
        targetNode => item(target_node_list,i) !grabing a node !using FOX function
        nodename = getLocalName(targetNode)  !finding its name !using FOX function
        if (nodename == targetNodeName) then !found our target node
            target_node_exists = .true.
            if (present(read_flag)) then
                read_flag =.true.
            endif
            exit
        endif
    enddo
    if (target_node_exists .eqv. .false.) then
        nullify(targetNode)
        if(present(mandatory)) then
            if (mandatory.eqv..false.) then
                outext='Could not find any node called "'//targetNodeName//'" in the xml file, ignoring'
                call Log%put(outext)
                if (present(read_flag)) then
                    read_flag =.false.
                endif
            else
                outext='Could not find any node called "'//targetNodeName//'" in the xml file, stoping'
                call Log%put(outext)
                stop
            endif
        else
            outext='Could not find any node called "'//targetNodeName//'" in the xml file, stoping'
            call Log%put(outext)
            stop
        endif
    endif
    end subroutine gotoNode

    end module simulation_xmlparser_mod
