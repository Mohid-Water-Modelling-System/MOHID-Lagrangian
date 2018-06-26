
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
    !> Module with the simulation xml parsing related definitions and routines.
    !------------------------------------------------------------------------------

    module simulation_xmlparser_mod

    use FoX_dom
    use common_modules

    implicit none
    private

    !Public access procedures
    public :: ReadXMLatt, ReadXMLvector, GotoChildNode

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private attribute xml parser routine. In the format <Tag att_name="att_value"
    !
    !> @param[in] xmlnode, tag, vec, mandatory
    !---------------------------------------------------------------------------
    subroutine ReadXMLatt(xmlnode, tag, att_name, att_value, read_flag, mandatory)
    implicit none
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
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info)
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node
        nodename = getLocalName(nodedetail)  !finding its name
        if (nodename == tag) then
            validtag=.true.
            exit
        endif
    enddo
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name
        nodedetail => item(target_node_list, 0)
        call extractDataAttribute(nodedetail, att_name%chars(), att_value_chars)
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

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private vector xml parser routine. Vector must be in format <Tag x="vec%x" y="vec%y" z="vec%z" />
    !
    !> @param[in] xmlnode, tag, vec, mandatory
    !---------------------------------------------------------------------------
    subroutine ReadXMLvector(xmlnode, tag, vec, read_flag, mandatory)
    implicit none
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
    nodeChildren => getChildNodes(xmlnode) !getting all of the nodes bellow the main source node (all of it's private info)
    do i=0, getLength(nodeChildren)-1
        nodedetail => item(nodeChildren,i) !grabing a node
        nodename = getLocalName(nodedetail)  !finding its name
        if (nodename == tag) then
            validtag =.true.
            exit
        endif
    enddo
    if (validtag) then
        target_node_list => getElementsByTagname(xmlnode, tag%chars())   !searching for tags with the given name
        nodedetail => item(target_node_list, 0)
        call extractDataAttribute(nodedetail, "x", vec%x)
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
    end subroutine


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private routine to retrieve a node within a node.
    !> Returns a nullifyed pointer if not found, stops if mandatory.
    !
    !> @param[in] currentNode, targetNode, targetNodeName, mandatory
    !---------------------------------------------------------------------------
    subroutine GotoChildNode(currentNode, targetNode, targetNodeName, read_flag, mandatory)
    implicit none
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
    target_node_list => getChildNodes(currentNode)
    do i=0, getLength(target_node_list)-1
        targetNode => item(target_node_list,i) !grabing a node
        nodename = getLocalName(targetNode)  !finding its name
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

    end subroutine

  end module simulation_xmlparser_mod
