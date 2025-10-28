
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
    use ModuleWaterQuality
    use geometry_mod

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
    type(Node), intent(in), pointer :: linksNode
    type(NodeList), pointer :: linkList
    type(Node), pointer :: anode
    type(Node), pointer :: aProp
    type(Node), pointer :: xmlProps                   !< .xml file handle
    type(NodeList), pointer :: propertyList        !< Node list for parameters
    type(string) :: xmlfilename, outext
    integer :: i, p, j
    type(string) :: att_name, att_val, tag, propKey, propValue
    type(string) :: sourceid, sourcetype, sourceprop

    linkList => getElementsByTagname(linksNode, "type")
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
        tag = tempSources%src(i)%prop%propertyType
        if (tag .ne. 'base') then
            call XMLReader%gotoNode(xmlProps,anode,tag) !finding the material type node
            tag = tempSources%src(i)%prop%propertySubType
            call XMLReader%gotoNode(anode,anode,tag)     !finding the actual material node
            propertyList => getElementsByTagname(anode, "property")        !searching for tags with the 'property' name
            call tempSources%src(i)%setPropertyNumber(getLength(propertyList))
            do j = 0, getLength(propertyList) - 1                          !extracting property tags one by one
                aProp => item(propertyList, j)
                att_name="key"
                call XMLReader%getLeafAttribute(aProp,att_name,propKey)
                att_name="value"
                call XMLReader%getLeafAttribute(aProp,att_name,propValue)
                call tempSources%src(i)%setPropertyAtribute(j+1, propKey, propValue)
            end do
            tag = 'particulate'
            att_name="value"
            call XMLReader%getNodeAttribute(anode, tag, att_name, att_val)
            call tempSources%src(i)%setPropertyBaseAtribute(tag, att_val)
            tag = 'density'
            call XMLReader%getNodeAttribute(anode, tag, att_name, att_val)
            call tempSources%src(i)%setPropertyBaseAtribute(tag, att_val)
            tag = 'radius'
            call XMLReader%getNodeAttribute(anode, tag, att_name, att_val)
            call tempSources%src(i)%setPropertyBaseAtribute(tag, att_val)
            tag = 'volume'
            call XMLReader%getNodeAttribute(anode, tag, att_name, att_val)
            call tempSources%src(i)%setPropertyBaseAtribute(tag, att_val)
            tag = 'area'
            call XMLReader%getNodeAttribute(anode, tag, att_name, att_val)
            call tempSources%src(i)%setPropertyBaseAtribute(tag, att_val)
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
    type(Node), intent(in), pointer :: case_node

    type(Node), pointer :: props_node          !< Single properties block to process
    type(string) :: outext
    type(string) :: tag, att_name

    tag="sourceTypes"    !the node we want
    call XMLReader%gotoNode(case_node,props_node,tag,mandatory =.false.)
    if (associated(props_node)) then
        tag="file"
        att_name="name"
        call XMLReader%getNodeAttribute(props_node, tag, att_name, Globals%Names%propsxmlfilename) !getting the file name from that tag
        outext='-->Properties to link to Sources found at '//Globals%Names%propsxmlfilename
        call Log%put(outext,.false.)
        tag="types"
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
    !> Sets the global variables responsible for controling field outputs.
    !> reads options from a xml file and adds a variable field to the print poll
    !> accordingly
    !> @param[in] exeNode
    !---------------------------------------------------------------------------
    subroutine setOutputFields(exeNode)
    type(Node), intent(in), pointer :: exeNode
    type(Node), pointer :: fileNode
    type(Node), pointer :: fieldNode
    type(Node), pointer :: outputFieldsFile
    type(NodeList), pointer :: outputFieldsList
    type(string) :: outext
    type(string) :: tag, att_name
    type(string) :: outputFieldsFilename
    type(string) :: fieldName, fieldOption
    type(string), dimension(:), allocatable :: fieldNameArray
    logical, dimension(:), allocatable :: toOutput
    integer :: i
    
    tag="outputFields" 
    call XMLReader%gotoNode(exeNode, fileNode, tag, mandatory =.false.)
    if (associated(fileNode)) then
        tag = "file"
        call XMLReader%gotoNode(fileNode, fileNode, tag)
        att_name="name"
        call XMLReader%getLeafAttribute(fileNode, att_name, outputFieldsFilename)
        !reading the file and building print/noprint array
        call XMLReader%getFile(outputFieldsFile, outputFieldsFilename)
        tag = "output"
        call XMLReader%gotoNode(outputFieldsFile ,outputFieldsFile, tag)
        outputFieldsList => getElementsByTagname(outputFieldsFile, "field")
        allocate(fieldNameArray(getLength(outputFieldsList)))
        allocate(toOutput(getLength(outputFieldsList)))
        toOutput = .false.
        do i = 0, getLength(outputFieldsList) - 1
            fieldNode => item(outputFieldsList, i)
            att_name="name"
            call XMLReader%getLeafAttribute(fieldNode, att_name, fieldName)
            att_name="output"
            call XMLReader%getLeafAttribute(fieldNode, att_name, fieldOption)
            fieldNameArray(i+1) = fieldName
            if (fieldOption == 'yes') toOutput(i+1) = .true.
        end do
    else
        outext='-->No output fields user override file, assuming basic settings for field output'
        call Log%put(outext,.false.)
    endif
    !calling the globals method to set the output variable field list
    if (.not.allocated(fieldNameArray)) allocate(fieldNameArray(0))
    if (.not.allocated(toOutput)) then
        allocate(toOutput(0))
        toOutput = .false.
    end if
    call Globals%Output%setOutputFields(fieldNameArray, toOutput)

    end subroutine setOutputFields
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> naming xml parser routine. Reads the naming file(s), opens the file(s) and
    !> stores the naming conventions for input files
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_naming(case_node)
    type(Node), intent(in), pointer :: case_node
    type(Node), pointer :: naming_node          !< Single naming block to process
    type(Node), pointer :: temp
    type(NodeList), pointer :: namingfileList
    type(string) :: outext
    type(string) :: tag, att_name
    type(string), allocatable, dimension(:) :: namingFilename
    integer :: i

    tag="variableNaming"    !the node we want
    call XMLReader%gotoNode(case_node,naming_node,tag,mandatory =.false.)
    if (associated(naming_node)) then
        namingfileList => getElementsByTagname(naming_node, "file")       !searching for tags with the 'namingfile' name
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
    !> @author Joao Sobrinho
    !> @brief
    !> Private external definitions parser routine. Builds the waterquality module definitions from the input .dat file.
    !> @param[in] case_node
    !---------------------------------------------------------------------------
    subroutine init_extDefFiles(case_node)
    type(Node), intent(in), pointer :: case_node
    type(Node), pointer :: execution_node       !< Single simdefs block to process
    type(string) :: outext
    integer :: i, id
    type(string) :: tag, att_name, att_val
    real(prec) :: waterQualityTimeStep
    integer :: STAT_CALL
    outext='-->Reading external definition files'
    call Log%put(outext,.false.)

    tag="externalDefinitionFiles"    !the node we want
    call XMLReader%gotoNode(case_node,execution_node,tag, mandatory =.false.)
    if (associated(execution_node)) then
        !Search for a water quality file
        tag="WQfile"
        att_name="name"
        call XMLReader%getNodeAttribute(execution_node, tag, att_name, att_val)

        call Globals%ExtImpFiles%setWaterQualityFileName(att_val) !Naming will need to change when other files are added
        
        if (Globals%ExtImpFiles%waterQualityFileName_hasValue) then
            !Start mohid water quality module.
            call StartWaterQuality(id, Globals%ExtImpFiles%waterQualityFileName%chars(), STAT = STAT_CALL, MohidLagr = .true.)
            if (id /= 0 .or. STAT_CALL /= 0) then
                outext='Failed to start mohid waterquality module - id was not 0 for some reason. stopping'
                call Log%put(outext)
                stop
            endif
            
            call GetDTWQM(0, DTSecond = waterQualityTimeStep)
            call Globals%SimDefs%setWqDt(waterQualityTimeStep) !Set waterquality module time step
        endif
        
        call Globals%ExtImpFiles%print()
    else
        outext='-->No mohid water quality file implementation found'
        call Log%put(outext,.false.)
    endif
    
    end subroutine init_extDefFiles
    
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
    type(NodeList), pointer :: pointList
    type(Node), pointer :: currPointNode
    type(string) :: outext
    type(string) :: tag, att_name, geoFileName, zMin, zMax
    integer :: i
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
    class is (polyline)
        pointList => getElementsByTagname(source_detail, "point")
        allocate(source_shape%point(getLength(pointList)-1))
        currPointNode => item(pointList, 0)
        call XMLReader%getLeafVector(currPointNode, source_shape%pt)
        do i = 1, getLength(pointList) - 1
            currPointNode => item(pointList, i)
            call XMLReader%getLeafVector(currPointNode, source_shape%point(i))
        end do        
    class is (sphere)
        tag='point'
        call XMLReader%getNodeVector(source_detail,tag,source_shape%pt)
        call extractDataAttribute(source_detail, "radius", source_shape%radius)
    class is (polygon)
        tag='file'
        att_name = 'name'
        call XMLReader%getNodeAttribute(source_detail, tag, att_name, geoFileName)
        zMin = notRead
        zMax = notRead
        tag='verticalBoundingBox'
        att_name = 'min'
        call XMLReader%getNodeAttribute(source_detail, tag, att_name, zMin, mandatory = .false.)
        att_name = 'max'
        call XMLReader%getNodeAttribute(source_detail, tag, att_name, zMax, mandatory = .false.)
        call Geometry%setPolygon(source_shape, geoFileName, zMin, zMax)
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
    type(NodeList), pointer :: activeList
    type(Node), pointer :: activeIntervalNode
    type(Node), pointer :: sourcedef
    type(Node), pointer :: source_node              !< Single source block to process
    type(Node), pointer :: source_detail
    type(Node), pointer :: source_ratefile, ratefileLeaf
    type(Node), pointer :: source_posifile, posifileLeaf
    integer :: i, j, k
    logical :: readflag
    integer :: id
    type(string) :: name, source_geometry, tag, att_name, att_val, att_val2, rate_file, posi_file
    real(prec) :: emitting_rate, rateScale, start, finish, temp, bottom_emission_depth, biofouling_rate, biofouling_start_after
    real(prec) :: tracer_volume
    logical :: emitting_fixed, posi_fixed, rateRead
    class(shape), allocatable :: source_shape
    type(vector) :: res
    real(prec), dimension(:,:), allocatable :: activeTimes
    
    readflag = .false.
    outext='-->Reading case Sources'
    call Log%put(outext,.false.)
    tag="sourceDefinitions"    !the node we want
    call XMLReader%gotoNode(case_node,sourcedef,tag)
    sourceList => getElementsByTagname(sourcedef, "source")
    !allocating the temporary source objects
    call tempSources%initialize(getLength(sourceList))

    do j = 0, getLength(sourceList) - 1
        rate_file = notSet
        rateRead = .false.
        posi_file = notSet
        posi_fixed = .true.
        source_node => item(sourceList,j)
        tag="setsource"
        att_name="id"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val)
        id=att_val%to_number(kind=1_I1P)
        att_name="name"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, name)
        !reading possible custom resolution
        tag="resolution"
        att_name="dp"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            res = att_val%to_number(kind=1._R8P)
        else
            call XMLReader%getNodeVector(source_node, tag, res, readflag, .false.)
            if (.not.readflag) res = 0.0
        end if
        !reading emission rate, need to check for options
        readflag = .false.
        tag="rate_dt"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            rateRead = .true.
            emitting_rate = 1.0/(att_val%to_number(kind=1._R8P)*Globals%SimDefs%dt)
            emitting_fixed = .true.
        end if
        tag="rate_seconds"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            tag="rate_trcPerEmission"
            att_name="value"
            readflag = .false.
            call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val2, readflag, .false.)
            if (readflag) then !Recieved number of tracers per rate_seconds
                rateRead = .true.
                emitting_rate = 1.0/att_val%to_number(kind=1._R8P) * att_val2%to_number(kind=1._R8P) !if 3600s, 5tracers and dt = 60 => 0.08 tracers per DT
                emitting_fixed = .true.
            else !if rate_seconds is used, then rate_trcPerEmission is needed
                rateRead = .false.
            endif
        else ! recieved rate in dt steps
            rateRead = .true.
            emitting_rate = 1.0/(att_val%to_number(kind=1._R8P)*Globals%SimDefs%dt)
            emitting_fixed = .true.
        end if
        tag="rate"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            rateRead = .true.
            emitting_rate = att_val%to_number(kind=1._R8P)
            emitting_fixed = .true.
        end if
        !Reading volume of each tracer (to be used with flow rate time series)
        tag="tracer_volume"
        att_name="tracer_volume"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            rateRead = .true.
            tracer_volume = att_val%to_number(kind=1._R8P)
        else
            tracer_volume = 0
        end if
        
        tag="rateTimeSeries"
        call XMLReader%gotoNode(source_node, source_ratefile, tag, mandatory =.false.)
        if (associated(source_ratefile)) then
            tag = "file"
            call XMLReader%gotoNode(source_ratefile, ratefileLeaf, tag)
            att_name="name"
            call XMLReader%getLeafAttribute(ratefileLeaf, att_name, att_val)
            rateRead = .true.
            rate_file = att_val
            emitting_fixed = .false.
            tag = "scale"
            att_name="value"
            call XMLReader%getNodeAttribute(source_ratefile, tag, att_name, att_val, readflag, mandatory = .false.)            
            if (readflag) then
                rateScale = att_val%to_number(kind=1._R8P)
            else
                rateScale = 1.0
            end if
        end if
        if (.not.rateRead) then
            outext='-->Source '//name//' (id = '//id// ') doesn''t have emission rate information. Possible options are [rate_dt, rate, rate_file]. Stoping'
            call Log%put(outext)
            stop
        end if
        !reading variable position tags, if any
        tag="positionTimeSeries"
        call XMLReader%gotoNode(source_node, source_posifile, tag, mandatory =.false.)
        if (associated(source_posifile)) then
            tag = "file"
            call XMLReader%gotoNode(source_posifile, posifileLeaf, tag)
            att_name="name"
            call XMLReader%getLeafAttribute(posifileLeaf, att_name, att_val)
            posi_file = att_val
            posi_fixed = .false.
        end if
        !reading bottom emission type of source
        readflag = .false.
        tag="bottom_emission_depth"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            bottom_emission_depth = att_val%to_number(kind=1._R8P)
        else
            bottom_emission_depth = 0
        end if
        !reading biofouling growth rate of this source
        readflag = .false.
        tag="biofouling_rate"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            biofouling_rate = att_val%to_number(kind=1._R8P)
        else
            biofouling_rate = 0
        end if
        
        !reading biofouling incubation period for this source
        readflag = .false.
        tag="biofouling_start_after"
        att_name="value"
        call XMLReader%getNodeAttribute(source_node, tag, att_name, att_val, readflag, .false.)
        if (readflag) then
            biofouling_start_after = att_val%to_number(kind=1._R8P)
        else
            biofouling_start_after = 0
        end if
        
        !Possible list of active intervals, and these might be in absolute dates or relative to sim time
        activeList => getElementsByTagname(source_node, "active")
        allocate(activeTimes(getLength(activeList),2))
        do k = 0, getLength(activeList) - 1
            activeIntervalNode => item(activeList, k)
            att_name="start"
            call XMLReader%getLeafAttribute(activeIntervalNode,att_name,att_val)
            temp = Utils%getRelativeTimeFromString(att_val, Globals%SimTime%StartDate)
            activeTimes(k+1,1) = temp
            att_name="end"
            call XMLReader%getLeafAttribute(activeIntervalNode,att_name,att_val)
            if (att_val == 'end') then
                temp = Globals%Parameters%TimeMax
            else
                temp = Utils%getRelativeTimeFromString(att_val, Globals%SimTime%StartDate)
            end if
            activeTimes(k+1,2) = temp
        end do
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
        call tempSources%src(j+1)%initialize(id, name, emitting_rate, emitting_fixed, rate_file, rateScale, &
                                             posi_fixed, posi_file, bottom_emission_depth, biofouling_rate, &
                                             biofouling_start_after, tracer_volume, activeTimes, source_geometry, &
                                             source_shape, res)
       
        deallocate(activeTimes)
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
    type(Node), intent(in), pointer :: case_node
    type(NodeList), pointer :: defsList       !< Node list for simdefs
    type(NodeList), pointer :: beachingChildren       !< Node list for source node children nodes
    type(NodeList), pointer :: beachingAreasList       !< Node list for source node children nodes
    type(Node), pointer :: simdefs_node       !< Single simdefs block to process
    type(Node), pointer :: beaching_node       !< Single beaching block to process
    type(Node), pointer :: beachingArea_node   !< Single beaching block to process
    type(Node), pointer :: beaching_detail
    type(string) :: outext
    integer :: i, id, j, coastType, unbeach, runUpEffect, runUpEffectUnbeach
    type(string) :: pts(2), tag, att_name, att_val, beaching_geometry
    real(prec) :: probability, waterColumnThreshold, beachTimeScale, unbeachTimeScale, beachSlope
    class(shape), allocatable :: beaching_shape
    type(vector) :: coords
    logical :: read_flag
    read_flag = .false.
    coords = 0.0
    outext='-->Reading case simulation definitions'
    call Log%put(outext,.false.)

    tag="simulation"    !the node we want
    call XMLReader%gotoNode(case_node,simdefs_node,tag)
    tag="resolution"
    att_name="dp"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        coords = att_val%to_number(kind=1._R8P)
    else
        call XMLReader%getNodeVector(simdefs_node, tag, coords)
    end if
    call Globals%SimDefs%setdp(coords)
    tag="timestep"
    att_name="dt"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val)
    call Globals%SimDefs%setdt(att_val)
    call Globals%Constants%setSmallDt(Globals%SimDefs%dt)
    pts=(/ 'BoundingBoxMin', 'BoundingBoxMax'/) !strings to search for
    do i=1, size(pts)
        call XMLReader%getNodeVector(simdefs_node, pts(i), coords)
        call Globals%SimDefs%setboundingbox(pts(i), coords)
    enddo
    
    tag="VerticalVelMethod"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setVerticalVelMethod(att_val)
    endif
    
    tag="DiffusionMethod"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setDiffusionMethod(att_val)
    endif

    tag="ResuspensionCriticalShearMethod"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setResuspensionCriticalShearMethod(att_val)
    endif
    
    tag="RemoveLandTracer"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setRemoveLandTracer(att_val)
    endif
    
    tag="BathyminNetcdf"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setBathyminNetcdf(att_val)
    endif

    tag="RugosityinNetcdf"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setRugosityinNetcdf(att_val)
    endif
    
    tag="TracerMaxAge"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%settracerMaxAge(att_val)
    endif
    
    tag="Temperature_add_offset"
    att_name="value"
    call XMLReader%getNodeAttribute(simdefs_node, tag, att_name, att_val, read_flag, .false.)
    if (read_flag) then
        call Globals%SimDefs%setTemperature_add_offset(att_val)
    endif
    
    tag="Beaching"
    call XMLReader%gotoNode(simdefs_node,beaching_node,tag,read_flag, .false.)
    
    if (read_flag) then
        tag="BeachingEnabled"
        att_name="value"
        call XMLReader%getNodeAttribute(beaching_node, tag, att_name, att_val, read_flag, .true.)
        if (read_flag) then
            call Globals%SimDefs%setFreeLitterAtBeaching(att_val)
        endif
        
        !Now get all the beaching areas
        beachingAreasList => getElementsByTagname(beaching_node, "BeachingArea")
        !allocating the temporary BeachingAreas objects
        call Globals%BeachingAreas%initialize(getLength(beachingAreasList))
        
        do j = 0, getLength(beachingAreasList) - 1
            beachingArea_node => item(beachingAreasList,j)
            tag="BeachAreaID"
            att_name="id"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val)
            id=att_val%to_number(kind=1_I1P)
            
            !reading possible custom resolution
            tag="CoastType"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            if (read_flag) then
                coastType = att_val%to_number(kind=1_I1P)
            else
                coastType = 1
            end if
            !reading emission rate, need to check for options
            
            tag="Probability"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            if (read_flag) then
                probability = att_val%to_number(kind=1._R8P)
            else
                probability = 0.0
            end if
            if (probability < 0.0 .or. probability > 1.0) then
                outext='-->Beaching Area id = '//id// ' has probability < 0 or > 1. Should be between 0 and 1. Stoping'
                call Log%put(outext)
                stop
            end if
            
            tag="WaterColumnThreshold"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            waterColumnThreshold = 0.1
            if (read_flag) then
                waterColumnThreshold = att_val%to_number(kind=1._R8P)
            end if
            
            
            tag="BeachTimeScale"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            beachTimeScale = 3600.0
            if (read_flag) then
                beachTimeScale = att_val%to_number(kind=1._R8P)
            end if
        
            tag="Unbeach"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            unbeach = 1
            if (read_flag) then
                unbeach = att_val%to_number(kind=1_I1P)
            end if
            
            tag="UnbeachTimeScale"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            unbeachTimeScale = 3600.0
            if (read_flag) then
                unbeachTimeScale = att_val%to_number(kind=1._R8P)
            end if
            
            tag="RunUpEffect"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            runUpEffect = 0
            if (read_flag) then
                runUpEffect = att_val%to_number(kind=1_I1P)
            end if
            
            tag="BeachSlope"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            beachSlope = 0.01
            if (read_flag) then
                beachSlope = att_val%to_number(kind=1._R8P)
            end if
        
            tag="RunUpEffectUnbeach"
            att_name="value"
            call XMLReader%getNodeAttribute(beachingArea_node, tag, att_name, att_val, read_flag, .false.)
            runUpEffectUnbeach = 0
            if (read_flag) then
                runUpEffectUnbeach = att_val%to_number(kind=1_I1P)
            end if
        
            beachingChildren => getChildNodes(beachingArea_node) !getting all of the nodes bellow the main source node (all of it's private info)
            do i=0, getLength(beachingChildren)-1
                beaching_detail => item(beachingChildren,i) !grabing a node
                beaching_geometry = getLocalName(beaching_detail)  !finding its name
                if (Geometry%isPolygon(beaching_geometry)) then  !if the node is a polygon
                    call Geometry%allocateShape(beaching_geometry,beaching_shape)                
                    call read_xml_geometry(beachingArea_node,beaching_detail,beaching_shape)
                    exit
                end if
            end do
            
            !initializing Beach area
            call Globals%BeachingAreas%beachArea(j+1)%initialize(id, coastType, probability, waterColumnThreshold, &
                                                 beachTimeScale, unbeach, unbeachTimeScale, runUpEffect, &
                                                 beachSlope, runUpEffectUnbeach, beaching_geometry, beaching_shape)
       
            deallocate(beaching_shape)
        enddo
        
    endif
    
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

    tag="constants"    !the node we want
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
        tag="BeachingLevel"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBeachingLevel(att_val)
        endif
        tag="BeachingStopProb"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBeachingStopProb(att_val)
        endif
        tag="DiffusionCoeff"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setDiffusionCoeff(att_val)
        endif
        tag="WindDragCoeff"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setWindDragCoeff(att_val)
        endif
        tag="ResuspensionCoeff"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setResuspensionCoeff(att_val)
        endif                
        tag="MeanDensity"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setMeanDensity(att_val)
        endif
        tag="MeanKViscosity"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setMeanKVisco(att_val)
        endif
        tag="AddBottomCell"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setAddBottomCell(att_val)
        endif
        tag="Rugosity"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setRugosity(att_val)
        endif
		tag="BedLoadThickness"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBedLoadThickness(att_val)
        endif
		tag="ResuspensionProb"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setResuspensionProb(att_val)
        endif
		tag="ResuspensionResidenceTime"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setResuspensionResidenceTime(att_val)
        endif
		tag="ResuspensionCriticalShear1"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setResuspensionCriticalShear1(att_val)
        endif
		tag="ResuspensionCriticalShear2"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setResuspensionCriticalShear2(att_val)
        endif
        tag="CriticalShearErosion"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setCritical_Shear_Erosion(att_val)
        endif
        tag="ToptBMin"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setTOptBacteriaMin(att_val)
        endif
        tag="ToptBMax"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setTOptBacteriaMax(att_val)
        endif
        tag="TBMin"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setTBacteriaMin(att_val)
        endif
        tag="TBMax"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setTBacteriaMax(att_val)
        endif
        tag="BK1"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBK1(att_val)
        endif
        tag="BK2"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBK2(att_val)
        endif
        tag="BK3"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBK3(att_val)
        endif
        tag="BK4"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setBK4(att_val)
        endif
        tag="Lin0_DegradationRate"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setLin0_DegradationRate(att_val)
        endif
        tag="Exp0_DegradationRate"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setExp0_DegradationRate(att_val)
        endif
        tag="Exp1_DegradationRate"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setExp1_DegradationRate(att_val)
        endif
        tag="Exp2_DegradationRate"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setExp2_DegradationRate(att_val)
        endif
        tag="ConstantMixingLength"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setMixingLength(att_val)
        endif
        tag="VarVelHX"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setVarVelHX(att_val)
        endif
        tag="VarVelH"
        att_name="value"
        call XMLReader%getNodeAttribute(constants_node, tag, att_name, att_val,readflag,.false.)
        if (readflag) then
            call Globals%Constants%setVarVelH(att_val)
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
    Globals%Names%inputsXmlFilename = Globals%Names%outpath//xmlfilename%basename(extension='.xml')//'_inputs.xml'
    outext='->Input files index file name is '//Globals%Names%inputsXmlFilename
    call Log%put(outext)

    tag="case"          !base document node
    call XMLReader%gotoNode(xmldoc,execution_node,tag)
    tag="execution"     !finding execution node
    call XMLReader%gotoNode(execution_node,execution_node,tag)
    tag="case"          !base document node
    call XMLReader%gotoNode(xmldoc,case_node,tag)
    tag="caseDefinitions"     !finding execution node
    call XMLReader%gotoNode(case_node,case_node,tag)

    ! building the simulation basic structures according to the case definition file
    ! every other structure in the simulation is built from these, i.e., not defined by the user directly
    call init_parameters(execution_node)
    call init_naming(execution_node)
    call init_extDefFiles(execution_node)
    call setOutputFields(execution_node)
    call init_caseconstants(case_node)
    call init_simdefs(case_node)
    call init_sources(case_node)
    call init_properties(case_node)    

    !setting the number of blocks to the correct ammount of selected threads
    Globals%SimDefs%numblocks = Globals%Parameters%numOPMthreads

    call XMLReader%closeFile(xmldoc)

    end subroutine InitFromXml

    end module simulationInitialize_mod
