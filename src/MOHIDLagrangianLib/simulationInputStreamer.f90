    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_input_streamer
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an input file reader class with an object exposable to the Simulation
    !> This class is in charge of selectig the correct reader for the selected input
    !> file format and controling the respective reader.
    !------------------------------------------------------------------------------

    module simulationInputStreamer_mod

    use common_modules
    use xmlParser_mod
    use netcdfParser_mod
    use fieldTypes_mod
    use background_mod
    use blocks_mod
    use boundingbox_mod

    use FoX_dom

    implicit none
    private

    type :: inputFileModel_class !< Input file model class
        type(string) :: name        !< name of the file
        real(prec) :: startTime     !< starting time of the data on the file
        real(prec) :: endTime       !< ending time of the data on the file
        logical :: used             !< flag that indicates the file is no longer to be read
        logical :: toRead
    contains
    procedure, private :: setReadStatus
    end type inputFileModel_class

    type :: input_streamer_class        !< Input Streamer class
        logical :: useInputFiles
        type(inputFileModel_class), allocatable, dimension(:) :: currentsInputFile      !< array of input file metadata for currents
        type(inputFileModel_class), allocatable, dimension(:) :: windsInputFile         !< array of input file metadata for winds
        type(inputFileModel_class), allocatable, dimension(:) :: wavesInputFile         !< array of input file metadata for waves
        type(inputFileModel_class), allocatable, dimension(:) :: waterPropsInputFile    !< array of input file metadata for water properties
        real(prec) :: lastCurrentsReadTime
        real(prec) :: lastWindsReadTime
        real(prec) :: lastWavesReadTime
        real(prec) :: lastWaterPropsReadTime
        integer :: nFileTypes
        real(prec) :: bufferSize                                               !< half of the biggest tail of data behind current time
        integer :: currentsBkgIndex, windsBkgIndex, wavesBkgIndex, waterPropsBkgIndex
    contains
    procedure :: initialize => initInputStreamer
    procedure :: loadDataFromStack
    procedure, private :: getCurrentsFile
    procedure, private :: getWindsFile
    procedure, private :: getWavesFile
    procedure, private :: getWaterPropsFile
    procedure, private :: resetReadStatus
    procedure :: print => printInputStreamer
    end type input_streamer_class

    !Public access vars
    public :: input_streamer_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Checks if file should be read and sets appropriate flag
    !> @param[in] self, lastReadTime, bufferSize
    !---------------------------------------------------------------------------
    subroutine setReadStatus(self, lastReadTime, bufferSize)
    class(inputFileModel_class), intent(inout) :: self
    real(prec), intent(in) :: lastReadTime
    real(prec), intent(in) :: bufferSize
    if (self%endTime >= Globals%SimTime%CurrTime) then
        if (self%startTime <= Globals%SimTime%CurrTime + bufferSize) then
            if (.not.self%used) self%toRead = .true.
            if (lastReadTime >= self%endTime) self%toRead = .false.
        end if
    end if
    end subroutine setReadStatus

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a
    !> currents input file
    !> @param[in] self, fileName
    !---------------------------------------------------------------------------
    type(background_class) function getCurrentsFile(self, fileName)
    class(input_streamer_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), allocatable, dimension(:) :: varList
    logical, allocatable, dimension(:) :: syntecticVar
    type(ncReader_class) :: ncReader

    allocate(varList(6))
    allocate(syntecticVar(6))
    varList(1) = Globals%Var%u
    syntecticVar(1) = .false.
    varList(2) = Globals%Var%v
    syntecticVar(2) = .false.
    varList(3) = Globals%Var%w
    syntecticVar(3) = .false.
    varList(4) = Globals%Var%landIntMask
    syntecticVar(4) = .true. 
    varList(5) = Globals%Var%resolution
    syntecticVar(5) = .true.
    varList(6) = Globals%Var%bathymetry
    syntecticVar(6) = .true.

    !need to send to different readers here if different file formats
    getCurrentsFile = ncReader%getFullFile(fileName, varList, syntecticVar)
    call getCurrentsFile%makeLandMaskField()
    call getCurrentsFile%makeResolutionField()
    call getCurrentsFile%makeBathymetryField()

    end function getCurrentsFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a
    !> meteorology input file (winds)
    !> @param[in] self, fileName
    !---------------------------------------------------------------------------
    type(background_class) function getWindsFile(self, fileName)
    class(input_streamer_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), allocatable, dimension(:) :: varList
    logical, allocatable, dimension(:) :: syntecticVar
    type(ncReader_class) :: ncReader

    allocate(varList(3))
    allocate(syntecticVar(3))
    varList(1) = Globals%Var%u10
    syntecticVar(1) = .false.
    varList(2) = Globals%Var%v10
    syntecticVar(2) = .false.
    varList(3) = Globals%Var%rad
    syntecticVar(3) = .false.

    !need to send to different readers here if different file formats
    getWindsFile = ncReader%getFullFile(fileName, varList, syntecticVar)

    end function getWindsFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a
    !> waves input file (stokes drift velocity)
    !> @param[in] self, fileName
    !---------------------------------------------------------------------------
    type(background_class) function getWavesFile(self, fileName)
    class(input_streamer_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), allocatable, dimension(:) :: varList
    logical, allocatable, dimension(:) :: syntecticVar
    type(ncReader_class) :: ncReader

    allocate(varList(2))
    allocate(syntecticVar(2))
    varList(1) = Globals%Var%vsdx
    syntecticVar(1) = .false.
    varList(2) = Globals%Var%vsdy
    syntecticVar(2) = .false.

    !need to send to different readers here if different file formats
    getWavesFile = ncReader%getFullFile(fileName, varList, syntecticVar)

    end function getWavesFile
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a
    !> water properties input file (temperature, density, salinity)
    !> @param[in] self, fileName
    !---------------------------------------------------------------------------
    type(background_class) function getWaterPropsFile(self, fileName)
    class(input_streamer_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), allocatable, dimension(:) :: varList
    logical, allocatable, dimension(:) :: syntecticVar
    type(ncReader_class) :: ncReader

    allocate(varList(3))
    allocate(syntecticVar(3))
    varList(1) = Globals%Var%temp
    syntecticVar(1) = .false.
    varList(2) = Globals%Var%sal
    syntecticVar(2) = .false.
    varList(3) = Globals%Var%density
    syntecticVar(3) = .false.

    !need to send to different readers here if different file formats
    getWaterPropsFile = ncReader%getFullFile(fileName, varList, syntecticVar)

    end function getWaterPropsFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> loads data from files and populates the backgrounds accordingly
    !> @param[in] self, bbox, blocks
    !---------------------------------------------------------------------------
    subroutine loadDataFromStack(self, bBox, blocks)
    class(input_streamer_class), intent(inout) :: self
    type(boundingbox_class), intent(in) :: bBox            !< Case bounding box
    type(block_class), dimension(:), intent(inout) :: blocks  !< Case Blocks
    type(background_class) :: tempBkgd, tempBkgd2, tempBkgd3
    integer :: i, j, k
    integer :: fNumber
    real(prec) :: tempTime(2)
    logical :: appended

    if (self%useInputFiles) then
        call self%resetReadStatus()
        !check what files on the stack are to read to backgrounds
        if (allocated(self%currentsInputFile)) then
            do i=1, size(self%currentsInputFile)
                call self%currentsInputFile(i)%setReadStatus(self%lastCurrentsReadTime, self%BufferSize)
            end do
        end if
        if (allocated(self%windsInputFile)) then
            do i=1, size(self%windsInputFile)
                call self%windsInputFile(i)%setReadStatus(self%lastWindsReadTime, self%BufferSize)
            end do
        end if
        if (allocated(self%wavesInputFile)) then
            do i=1, size(self%wavesInputFile)
                call self%wavesInputFile(i)%setReadStatus(self%lastWavesReadTime, self%BufferSize)
            end do
        end if

        if (allocated(self%waterPropsInputFile)) then
            do i=1, size(self%waterPropsInputFile)
                call self%waterPropsInputFile(i)%setReadStatus(self%lastWaterPropsReadTime, self%BufferSize)
            end do
        end if
        !read selected files
        if (allocated(self%currentsInputFile)) then
            do i=1, size(self%currentsInputFile)
                if (self%currentsInputFile(i)%toRead) then
                    !import data to temporary background
                    tempBkgd = self%getCurrentsFile(self%currentsInputFile(i)%name)
                    self%currentsInputFile(i)%used = .true.
                    do j=1, size(blocks)
                        !slice data by block and either join to existing background or add a new one                        
                        if (blocks(j)%Background(self%currentsBkgIndex)%initialized) then
                            tempBkgd2 = tempBkgd%getHyperSlab(blocks(j)%extents)
                            call blocks(j)%Background(self%currentsBkgIndex)%append(tempBkgd2, appended)
                            call tempBkgd2%finalize()
                        end if
                        if (.not.blocks(j)%Background(self%currentsBkgIndex)%initialized) blocks(j)%Background(self%currentsBkgIndex) = tempBkgd%getHyperSlab(blocks(j)%extents)
                        !save last time already loaded
                        tempTime = blocks(j)%Background(self%currentsBkgIndex)%getDimExtents(Globals%Var%time)
                        if (self%lastCurrentsReadTime == -1.0) self%lastCurrentsReadTime = tempTime(2)
                        self%lastCurrentsReadTime = min(tempTime(2), self%lastCurrentsReadTime)
                    end do
                    !clean out the temporary background data
                    call tempBkgd%finalize()
                end if
            end do
        end if
        if (allocated(self%windsInputFile)) then
            do i=1, size(self%windsInputFile)
                if (self%windsInputFile(i)%toRead) then
                    !import data to temporary background
                    tempBkgd = self%getWindsFile(self%windsInputFile(i)%name)
                    self%windsInputFile(i)%used = .true.
                    do j=1, size(blocks)
                        !slice data by block and either join to existing background or add a new one
                        if (blocks(j)%Background(self%windsBkgIndex)%initialized) then
                            tempBkgd2 = tempBkgd%getHyperSlab(blocks(j)%extents)
                            call blocks(j)%Background(self%windsBkgIndex)%append(tempBkgd2, appended)
                            call tempBkgd2%finalize()
                        end if
                        if (.not.blocks(j)%Background(self%windsBkgIndex)%initialized) blocks(j)%Background(self%windsBkgIndex) = tempBkgd%getHyperSlab(blocks(j)%extents)
                        !save last time already loaded
                        tempTime = blocks(j)%Background(self%windsBkgIndex)%getDimExtents(Globals%Var%time)
                        if (self%lastWindsReadTime == -1.0) self%lastWindsReadTime = tempTime(2)
                        self%lastWindsReadTime = min(tempTime(2), self%lastWindsReadTime)
                    end do
                    !clean out the temporary background data
                    call tempBkgd%finalize()
                end if
            end do
        end if
        if (allocated(self%wavesInputFile)) then
            do i=1, size(self%wavesInputFile)
                if (self%wavesInputFile(i)%toRead) then
                    !import data to temporary background
                    tempBkgd = self%getWavesFile(self%wavesInputFile(i)%name)
                    self%wavesInputFile(i)%used = .true.
                    do j=1, size(blocks)
                        !slice data by block and either join to existing background or add a new one
                        if (blocks(j)%Background(self%wavesBkgIndex)%initialized) then
                            tempBkgd2 = tempBkgd%getHyperSlab(blocks(j)%extents)
                            call blocks(j)%Background(self%wavesBkgIndex)%append(tempBkgd2, appended)
                            call tempBkgd2%finalize()
                        end if
                        if (.not.blocks(j)%Background(self%wavesBkgIndex)%initialized) blocks(j)%Background(self%wavesBkgIndex) = tempBkgd%getHyperSlab(blocks(j)%extents)
                        !save last time already loaded
                        tempTime = blocks(j)%Background(self%wavesBkgIndex)%getDimExtents(Globals%Var%time)
                        if (self%lastWavesReadTime == -1.0) self%lastWavesReadTime = tempTime(2)
                        self%lastWavesReadTime = min(tempTime(2), self%lastWavesReadTime)
                    end do
                    !clean out the temporary background data
                    call tempBkgd%finalize()
                end if
            end do
        end if
        if (allocated(self%waterPropsInputFile)) then
            do i=1, size(self%waterPropsInputFile)
                if (self%waterPropsInputFile(i)%toRead) then
                    !import data to temporary background
                    tempBkgd = self%getWaterPropsFile(self%waterPropsInputFile(i)%name)
                    self%waterPropsInputFile(i)%used = .true.
                    do j=1, size(blocks)
                        !slice data by block and either join to existing background or add a new one
                        if (blocks(j)%Background(self%waterPropsBkgIndex)%initialized) then
                            tempBkgd2 = tempBkgd%getHyperSlab(blocks(j)%extents)
                            call blocks(j)%Background(self%waterPropsBkgIndex)%append(tempBkgd2, appended)
                            call tempBkgd2%finalize()
                        end if
                        if (.not.blocks(j)%Background(self%waterPropsBkgIndex)%initialized) blocks(j)%Background(self%waterPropsBkgIndex) = tempBkgd%getHyperSlab(blocks(j)%extents)
                        !save last time already loaded
                        tempTime = blocks(j)%Background(self%waterPropsBkgIndex)%getDimExtents(Globals%Var%time)
                        if (self%lastWaterPropsReadTime == -1.0) self%lastWaterPropsReadTime = tempTime(2)
                        self%lastWaterPropsReadTime = min(tempTime(2), self%lastWaterPropsReadTime)
                    end do
                    !clean out the temporary background data
                    call tempBkgd%finalize()
                end if
            end do
        end if
    end if

    end subroutine loadDataFromStack

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the input writer object, imports metadata on input files
    !> @param[in] self, blocks
    !---------------------------------------------------------------------------
    subroutine initInputStreamer(self, blocks)
    class(input_streamer_class), intent(inout) :: self
    type(block_class), dimension(:), intent(inout) :: blocks  !< Case Blocks
    type(Node), pointer :: xmlInputs           !< .xml file handle
    type(Node), pointer :: typeNode
    type(Node), pointer :: fileNode
    type(NodeList), pointer :: fileList
    type(string) :: tag, att_name, att_val
    type(string), allocatable, dimension(:) :: fileNames
    integer :: i, nBkg

    self%bufferSize = Globals%Parameters%BufferSize
    self%lastCurrentsReadTime = -1.0
    self%lastWindsReadTime = -1.0
    self%lastWavesReadTime = -1.0
    self%nFileTypes = 0
    self%currentsBkgIndex = 0
    self%windsBkgIndex = 0
    self%wavesBkgIndex = 0
    self%waterPropsBkgIndex = 0
    nBkg = 0

    call XMLReader%getFile(xmlInputs,Globals%Names%inputsXmlFilename, mandatory = .false.)
    if (associated(xmlInputs)) then
        self%useInputFiles = .true.
        !Go to the file_collection node
        tag = "file_collection"
        call XMLReader%gotoNode(xmlInputs,xmlInputs,tag)

        !For currents data
        tag = Globals%DataTypes%currents
        call XMLReader%gotoNode(xmlInputs,typeNode,tag, mandatory=.false.)
        if (associated(typeNode)) then
            fileList => getElementsByTagname(typeNode, "file")
            allocate(fileNames(getLength(fileList)))
            allocate(self%currentsInputFile(getLength(fileList)))
            do i = 0, getLength(fileList) - 1
                fileNode => item(fileList, i)
                tag="name"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, fileNames(i+1))
                self%currentsInputFile(i+1)%name = fileNames(i+1)
                tag="startTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%currentsInputFile(i+1)%startTime = att_val%to_number(kind=1._R8P)
                tag="endTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%currentsInputFile(i+1)%endTime = att_val%to_number(kind=1._R8P)
                self%currentsInputFile(i+1)%used = .false.
            end do
            deallocate(fileNames)
            nBkg = nBkg + 1
            self%currentsBkgIndex = nBkg
        end if

        !For wind data
        tag = Globals%DataTypes%winds
        call XMLReader%gotoNode(xmlInputs,typeNode,tag, mandatory=.false.)
        if (associated(typeNode)) then
            fileList => getElementsByTagname(typeNode, "file")
            allocate(fileNames(getLength(fileList)))
            allocate(self%windsInputFile(getLength(fileList)))
            do i = 0, getLength(fileList) - 1
                fileNode => item(fileList, i)
                tag="name"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, fileNames(i+1))
                self%windsInputFile(i+1)%name = fileNames(i+1)
                tag="startTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%windsInputFile(i+1)%startTime = att_val%to_number(kind=1._R8P)
                tag="endTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%windsInputFile(i+1)%endTime = att_val%to_number(kind=1._R8P)
                self%windsInputFile(i+1)%used = .false.
            end do
            deallocate(fileNames)
            nBkg = nBkg + 1
            self%windsBkgIndex = nBkg
        end if

        !For wave data
        tag = Globals%DataTypes%waves
        call XMLReader%gotoNode(xmlInputs,typeNode,tag, mandatory=.false.)
        if (associated(typeNode)) then
            fileList => getElementsByTagname(typeNode, "file")
            allocate(fileNames(getLength(fileList)))
            allocate(self%wavesInputFile(getLength(fileList)))
            do i = 0, getLength(fileList) - 1
                fileNode => item(fileList, i)
                tag="name"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, fileNames(i+1))
                self%wavesInputFile(i+1)%name = fileNames(i+1)
                tag="startTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%wavesInputFile(i+1)%startTime = att_val%to_number(kind=1._R8P)
                tag="endTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%wavesInputFile(i+1)%endTime = att_val%to_number(kind=1._R8P)
                self%wavesInputFile(i+1)%used = .false.
            end do
            nBkg = nBkg + 1
            self%wavesBkgIndex = nBkg
        end if
        
        !For water properties data
        tag = Globals%DataTypes%waterProps
        call XMLReader%gotoNode(xmlInputs,typeNode,tag, mandatory=.false.)
        if (associated(typeNode)) then
            fileList => getElementsByTagname(typeNode, "file")
            allocate(fileNames(getLength(fileList)))
            allocate(self%waterPropsInputFile(getLength(fileList)))
            do i = 0, getLength(fileList) - 1
                fileNode => item(fileList, i)
                tag="name"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, fileNames(i+1))
                self%waterPropsInputFile(i+1)%name = fileNames(i+1)
                tag="startTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%waterPropsInputFile(i+1)%startTime = att_val%to_number(kind=1._R8P)
                tag="endTime"
                att_name="value"
                call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
                self%waterPropsInputFile(i+1)%endTime = att_val%to_number(kind=1._R8P)
                self%waterPropsInputFile(i+1)%used = .false.
            end do
            nBkg = nBkg + 1
            self%waterPropsBkgIndex = nBkg
        end if
        
    else
        self%useInputFiles = .false.
    end if
    !allocating the necessary background array in every block
    do i=1, size(blocks)
        allocate(blocks(i)%Background(nBkg))
    end do
    end subroutine initInputStreamer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> resets input files read status
    !---------------------------------------------------------------------------
    subroutine resetReadStatus(self)
    class(input_streamer_class), intent(inout) :: self
    integer :: i
    if (allocated(self%currentsInputFile)) then
        do i=1, size(self%currentsInputFile)
            self%currentsInputFile(i)%toRead = .false.
        end do
    end if
    if (allocated(self%windsInputFile)) then
        do i=1, size(self%windsInputFile)
            self%windsInputFile(i)%toRead = .false.
        end do
    end if
    if (allocated(self%wavesInputFile)) then
        do i=1, size(self%wavesInputFile)
            self%wavesInputFile(i)%toRead = .false.
        end do
    end if
    if (allocated(self%waterPropsInputFile)) then
        do i=1, size(self%waterPropsInputFile)
            self%waterPropsInputFile(i)%toRead = .false.
        end do
    end if
    end subroutine resetReadStatus

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Prints the input writer object and metadata on input files
    !---------------------------------------------------------------------------
    subroutine printInputStreamer(self)
    class(input_streamer_class), intent(in) :: self
    type(string) :: outext, temp_str
    integer :: i
    logical :: written
    written = .false.
    outext = '-->Input streamer stack:'//new_line('a')
    if (allocated(self%currentsInputFile)) then
        outext = outext//'--->'//Globals%DataTypes%currents%startcase()//' data '
        do i=1, size(self%currentsInputFile)
            outext = outext//new_line('a')
            outext = outext//'---->File '//self%currentsInputFile(i)%name
        end do
        written = .true.
    end if
    if (allocated(self%windsInputFile)) then
        if (written) outext = outext//new_line('a')
        outext = outext//'--->'//Globals%DataTypes%winds%startcase()//' data '
        do i=1, size(self%windsInputFile)
            outext = outext//new_line('a')
            outext = outext//'---->File '//self%windsInputFile(i)%name
        end do
        written = .true.
    end if
    if (allocated(self%wavesInputFile)) then
        if (written) outext = outext//new_line('a')
        outext = outext//'--->'//Globals%DataTypes%waves%startcase()//' data '
        do i=1, size(self%wavesInputFile)
            outext = outext//new_line('a')
            outext = outext//'---->File '//self%wavesInputFile(i)%name
        end do
    end if
    if (allocated(self%waterPropsInputFile)) then
        if (written) outext = outext//new_line('a')
        outext = outext//'--->'//Globals%DataTypes%waterProps%startcase()//' data '
        do i=1, size(self%waterPropsInputFile)
            outext = outext//new_line('a')
            outext = outext//'---->File '//self%waterPropsInputFile(i)%name
        end do
    end if
    if (.not.self%useInputFiles) outext = '-->Input streamer stack is empty, no input data'
    call Log%put(outext,.false.)
    end subroutine printInputStreamer

    end module simulationInputStreamer_mod