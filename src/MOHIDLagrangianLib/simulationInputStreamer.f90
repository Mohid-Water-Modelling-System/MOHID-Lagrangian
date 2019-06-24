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
    end type inputFileModel_class

    type :: input_streamer_class        !< Input Streamer class
        logical :: useInputFiles
        type(inputFileModel_class), allocatable, dimension(:) :: currentsInputFile !< array of input file metadata for currents
        type(inputFileModel_class), allocatable, dimension(:) :: windsInputFile !< array of input file metadata for currents
        type(inputFileModel_class), allocatable, dimension(:) :: wavesInputFile !< array of input file metadata for currents
        real(prec) :: bufferSize                                               !< half of the biggest tail of data behind current time
        real(prec) :: lastReadTime
    contains
    procedure :: initialize => initInputStreamer
    procedure :: loadDataFromStack
    procedure :: getFullFile
    procedure, private :: resetReadStatus
    procedure :: print => printInputStreamer
    end type input_streamer_class

    !Public access vars
    public :: input_streamer_class

    contains

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
    type(background_class) :: tempBkgd
    integer :: i, j
    integer :: fNumber
    real(prec) :: tempTime(2)
    logical :: needToRead, appended

    needToRead = .false.
    if (self%useInputFiles) then
        !check if we need to import data (current time and buffer size)
        if (self%lastReadTime <= Globals%SimTime%CurrTime + self%BufferSize/4.0) needToRead = .true.
        if (self%lastReadTime >= Globals%SimTime%TimeMax) needToRead = .false.
        if (needToRead) then
            call self%resetReadStatus()
            !check what files on the stack are to read to backgrounds
            do i=1, size(self%currentsInputFile)
                if (self%currentsInputFile(i)%endTime >= Globals%SimTime%CurrTime) then
                    if (self%currentsInputFile(i)%startTime <= Globals%SimTime%CurrTime + self%BufferSize) then
                        if (.not.self%currentsInputFile(i)%used) self%currentsInputFile(i)%toRead = .true.
                    end if
                end if
            end do
            !read selected files
            do i=1, size(self%currentsInputFile)
                if (self%currentsInputFile(i)%toRead) then
                    !import data to temporary background
                    tempBkgd = self%getFullFile(self%currentsInputFile(i)%name)
                    self%currentsInputFile(i)%used = .true.
                    do j=1, size(blocks)
                        if (.not.allocated(blocks(j)%Background)) allocate(blocks(j)%Background(1))
                        !slice data by block and either join to existing background or add a new one
                        if (blocks(j)%Background(1)%initialized) call blocks(j)%Background(1)%append(tempBkgd%getHyperSlab(blocks(j)%extents), appended)
                        if (.not.blocks(j)%Background(1)%initialized) blocks(j)%Background(1) = tempBkgd%getHyperSlab(blocks(j)%extents)

                        !save last time already loaded
                        tempTime = blocks(j)%Background(1)%getDimExtents(Globals%Var%time)
                        self%lastReadTime = tempTime(2)

                    end do
                    !clean out the temporary background data (this structure, even tough it is a local variable, has pointers inside)
                    call tempBkgd%finalize()
                end if
            end do
        end if
    end if

    end subroutine loadDataFromStack

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a NC file
    !> @param[in] self, nfile
    !---------------------------------------------------------------------------
    type(background_class) function getFullFile(self, fileName)
    class(input_streamer_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(ncfile_class) :: ncFile
    type(scalar1d_field_class), allocatable, dimension(:) :: backgrounDims
    type(generic_field_class), allocatable, dimension(:) :: gfield
    type(string) :: name
    type(box) :: extents
    type(vector) :: pt
    real(prec), dimension(3,2) :: dimExtents
    integer :: i
    type(string) :: outext

    allocate(gfield(4))

    outext = '->Reading '//fileName
    call Log%put(outext,.false.)

    call ncFile%initialize(fileName)
    call ncFile%getVarDimensions(Globals%Var%u, backgrounDims)
    call ncFile%getVar(Globals%Var%u, gfield(1))
    call ncFile%getVar(Globals%Var%v, gfield(2))
    call ncFile%getVar(Globals%Var%w, gfield(3))
    !reading a field to use later as land mask
    call ncFile%getVar(Globals%Var%u, gfield(4), .true.)
    gfield(4)%name = Globals%Var%landMask
    call ncFile%finalize()

    dimExtents = 0.0
    do i = 1, size(backgrounDims)
        if (backgrounDims(i)%name == Globals%Var%lon) then
            dimExtents(1,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(1,2) = backgrounDims(i)%getFieldMaxBound()
        else if (backgrounDims(i)%name == Globals%Var%lat) then
            dimExtents(2,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(2,2) = backgrounDims(i)%getFieldMaxBound()
        else if (backgrounDims(i)%name == Globals%Var%level) then
            dimExtents(3,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(3,2) = backgrounDims(i)%getFieldMaxBound()
        end if
    end do
    extents%pt = dimExtents(1,1)*ex + dimExtents(2,1)*ey + dimExtents(3,1)*ez
    pt = dimExtents(1,2)*ex + dimExtents(2,2)*ey + dimExtents(3,2)*ez
    extents%size = pt - extents%pt

    name = fileName%basename(strip_last_extension=.true.)
    getFullFile = Background(1, name, extents, backgrounDims)
    do i =1, size(gfield)
        call getFullFile%add(gfield(i))
        call gfield(i)%finalize()
    end do

    !call getFullFile%print()

    end function getFullFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the input writer object, imports metadata on input files
    !---------------------------------------------------------------------------
    subroutine initInputStreamer(self)
    class(input_streamer_class), intent(inout) :: self
    type(Node), pointer :: xmlInputs           !< .xml file handle
    type(Node), pointer :: fileNode
    type(NodeList), pointer :: fileList
    type(string) :: tag, att_name, att_val
    type(string), allocatable, dimension(:) :: fileNames
    integer :: i

    self%bufferSize = Globals%Parameters%BufferSize
    self%lastReadTime = -1.0

    call XMLReader%getFile(xmlInputs,Globals%Names%inputsXmlFilename, mandatory = .false.)
    if (associated(xmlInputs)) then
        self%useInputFiles = .true.
        !Go to the file_collection node
        tag = "file_collection"
        call XMLReader%gotoNode(xmlInputs,xmlInputs,tag)
        tag = "currents"
        call XMLReader%gotoNode(xmlInputs,xmlInputs,tag)
        fileList => getElementsByTagname(xmlInputs, "file")       !searching for tags with the 'namingfile' name
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
            self%currentsInputFile(i+1)%startTime = att_val%to_number(kind=1._R4P)
            tag="endTime"
            att_name="value"
            call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
            self%currentsInputFile(i+1)%endTime = att_val%to_number(kind=1._R4P)
            self%currentsInputFile(i+1)%used = .false.
        end do
        call Globals%setInputFileNames(fileNames)
    else
        self%useInputFiles = .false.
    end if
    end subroutine initInputStreamer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> resets input files read status
    !---------------------------------------------------------------------------
    subroutine resetReadStatus(self)
    class(input_streamer_class), intent(inout) :: self
    integer :: i
    do i=1, size(self%currentsInputFile)
        self%currentsInputFile(i)%toRead = .false.
    end do
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
    outext = '-->Input streamer stack:'//new_line('a')
    outext = outext//'--->Currents data '
    do i=1, size(self%currentsInputFile)
        outext = outext//new_line('a')
        outext = outext//'---->File '//self%currentsInputFile(i)%name!//new_line('a')
        !temp_str=self%currentsInputFile(i)%startTime
        !outext = outext//'      Starting time is '//temp_str//' s'//new_line('a')
        !temp_str=self%currentsInputFile(i)%endTime
        !outext = outext//'      Ending time is   '//temp_str//' s'
    end do
    if (.not.self%useInputFiles) outext = '-->Input streamer stack is empty, no input data'
    call Log%put(outext,.false.)
    end subroutine printInputStreamer

    end module simulationInputStreamer_mod