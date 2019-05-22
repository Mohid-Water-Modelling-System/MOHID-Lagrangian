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
    end type inputFileModel_class

    type :: input_streamer_class        !< Input Streamer class
        logical :: useInputFiles
        type(inputFileModel_class), allocatable, dimension(:) :: currentsInputFile !< array of input file metadata for currents
        type(inputFileModel_class), allocatable, dimension(:) :: windsInputFile !< array of input file metadata for currents
        type(inputFileModel_class), allocatable, dimension(:) :: wavesInputFile !< array of input file metadata for currents
        real(prec) :: buffer_size                                               !< half of the biggest tail of data behind current time
    contains
    procedure :: initialize => initInputStreamer
    procedure :: loadDataFromStack
    procedure :: getFullFile

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
    class(boundingbox_class), intent(in) :: bBox            !< Case bounding box
    class(block_class), dimension(:), intent(inout) :: blocks  !< Case Blocks
    integer :: i
    integer :: fNumber

    if (self%useInputFiles) then

        do i=1, size(self%currentsInputFile)
            if (self%currentsInputFile(i)%endTime > Globals%SimTime%CurrTime) then
                if (self%currentsInputFile(i)%startTime <= Globals%SimTime%CurrTime) then
                    fNumber = i
                    exit
                end if
            end if
        end do

        do i=1, size(blocks)
            if (Globals%Sim%getnumdt() == 1 ) then
                allocate(blocks(i)%Background(1))
                blocks(i)%Background(1) = self%getFullFile(fNumber)
            end if
        end do

    end if

    end subroutine loadDataFromStack

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a NC file
    !> @param[in] self, nfile
    !---------------------------------------------------------------------------
    type(background_class) function getFullFile(self, nfile)
    class(input_streamer_class), intent(in) :: self
    integer, intent(in) :: nfile
    type(ncfile_class) :: ncFile
    type(scalar1d_field_class), allocatable, dimension(:) :: backgrounDims
    type(generic_field_class) :: gfield1, gfield2, gfield3
    type(string) :: name
    type(box) :: extents
    type(vector) :: pt
    real(prec), dimension(3,2) :: dimExtents
    integer :: i

    call ncFile%initialize(self%currentsInputFile(nfile)%name)
    call ncFile%getVarDimensions(Globals%Var%u, backgrounDims)
    call ncFile%getVar(Globals%Var%u, gfield1)
    call ncFile%getVar(Globals%Var%v, gfield2)
    call ncFile%getVar(Globals%Var%w, gfield3)
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

    name = self%currentsInputFile(nfile)%name%basename(strip_last_extension=.true.)
    getFullFile = Background(nfile, name, extents, backgrounDims)
    call getFullFile%add(gfield1)
    call getFullFile%add(gfield2)
    call getFullFile%add(gfield3)

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
    type(string) :: tag, att_name, att_val, outext
    type(string), allocatable, dimension(:) :: fileNames
    integer :: i

    self%buffer_size = 3600*24*2 !seconds/hour*hours*days

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