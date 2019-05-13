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

    use FoX_dom

    implicit none
    private

    type :: inputFileModel_class !< Input file model class
        type(string) :: name        !< name of the file
        real(prec) :: startTime     !< starting time of the data on the file
        real(prec) :: endTime       !< ending time of the data on the file
    end type inputFileModel_class

    type :: input_streamer_class        !< Input Streamer class
        type(inputFileModel_class), allocatable, dimension(:) :: inputFileModel !< array of input file metadata
        real(prec) :: buffer_size                                               !< half of the biggest tail of data behind current time
    contains
    procedure :: initialize => initInputStreamer
    procedure :: print => printInputStreamer
    end type input_streamer_class

    !Public access vars
    public :: input_streamer_class

    contains

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

    call XMLReader%getFile(xmlInputs,Globals%Names%inputsXmlFilename)
    !Go to the file_collection node
    tag = "file_collection"
    call XMLReader%gotoNode(xmlInputs,xmlInputs,tag)
    fileList => getElementsByTagname(xmlInputs, "file")       !searching for tags with the 'namingfile' name
    allocate(fileNames(getLength(fileList)))
    allocate(self%inputFileModel(getLength(fileList)))
    do i = 0, getLength(fileList) - 1
        fileNode => item(fileList, i)
        tag="name"
        att_name="value"
        call XMLReader%getNodeAttribute(fileNode, tag, att_name, fileNames(i+1))
        self%inputFileModel(i+1)%name = fileNames(i+1)
        tag="startTime"
        att_name="value"
        call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
        self%inputFileModel(i+1)%startTime = att_val%to_number(kind=1._R4P)
        tag="endTime"
        att_name="value"
        call XMLReader%getNodeAttribute(fileNode, tag, att_name, att_val)
        self%inputFileModel(i+1)%endTime = att_val%to_number(kind=1._R4P)
    end do
    call Globals%setInputFileNames(fileNames)
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
    outext = '-->Input Streamer stack:'
    do i=1, size(self%inputFileModel)
        outext = outext//new_line('a')
        outext = outext//'--->File '//self%inputFileModel(i)%name//new_line('a')
        temp_str=self%inputFileModel(i)%startTime
        outext = outext//'      Starting time is '//temp_str//' s'//new_line('a')
        temp_str=self%inputFileModel(i)%endTime
        outext = outext//'      Ending time is   '//temp_str//' s'
    end do
    call Log%put(outext,.false.)
    end subroutine printInputStreamer

    end module simulationInputStreamer_mod