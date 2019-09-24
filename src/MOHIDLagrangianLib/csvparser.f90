
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : csvparser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module defining a csv parsing and emitting class and methods, encapsulates the
    !> fortran-csv-module library. https://github.com/jacobwilliams/fortran-csv-module
    !------------------------------------------------------------------------------

    module csvParser_mod

    use csv_module
    use stringifor

    use common_modules
    use simulationLogger_mod
    use simulationPrecision_mod

    implicit none
    private

    type :: csvparser_class  !< The .csv parser class
        private
        type(string) :: filename
        type(string), dimension(:), allocatable :: varName
        real(prec), dimension(:,:), allocatable :: dataMatrix
    contains
    procedure :: getFileToFileObject
    procedure :: getFile
    procedure :: getDataByLabel
    procedure :: getColumn
    procedure :: closeFile
    end type csvparser_class

    !Public access vars
    public :: csvparser_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a csv file and returns the file with the data.
    !> Gets and trashes the header.
    !> @param[in] self, csvfilename, header_row
    !---------------------------------------------------------------------------
    subroutine getFile(self, csvfilename, header_row)
    class(csvparser_class), intent(inout) :: self
    type(string), intent(in) :: csvfilename
    integer, intent(in) :: header_row
    character(len=30), dimension(:), allocatable :: header_csv
    type(string) :: header, temp
    type(csv_file) :: csvdoc
    logical :: status
    integer :: nColumns, nVals, i
    real(prec), dimension(:), allocatable :: test_data
    type(string) :: outext

    outext = '-> Reading .csv file '// csvfilename
    call Log%put(outext)
    call csvdoc%read(csvfilename%chars(), header_row=header_row, status_ok=status)
    if (status .eqv. .false.) then
        outext = '[csvparser_class::getFile]: Cannot open .csv file, supposedly at '// csvfilename //', stoping'
        call Log%put(outext)
        stop
    end if
    ! get the header
    call csvdoc%get_header(header_csv,status)

    !cleaning up the header and getting individual data labels
    do i=1, size(header_csv)
        header = header // header_csv(i)
    end do
    temp = header%replace(old='!', new='')
    temp = temp%replace(old='  ', new=' ')
    do while (temp/= header)
        header = temp
        temp = header%replace(old='  ', new=' ')
    end do
    !spliting and storing in parser object
    call header%split(tokens=self%varName, sep=' ')
    nColumns = size(self%varName)
    call csvdoc%get(1, test_data, status)
    if (status .eqv. .false.) then
        outext = 1
        outext = '[csvparser_class::getFile]: Cannot get data from column'// outext //' .csv file '// csvfilename //', stoping'
        call Log%put(outext)
        stop
    end if
    nVals = size(test_data)
    allocate(self%dataMatrix(nVals, nColumns))
    do i=1, size(self%dataMatrix,2)
        call csvdoc%get(i, test_data, status)
        if (status .eqv. .false.) then
            outext = i
            outext = '[csvparser_class::getFile]: Cannot get data from column'// outext //' .csv file '// csvfilename //', stoping'
            call Log%put(outext)
            stop
        end if
        self%dataMatrix(:,i) = test_data
    end do
    self%filename = csvfilename

    ! print*, self%filename%chars()
    ! do i=1, size(self%varName)
    !     print*, self%varName(i)%chars()
    ! end do
    ! print*, self%dataMatrix

    end subroutine getFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that gets data from a column of a csV file, given
    !> a string with a variable name
    !> @param[in] self, dataName, data_out
    !---------------------------------------------------------------------------
    subroutine getDataByLabel(self, dataName, data_out)
    class(csvparser_class), intent(in) :: self
    type(string), intent(in) :: dataName
    real(prec), dimension(:), allocatable, intent(out) :: data_out
    integer :: idx
    type(string) :: outext

    idx = Utils%find_str(self%varName, dataName)
    if (idx == MV_INT) then
        outext = '[csvparser_class::getDataByLabel]: File '//self%filename//' doesn''t list variable representing '//dataName//', stoping'
        call Log%put(outext)
        stop
    end if
    allocate(data_out(size(self%dataMatrix,1)))
    data_out = self%dataMatrix(:,idx)

    end subroutine getDataByLabel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a csv file and returns the file with the data.
    !> Gets and trashes the header.
    !> @param[in] self, csvdoc, csvfilename, header_row
    !---------------------------------------------------------------------------
    subroutine getFileToFileObject(self, csvdoc, csvfilename, header_row)
    class(csvparser_class), intent(in) :: self
    type(csv_file), intent(inout) :: csvdoc
    type(string), intent(in) :: csvfilename
    integer, intent(in) :: header_row
    character(len=30), dimension(:), allocatable :: header
    logical :: status
    type(string) :: outext
    outext = '-> Reading .csv file '// csvfilename
    call Log%put(outext)
    call csvdoc%read(csvfilename%chars(), header_row=header_row, status_ok=status)
    if (status .eqv. .false.) then
        outext = '[csvparser_class::getFile]: Cannot open .csv file, supposedly at '// csvfilename //', stoping'
        call Log%put(outext)
        stop
    end if
    ! get the header
    call csvdoc%get_header(header,status)
    end subroutine getFileToFileObject

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that gets data from a column of a .csv file.
    !> @param[in] self, csvdoc, csvfilename, column, data_out
    !---------------------------------------------------------------------------
    subroutine getColumn(self, csvdoc, csvfilename, column, data_out)
    class(csvparser_class), intent(in) :: self
    type(csv_file), intent(inout) :: csvdoc
    type(string), intent(in) :: csvfilename
    integer, intent(in) :: column
    real(prec), dimension(:), allocatable :: data_out
    logical :: status
    type(string) :: outext
    call csvdoc%get(column, data_out, status)
    if (status .eqv. .false.) then
        outext = column
        outext = '[csvparser_class::getFile]: Cannot get data from column'// outext //' .csv file '// csvfilename //', stoping'
        call Log%put(outext)
        stop
    end if
    end subroutine getColumn

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that closes a .csv file.
    !> @param[in] self, csvdoc
    !---------------------------------------------------------------------------
    subroutine closeFile(self, csvdoc)
    class(csvparser_class), intent(in) :: self
    type(csv_file), intent(inout) :: csvdoc
    call csvdoc%destroy()
    end subroutine closeFile


    end module csvParser_mod
