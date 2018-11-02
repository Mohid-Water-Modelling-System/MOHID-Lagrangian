
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

    module csvparser_mod

    use csv_module
    use common_modules

    implicit none
    private

    type :: csvparser_class  !< The .xml parser class
    contains
    procedure :: getFile
    procedure :: getColumn
    procedure :: closeFile
    end type csvparser_class

    !Simulation variables
    type(csvparser_class) :: CSVReader

    !Public access vars
    public :: CSVReader

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a csv file and returns the file with the data.
    !> Gets and trashes the header.
    !> @param[in] self, csvdoc, csvfilename, header_row
    !---------------------------------------------------------------------------
    subroutine getFile(self, csvdoc, csvfilename, header_row)
    class(csvparser_class), intent(in) :: self
    type(csv_file), intent(inout) :: csvdoc
    type(string), intent(in) :: csvfilename
    integer, intent(in) :: header_row
    character(len=30), dimension(:), allocatable :: header
    logical :: status
    type(string) :: outext
    call csvdoc%read(csvfilename%chars(), header_row=header_row, status_ok=status)
    if (status .eqv. .false.) then
        outext = '[CSVReader::getFile]: Cannot open .csv file, supposedly at '// csvfilename //', stoping'
        call Log%put(outext)
        stop
    end if
    ! get the header
    call csvdoc%get_header(header,status)
    end subroutine getFile

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
        outext = '[CSVReader::getFile]: Cannot get data from column'// outext //' .csv file '// csvfilename //', stoping'
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


    end module csvparser_mod
