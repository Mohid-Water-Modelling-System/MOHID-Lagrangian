
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : MTimeSeriesParser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module defining a MOHID time series parsing class and methods, encapsulates the
    !> MOHID ModuleTimeSerie module, from Base 1.
    !------------------------------------------------------------------------------

    module MTimeSeriesParser_mod

    use ModuleTimeSerie
    use ModuleGlobalData
    use ModuleTime 
    use stringifor
    use datetime_module
    
    use simulationLogger_mod
    use simulationPrecision_mod
    use utilities_mod
    use simulationGlobals_mod

    implicit none
    private

    type :: mTimeSeriesParser_class  !< The .csv parser class
        type(string), dimension(:), allocatable :: varName
        real(prec), dimension(:,:), allocatable :: dataMatrix
    contains
    procedure :: getFile
    end type mTimeSeriesParser_class

    !Public access vars
    public :: mTimeSeriesParser_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a csv file and returns the file with the data.
    !> Gets and trashes the header.
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, filename)
    class(mTimeSeriesParser_class), intent(inout) :: self
    type(string), intent(in) :: filename
    integer :: id, STAT_CALL
    integer :: nColumns, nVals
    type(string) :: outext
    type(string) :: header, temp
    type(datetime) :: initialDateTime
    type(timedelta) :: dateTimeOffset
    real(prec) :: timeOffset
    character(len=line_length) :: mHeader
    real, dimension(:,:), pointer :: mDataMatrix
    type(T_Time) :: mInitialDate
    real :: mDate(6)
    
    outext = '-> Reading MOHID Time Series file '// filename
    call Log%put(outext)
    id = 0
    call StartTimeSerieInput(id, filename%chars(), STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = '[MTimeSeriesParser::getFile]: Cannot open '//filename//' file, stoping'
        call Log%put(outext)
        stop
    end if
    
    !getting data limits and header
    call GetTimeSerieDataColumns(id, nColumns, STAT = STAT_CALL)    
    call GetTimeSerieDataValues(id, nVals, STAT = STAT_CALL)
    call GetTimeSerieHeader(id, mHeader, STAT = STAT_CALL)
    !getting actual data matrix, already with time in seconds from date in file
    allocate(mDataMatrix(nVals, nColumns))
    call GetTimeSerieDataMatrix(id, mDataMatrix, STAT = STAT_CALL)
    !getting initial date and offseting 
    call GetTimeSerieInitialData(id, mInitialDate, STAT = STAT_CALL)
    call ExtractDate(mInitialDate, mDate(1), mDate(2), mDate(3), mDate(4), mDate(5), mDate(6))
    initialDateTime = Utils%getDateTimeFromDate(int(mDate))
    dateTimeOffset = initialDateTime - Globals%SimTime%StartDate
    timeOffset = dateTimeOffset%total_seconds()
    print*, mDataMatrix(1,:)
    mDataMatrix(:,1) = mDataMatrix(:,1) + timeOffset
    print*, mDataMatrix(1,:)
    !storing data in our object
    allocate(self%dataMatrix(nVals, nColumns))
    self%dataMatrix = mDataMatrix
    !cleaning up the header and getting individual data labels
    header = mHeader    
    temp = header%replace(old='!', new='')
    temp = temp%replace(old='  ', new=' ')
    do while (temp/= header)
        header = temp
        temp = header%replace(old='  ', new=' ')
    end do
    call header%split(tokens=self%varName, sep=' ')
    
    !destroying the MOHID Time Series module instance
    call KillTimeSerie(id, STAT = STAT_CALL)
    
    end subroutine getFile

    end module MTimeSeriesParser_mod
