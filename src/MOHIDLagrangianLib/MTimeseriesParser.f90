
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
        
    use common_modules

    implicit none
    private

    type :: mTimeSeriesParser_class  !< The .csv parser class
        private
        type(string) :: filename
        type(string), dimension(:), allocatable :: varName
        real(prec), dimension(:,:), allocatable :: dataMatrix
    contains
    procedure :: getFile
    procedure :: getDataByLabel
    end type mTimeSeriesParser_class

    !Public access vars
    public :: mTimeSeriesParser_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that parses a MOHID time series file and stores the data internaly.
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, filename, varList)
    class(mTimeSeriesParser_class), intent(inout) :: self
    type(string), intent(in) :: filename
    type(string), dimension(:), optional, intent(in) :: varList
    integer :: id, STAT_CALL, i, j
    integer :: nColumns, nVals
    type(string) :: outext
    type(string) :: header, temp
    type(datetime) :: initialDateTime
    type(timedelta) :: dateTimeOffset
    real(prec) :: timeOffset
    logical :: found
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
    !getting initial date and offseting to match simulation time dimension
    call GetTimeSerieInitialData(id, mInitialDate, STAT = STAT_CALL)
    call ExtractDate(mInitialDate, mDate(1), mDate(2), mDate(3), mDate(4), mDate(5), mDate(6))
    initialDateTime = Utils%getDateTimeFromDate(int(mDate))
    dateTimeOffset = initialDateTime - Globals%SimTime%StartDate
    timeOffset = dateTimeOffset%total_seconds()
    mDataMatrix(:,1) = mDataMatrix(:,1) + timeOffset !correcting for current time dimension
    !storing data in parser object
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
    !spliting and storing in parser object
    call header%split(tokens=self%varName, sep=' ')
    !searching and replacing the variable names with required simulation names
    if (present(varList)) then
        do i=1, size(varList)
            found = .false.
            do j=1,size(self%varName)
                if (Globals%Var%getVarSimName(self%varName(j)) == varList(i)) then
                    self%varName(j) = Globals%Var%getVarSimName(self%varName(j))
                    found = .true.
                    exit
                end if
            end do
            if (.not.found) then
                outext = '[MTimeSeriesParser::getFile]: File '//filename//' doesnt list variable representing '//varList(i)//', stoping'
                call Log%put(outext)
                stop
            end if
        end do
    end if
    self%filename = filename
    
    !destroying the MOHID Time Series module instance
    call KillTimeSerie(id, STAT = STAT_CALL)
    
    end subroutine getFile
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that gets data from a column of a MOHID time series file, given 
    !> a string with a variable name
    !> @param[in] self, dataName, data_out
    !---------------------------------------------------------------------------
    subroutine getDataByLabel(self, dataName, data_out)
    class(mTimeSeriesParser_class), intent(in) :: self
    type(string), intent(in) :: dataName
    real(prec), dimension(:), allocatable :: data_out
    integer :: idx
    type(string) :: outext
    
    idx = Utils%find_str(self%varName, dataName)
    if (idx == MV_INT) then
        outext = '[MTimeSeriesParser::getColumn]: File '//self%filename//' doesnt list variable representing '//dataName//', stoping'
        call Log%put(outext)
        stop
    end if
    allocate(data_out(size(self%dataMatrix,2)))
    data_out = self%dataMatrix(:,idx)
    
    end subroutine getDataByLabel
    

    end module MTimeSeriesParser_mod