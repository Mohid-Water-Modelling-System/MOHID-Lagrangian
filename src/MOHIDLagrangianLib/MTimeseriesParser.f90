
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
    use stringifor
    
    use simulationLogger_mod
    use simulationPrecision_mod

    implicit none
    private

    type :: mTimeSeriesParser_class  !< The .csv parser class
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
    class(mTimeSeriesParser_class), intent(in) :: self
    type(string), intent(in) :: filename
    integer :: id, STAT_CALL
    type(string) :: outext
    
    outext = '-> Reading MOHID Time Series file '// filename
    call Log%put(outext)
    id = 0
    call StartTimeSerieInput(id, filename%chars(), STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = '[MTimeSeriesParser::getFile]: Cannot open '//filename//' file, stoping'
        call Log%put(outext)
        stop
    end if
    
    end subroutine getFile

    end module MTimeSeriesParser_mod
