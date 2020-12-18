    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : about
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : Feb 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to print version, licence, preambles.
    !------------------------------------------------------------------------------

    module simulationAbout_mod

    use common_modules

    implicit none
    private

    !Public access procedures
    public :: PrintLicPreamble

    !version control
    type(string) :: version
    type(string) :: author
    type(string) :: date

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public licence and preamble printer routine.
    !---------------------------------------------------------------------------
    subroutine PrintLicPreamble
    implicit none
    type(string) :: outext

    version  ="v0.8"
    author   ="R. Birjukovs Canelas - D. Garaboa Paz"
    date     ="17-11-2020"

    outext = ' __  __  ___  _   _ ___ ____  _                                      _              '//new_line('a')//&
        ' |  \/  |/ _ \| | | |_ _|  _ \| |    __ _  __ _ _ __ __ _ _ __   __ _(_) __ _ _ __  '//new_line('a')//&
        ' | |\/| | | | | |_| || || | | | |   / _` |/ _` | __/  _` |  _ \ / _` | |/ _` |  _ \ '//new_line('a')//&
        ' | |  | | |_| |  _  || || |_| | |__| (_| | (_| | | | (_| | | | | (_| | | (_| | | | |'//new_line('a')//&
        ' |_|  |_|\___/|_| |_|___|____/|_____\__,_|\__, |_|  \__,_|_| |_|\__, |_|\__,_|_| |_|'//new_line('a')//&
        '                                          |___/                 |___/               '//new_line('a')//&

        '  <MOHIDLagrangian> Copyright (C) 2018 by'//new_line('a')//&
        '  R. Birjukovs Canelas'//new_line('a')//&
        '  MARETEC - Research Centre for Marine, Environment and Technology'//new_line('a')//&
        '  University of Lisbon - IST'//new_line('a')//&
        '  A. Daniel Garaboa Paz'//new_line('a')//&
        '  Non-Linear Physics Group - University of Santiago de Compostela, Spain'//new_line('a')//&
        ''//new_line('a')//&
        '  MOHID Lagrangian is free software: you can redistribute it and/or'//new_line('a')//&
        '  modify it under the terms of the GNU General Public License as'//new_line('a')//&
        '  published by the Free Software Foundation, either version 3 of'//new_line('a')//&
        '  the License, or (at your option) any later version.'//new_line('a')//&
        ''//new_line('a')//&
        '  MOHID Lagrangian is distributed WITHOUT ANY WARRANTY; without even'//new_line('a')//&
        '  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR'//new_line('a')//&
        '  PURPOSE. See the GNU General Public License for more details.'//new_line('a')//&
        ''//new_line('a')//&
        '  You should have received a copy of the GNU General Public License,'//new_line('a')//&
        '  along with MOHID Lagrangian. If not, see <http://www.gnu.org/licenses/>.,'//new_line('a')//&
        ''//new_line('a')//&
        ''//new_line('a')//&
        'MOHID Lagrangian '//version//' ('//author//') ('//date//')'//new_line('a')//&
        '====================================================================='

    call Log%put(outext,.false.)

    end subroutine

    end module simulationAbout_mod
