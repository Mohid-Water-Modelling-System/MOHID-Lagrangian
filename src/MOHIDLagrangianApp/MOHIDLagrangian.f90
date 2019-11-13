    !MARETEC - Research Centre for Marine, Environment and Technology
    !Copyright (C) 2018  Ricardo Birjukovs Canelas
    !
    !This program is free software: you can redistribute it and/or modify
    !it under the terms of the GNU General Public License as published by
    !the Free Software Foundation, either version 3 of the License, or
    !(at your option) any later version.
    !
    !This program is distributed in the hope that it will be useful,
    !but WITHOUT ANY WARRANTY; without even the implied warranty of
    !MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !GNU General Public License for more details.
    !
    !You should have received a copy of the GNU General Public License
    !along with this program.  If not, see <https://www.gnu.org/licenses/>.

    program MOHIDLagrangian

        use flap !Command line argument lib   
        use common_modules
        use simulation_mod
    
        implicit none
    
        !-----------------------------------------------------------------------------------------------------------------------------------
        type(command_line_interface) :: cli
        integer :: ierr
        
        character(len=CHAR_LEN)  :: casefilename_char
        character(len=CHAR_LEN)  :: outpath_char
        type(string) :: casefilename !< Simulation input case file
        type(string) :: outpath      !< Simulation output path
    
        type(simulation_class) :: Sim !< Simulation object
    
        !-----------------------------------------------------------------------------------------------------------------------------------
    
        ! Initialize command line arguments
        call cli%init(description = 'MOHID Lagrangian command line argument parser')
    
        ! Register keys for key/value pairs.
        call cli%add( switch = '--infile', switch_ab = '-i', &
            help = 'input definition file (.xml)', &
            required = .true., act = 'store', error = ierr )
        if (ierr/=0) stop
        call cli%add( switch = '--outpath', switch_ab = '-o', &
            help = 'output path', required = .true., &
            act = 'store', error = ierr )
        if (ierr/=0) stop
        call cli%get( switch = '-i', val = casefilename_char, error = ierr )
        if (ierr/=0) stop
        call cli%get( switch = '-o', val = outpath_char, error = ierr )
        if (ierr/=0) stop
        ! Store the arguments with each respective variable    
        casefilename=trim(casefilename_char)
        outpath=trim(outpath_char)
        ! Completing the outpath with the separator
        outpath=outpath//'/'
    
        !Our Lagrangian adventure starts here. Strap on. In case of emergency landing put your head between your knees.
        call Sim%initialize(casefilename,outpath)
    
        !Main time cycle
        call Sim%run()
    
        ! Finalization of the program - deallocation, file closing, etc
        call Sim%finalize()
    
        end program MOHIDLagrangian