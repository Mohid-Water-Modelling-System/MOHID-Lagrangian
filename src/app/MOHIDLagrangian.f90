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

    use CLA    !Command line argument module
    use commom_modules
    use sources_mod
    use simulation_mod

    implicit none

    !-----------------------------------------------------------------------------------------------------------------------------------
    character(len=STRLEN)  :: casefilename_char
    character(len=STRLEN)  :: outpath_char
    type(string) :: casefilename !< Simulation input case file
    type(string) :: outpath      !< Simulation output path
    
    type(simulation_class) :: Sim !< Simulation object
    
    !-----------------------------------------------------------------------------------------------------------------------------------

    ! Initialize command line arguments
    call cla_init
    ! Register keys for key/value pairs.
    call cla_register('-i','--infile','input definition file (xml)', cla_char, 'casefilename')
    call cla_register('-o','--outpath','output path', cla_char, 'outpath')
    ! Store the arguments with each respective variable
    call cla_get('-i',casefilename_char)
    call cla_get('-o',outpath_char)
    casefilename=trim(casefilename_char)
    outpath=trim(outpath_char)
    ! Completing the outpath with the separator
    outpath=outpath//'/'
    
    !Our Lagrangian adventure starts here. Strap on.
    call Sim%initialize(casefilename,outpath)

    !main time cycle
    do while (Globals%SimTime .LT. Globals%Parameters%TimeMax)



        !Do your Lagrangian things here :D


    Globals%SimTime = Globals%SimTime + Globals%SimDefs%dt
    enddo

    ! Finalization of the program - deallocation, file closing, etc
    !call finalizeMohidLagrangian
    call Sim%finalize()
    
    end program MOHIDLagrangian
