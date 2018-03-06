    program MOHIDLagrangian

    use source
    use CLA    !Command line argument module
    use commom_modules

    implicit none

    !-----------------------------------------------------------------------------------------------------------------------------------
    character(len=STRLEN)  :: defxmlfilename_char
    character(len=STRLEN)  :: outpath_char
    type(string) :: defxmlfilename
    type(string) :: outpath
    character(4) :: xmlextention = '.xml'
    !-----------------------------------------------------------------------------------------------------------------------------------

    ! Initialize command line arguments
    call cla_init

    ! Register keys for key/value pairs.
    call cla_register('-i','--infile','input definition file (xml)', cla_char, 'defxmlfilename')
    call cla_register('-o','--outpath','output path', cla_char, 'outpath')

    ! Store the arguments with each respective variable
    call cla_get('-i',defxmlfilename_char)
    call cla_get('-o',outpath_char)
    ! CLA uses chars, we use 'strings'
    defxmlfilename=trim(defxmlfilename_char)
    outpath=trim(outpath_char)

    ! Completing the names with file extensions
    defxmlfilename=defxmlfilename//xmlextention
    outpath=outpath//'\'

    ! Initialization routines to build the xml-defined case
    call init_mohidlagrangian(defxmlfilename)


    end program MOHIDLagrangian