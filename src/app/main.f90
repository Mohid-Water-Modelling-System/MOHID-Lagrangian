    program MOHIDLagrangian

    use tracer
    use PENF
    use CLA    !Command line argument module

    !-----------------------------------------------------------------------------------------------------------------------------------
    implicit none
    character(len=:), allocatable :: parsedxml      !< String containing the parsed XML data
    type(xml_file)                :: xfile          !< XML file handler
    integer                       :: xunit          !< XML file unit
    !-----------------------------------------------------------------------------------------------------------------------------------
    character(len=STRLEN)  :: defxmlfilename
    character(len=STRLEN)  :: outdefxmlfilename
    character(4) :: xmlextention = '.xml'
    character(len=STRLEN)  :: outpath
    !-----------------------------------------------------------------------------------------------------------------------------------


    ! Initialize command line arguments
    call cla_init

    ! Register your keys for key/value pairs.
    ! EACH KEY MUST BEGIN WITH -
    call cla_register('-i','--infile','input definition file (xml)', cla_char, 'defxmlfilename')
    call cla_register('-o','--outpath','output path', cla_char, 'outpath')

    ! Store the arguments with each respective variable
    call cla_get('-i',defxmlfilename)
    call cla_get('-o',outpath)
    !print *,' --infile     = ',trim(defxmlfilename)
    !print *,' --outpath    = ',trim(outpath)

    ! Completing the names with file extensions
    defxmlfilename=trim(defxmlfilename)//xmlextention
    outpath=trim(outpath)//'\'
    print "(A)",' defxmlfilename     = ',defxmlfilename
    print "(A)",' outpath     = ',outpath

    !Reading input XML file
    call xfile%parse(filename=defxmlfilename)
    parsedxml = xfile%stringify()
    !print "(A)", parsedxml
    if (parsedxml == '') then !file is empty, might as well stop now
        print *, 'Definition file', defxmlfilename, 'is empty. Stoping run.'
        stop
    endif

    !writting the xml file to the output path
    outdefxmlfilename=outpath//defxmlfilename
    outdefxmlfilename=trim(outdefxmlfilename)
    print "(A)",' outdefxmlfilename     = ',outdefxmlfilename
    open(newunit=xunit, file=outdefxmlfilename, access='STREAM', form='UNFORMATTED')
    write(unit=xunit)parsedxml
    close(unit=xunit)

    end program MOHIDLagrangian