    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_vtkwritter
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : July 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines a vtk writer class with an object exposable to the Output streamer
    !------------------------------------------------------------------------------

    module vtkwritter_mod

    use common_modules
    use vtk_fortran
    use boundingbox_mod
    use blocks_mod

    implicit none
    private
    
    type :: vtkwritter_class
        integer :: numVtkFiles
    contains
    procedure :: initialize => initVTKwritter
    procedure :: Domain
    procedure :: TracerSerial
    end type vtkwritter_class
    
    type(vtkwritter_class) :: vtkWritter

    !Public access vars
    public :: vtkWritter

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes a VTK writer object
    !---------------------------------------------------------------------------
    subroutine initVTKwritter(self)
    implicit none
    class(vtkwritter_class), intent(inout) :: self
    self%numVtkFiles = 0
    end subroutine initVTKwritter
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public Tracer writting routine. Writes Tracer data in binary XML VTK
    !> format using an unstructured grid. Serial writer for serial files.
    !> @parm[in] self, filename, blocks
    !---------------------------------------------------------------------------
    subroutine TracerSerial(self, filename, blocks)
    implicit none
    class(vtkwritter_class), intent(inout) :: self
    type(string), intent(in) :: filename
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    
    type(vtk_file) :: vtkfile
    type(string) :: fullfilename
    type(string) :: outext
    integer :: error, i
    integer :: np
    integer, parameter :: nc = 0           !< Number of cells
    integer(I1P), dimension(1:nc) :: cell_type  !< Cells type
    integer(I4P), dimension(1:nc) :: offset     !< Cells offset
    integer(I4P), dimension(:), allocatable :: connect    !< Connectivity

    fullfilename = filename%chars()//'.vtu'
    outext = '->Writting output file '//fullfilename
    call Log%put(outext)
    fullfilename = Globals%Names%outpath//'/'//fullfilename
    
    error = vtkfile%initialize(format='binary', filename=fullfilename%chars(), mesh_topology='UnstructuredGrid')
    !Write the data of each block
    do i = 1, size(blocks)
        if (blocks(i)%Tracer%numActive > 0) then
            np = blocks(i)%Tracer%numActive
            allocate(connect(np))
            error = vtkfile%xml_writer%write_piece(np=np, nc=nc)
            error = vtkfile%xml_writer%write_geo(np=np, nc=nc, x=blocks(i)%AoT%x, y=blocks(i)%AoT%y, z=blocks(i)%AoT%z)
            error = vtkfile%xml_writer%write_connectivity(nc=nc, connectivity=connect, offset=offset, cell_type=cell_type)
            error = vtkfile%xml_writer%write_dataarray(location='node', action='open')
            error = vtkfile%xml_writer%write_dataarray(data_name='id', x=blocks(i)%AoT%id)
            error = vtkfile%xml_writer%write_dataarray(data_name='velocity', x=blocks(i)%AoT%u, y=blocks(i)%AoT%v, z=blocks(i)%AoT%w)
            error = vtkfile%xml_writer%write_dataarray(location='node', action='close')
            error = vtkfile%xml_writer%write_piece()
            deallocate(connect)
        end if
    end do
    error = vtkfile%finalize()    
    self%numVtkFiles = self%numVtkFiles + 1
    
    end subroutine TracerSerial
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation domain writting routine. Writes binary XML VTK
    !> format using an unstructured grid.
    !> @parm[in] self, filename, bbox, npbbox, blocks
    !---------------------------------------------------------------------------
    subroutine Domain(self, filename, bbox, npbbox, blocks)
    implicit none
    class(vtkwritter_class), intent(inout) :: self
    type(string), intent(in) :: filename                    !< name of the case to add
    class(boundingbox_class), intent(in) :: bbox            !< Case bounding box
    class(block_class), dimension(:), intent(in) :: blocks  !< Case Blocks
    integer, intent(in) :: npbbox                           !< number of points of the bbox geometry
    
    type(vtk_file) :: vtkfile                              !< file object
    type(string) :: fullfilename
    type(string) :: outext
    integer :: error, i, b
    integer, parameter :: nc = 1                            !< Number of cells of the geometry
    real(prec), dimension(1:npbbox) :: xx, yy, zz           !< coordinates of the  geometry
    type(vector) :: pts(npbbox)                             !< coordinates of the  geometry
    integer, dimension(1:npbbox) :: connect, var            !< Connectivity array and a simple pointwise variable
    integer(I4P), dimension(1:nc) :: offset
    integer(I1P), dimension(1:nc) :: cell_type              !< VTK cell type. 1 for point and 12 for hexahedron
    
    offset = [8]
    cell_type = [12]
    
    !preparing file
    fullfilename = filename%chars()//'_BoundingBox.vtu'
    outext = '->Writting Bounding Box file '//fullfilename
    call Log%put(outext)
    fullfilename = Globals%Names%outpath//'/'//fullfilename
    
    error = vtkfile%initialize(format='binary', filename=fullfilename%chars(), mesh_topology='UnstructuredGrid')
    
    !Writting bounding box geometry
    pts = Geometry%getPoints(bbox)
    do i=1, npbbox
        xx(i) = pts(i)%x
        yy(i) = pts(i)%y
        zz(i) = pts(i)%z
        connect(i) = i-1
    end do    
    error = vtkfile%xml_writer%write_piece(np=npbbox, nc=nc)
    error = vtkfile%xml_writer%write_geo(np=npbbox, nc=nc, x=xx, y=yy, z=zz)
    error = vtkfile%xml_writer%write_connectivity(nc=nc, connectivity=connect, offset=offset, cell_type=cell_type)
    error = vtkfile%xml_writer%write_piece()
    
    !Closing file
    error = vtkfile%finalize()
    
    !preparing file
    fullfilename = filename%chars()//'_Blocks.vtu'
    outext = '->Writting Blocks file '//fullfilename
    call Log%put(outext)
    fullfilename = Globals%Names%outpath//'/'//fullfilename
    
    error = vtkfile%initialize(format='binary', filename=fullfilename%chars(), mesh_topology='UnstructuredGrid')
    
    !Writting block geometries
    do b=1, size(blocks)
        pts = Geometry%getPoints(blocks(b)%extents)
        do i=1, npbbox
            xx(i) = pts(i)%x
            yy(i) = pts(i)%y
            !zz(i) = pts(i)%z
            connect(i) = i-1
            var(i) = b
        end do    
        error = vtkfile%xml_writer%write_piece(np=npbbox, nc=nc)
        error = vtkfile%xml_writer%write_geo(np=npbbox, nc=nc, x=xx, y=yy, z=zz)
        error = vtkfile%xml_writer%write_connectivity(nc=nc, connectivity=connect, offset=offset, cell_type=cell_type)
        error = vtkfile%xml_writer%write_dataarray(location='node', action='open')
        error = vtkfile%xml_writer%write_dataarray(data_name='Block', x=var)
        error = vtkfile%xml_writer%write_dataarray(location='node', action='close')
        error = vtkfile%xml_writer%write_piece()
    end do
        
    !Closing file
    error = vtkfile%finalize()
    
    end subroutine Domain
    
    

  end module vtkwritter_mod
