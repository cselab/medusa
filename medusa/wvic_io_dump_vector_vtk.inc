!-------------------------------------------------------------------------------
!* filename: wvic_io_dump_vector_vtk                                        *!
!* project : ppm                                                              *!
!* purpose : dump a vector field to a vtk file                                *!
!*         :                                                                  *!
!* author  : Johannes Tophoj Rasmussen                                        *!
!*         : Department of Mechanical Engineering, Fluid Mechanics            *!
!*         : DTU                                                              *!
!*         : jtra@mek.dtu.dk                                                  *!
!*                                                                            *!
!* Output precision is changed (and hardcoded) through the                    *!
!* precision_vtk parameter in wvic_module_io                                  *!
!------------------------------------------------------------------------------!


SUBROUTINE wvic_io_dump_vector_vtk(vectfield, n_data, vectname, iter,info)

  USE module_wvic
  USE ppm_module_write
  USE ppm_module_data
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! arguments
  !-----------------------------------------------------------------------------
  CHARACTER(len=256), INTENT(IN)             :: vectname 
  REAL(mk), DIMENSION(:,:,:,:,:),    POINTER :: vectfield
  INTEGER, DIMENSION(:,:), POINTER           :: n_data 
  INTEGER, INTENT(in)                        :: iter 
  INTEGER, INTENT(inout)                     :: info

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  CHARACTER(len=256)       :: vtkfile, buffer,datatype,msg
  CHARACTER                :: lf
  INTEGER                  :: i,j,k
  INTEGER                  :: imin,jmin,kmin,imax,jmax,kmax
  INTEGER                  :: iming,jming,kming,imaxg,jmaxg,kmaxg
  INTEGER                  :: ioff0, nbytes0
  INTEGER                  :: isub, isubl
  INTEGER                  :: ivectname,idatatype

  isub = 1;
  isubl = isublist(isub)

  lf = char(10)

  !-----------------------------------------------------------------------------
  ! Local bounds in accordance with the ImageData format (vti)
  !-----------------------------------------------------------------------------
  imin = istart(1,isubl)-1
  jmin = istart(2,isubl)-1
  kmin = istart(3,isubl)-1
  ! For the upper bounds n_data(...)-2 should ensure that exactly all the data
  ! is saved, but paraview seems unable to connect points between topologies.
  ! Therefore one more slice is included and n_data(...)-1 is used.
  imax = istart(1,isubl)-1+n_data(1,isubl)-1
  jmax = istart(2,isubl)-1+n_data(2,isubl)-1
  kmax = istart(3,isubl)-1+n_data(3,isubl)-1

  !nsubs    = ndims(1)*ndims(2)*ndims(3)

  !-----------------------------------------------------------------------------
  ! Calculate data size (and offset in case of multiple fields in the file),
  ! Precision is set trough the precision_vtk parameter
  !-----------------------------------------------------------------------------
  IF (precision_vtk .EQ. 1) THEN
    datatype  = 'Float8'
    idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 2) THEN
    datatype  = 'Float16'
    idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 4) THEN
    datatype  = 'Float32'
    idatatype = LEN_TRIM(datatype)
  ELSEIF (precision_vtk .EQ. 8) THEN
    datatype  = 'Float64'
    idatatype = LEN_TRIM(datatype)
  END IF
  !PPM replace n_data
  nbytes0   = 3*n_data(1,isubl)*n_data(2,isubl)*n_data(3,isubl)*precision_vtk
  ioff0     = 0

  !-----------------------------------------------------------------------------
  ! Ensure that a proper handle has been specified for the field
  !-----------------------------------------------------------------------------
  ivectname = LEN_TRIM(vectname)
  IF (ivectname.EQ.0) THEN
    WRITE(msg,*) 'Error: Name for field must be specified'
    CALL ppm_write(rank,'wvic_dump_vector_vtk',msg,info)
    STOP
  END IF
  
  !-----------------------------------------------------------------------------
  ! If RANK 0 output the parallel pvti file
  !-----------------------------------------------------------------------------
  IF (rank .EQ. 0) THEN
    !---------------------------------------------------------------------------
    ! Global bounds in accordance with the ImageData format (vti)
    !---------------------------------------------------------------------------
    iming = 0
    jming = 0
    kming = 0
    imaxg = nx(1)-1
    jmaxg = nx(2)-1
    kmaxg = nx(3)-1

    WRITE(vtkfile,'(A,A,A,A,I5.5,A)')&
      & vectname(1:ivectname),'_',runtag(1:iruntag),'I',iter,'.pvti'

    OPEN(17,file=vtkfile,access='stream')

    buffer =  '<?xml version="1.0"?>'//lf
    WRITE(17) TRIM(buffer)
    buffer =  '<VTKFile type="PImageData" version="0.1" byte_order="LittleEndian">'//lf
    WRITE(17) TRIM(buffer)
    WRITE(buffer,'(A,A,6(I8),A,3(E20.12),A,3(E20.12),A)') '  <PImageData',&
      &  ' WholeExtent="',iming,imaxg,jming,jmaxg,kming,kmaxg, &
      & '" Ghostlevel="0" Origin="',min_physg(1),min_physg(2),min_physg(3), &
      & '" Spacing="',dx,dy,dz,'">'//lf
    WRITE(17) TRIM(buffer)
    buffer =  '    <PPointData Vectors="'//vectname(1:ivectname)//'">'//lf
    WRITE(17) TRIM(buffer)
    buffer =  '      <PDataArray type="Float32" '//&
           &'Name="'//vectname(1:ivectname)//'" NumberOfComponents="3"/>'//lf
    WRITE(17) TRIM(buffer)
    buffer =  '    </PPointData>'//lf
    WRITE(17) TRIM(buffer)

    !PPM replace n_data
    ! Loop over all topologies
    DO i=1,ppm_nsubs(ppm_topoid)
      WRITE(buffer,'(A,A,6(I8),A,A,A,A,A,I5.5,A,I5.5,A)') '    <Piece ', &
        & 'Extent="',istart(1,i)-1,istart(1,i)+n_data(1,i)-2,istart(2,i)-1, &
        & istart(2,i)+n_data(2,i)-2,istart(3,i)-1,istart(3,i)+n_data(3,i)-2, &
        & '" Source="', vectname(1:ivectname),'_',runtag(1:iruntag),'R', &
        & sub2proc(i),'I',iter, '.vti"/>'//lf
      WRITE(17) TRIM(buffer)
    END DO

    buffer =  '  </PImageData>'//lf
    WRITE(17) TRIM(buffer)
    buffer =  '</VTKFile>'//lf
    WRITE(17) TRIM(buffer)
    CLOSE(17)

  END IF
  ! End of outputting the parallel pvti file
  !-----------------------------------------------------------------------------

  WRITE(vtkfile,'(A,A,A,A,I5.5,A,I5.5,A)')&
    & vectname(1:ivectname),'_',runtag(1:iruntag),'R',rank,'I',iter,'.vti'

  OPEN(17,file=vtkfile,access='stream')

  !-----------------------------------------------------------------------------
  ! Write header
  !-----------------------------------------------------------------------------
  buffer =  '<?xml version="1.0"?>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
  WRITE(17) TRIM(buffer)
  WRITE(buffer,'(A,A,6(I8),A,3(E20.12),A,3(E20.12),A)') '  <ImageData',&
    &  ' WholeExtent="',imin,imax,jmin,jmax,kmin,kmax, &
    & '" Origin="',min_physg(1),min_physg(2),min_physg(3), &
    & '" Spacing="',dx,dy,dz,'">\n'
  WRITE(17) TRIM(buffer)
  WRITE(buffer,'(A,A,6(I8),A,3(I8),A)') '    <Piece ',&
    &  'Extent="',imin,imax,jmin,jmax,kmin,kmax,'">\n'
  WRITE(17) TRIM(buffer)
  buffer =  '      <PointData Vectors="'//vectname(1:ivectname)//'">'//lf
  WRITE(17) TRIM(buffer)
  WRITE(buffer,'(A,A,A,A,A,A,I8,A)') '        <DataArray type="', &
   & datatype(1:idatatype),'" Name="',vectname(1:ivectname), &
   & '" NumberOfComponents="3" Format="appended" ', &
!   & 'RangeMin="" RangeMax="" ', & !JTR has been omitted for now
   & 'offset="',ioff0,'" />\n'
  WRITE(17) TRIM(buffer)
  buffer =  '      </PointData>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '      <CellData>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '      </CellData>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '    </Piece>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '  </ImageData>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '  <AppendedData encoding="raw">'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '_'
  WRITE(17) TRIM(buffer)

  !-----------------------------------------------------------------------------
  ! Write field - finish XML file, and close
  !-----------------------------------------------------------------------------
  WRITE(17) nbytes0
  DO k=1,n_data(3,isubl)
    DO j=1,n_data(2,isubl)
      DO i=1,n_data(1,isubl)
        WRITE(17) REAL(vectfield(1,i,j,k,isub),precision_vtk)
        WRITE(17) REAL(vectfield(2,i,j,k,isub),precision_vtk)
        WRITE(17) REAL(vectfield(3,i,j,k,isub),precision_vtk)
      END DO
    END DO
  END DO

  buffer =  lf//'  </AppendedData>'//lf
  WRITE(17) TRIM(buffer)
  buffer =  '</VTKFile>'//lf
  WRITE(17) TRIM(buffer)
  CLOSE(17)


9999 CONTINUE
  RETURN
END SUBROUTINE wvic_io_dump_vector_vtk

