module Output_HDF5_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use Output_Aux_module
  use Output_Common_module
  
  use PFLOTRAN_Constants_module
  
  implicit none

  private
  
  PetscMPIInt, private, parameter :: ON=1, OFF=0

  public :: OutputHDF5GetH5Filename, &
            OutputHDF5GetGroupName_Time, &
            OutputHDF5FilenameID, &
            OutputHDF5DatasetStringArray, &
            OutputHDF5AttributeStringArray, &
            OutputHDF5OpenFile, &
            OutputHDF5CloseFile, &
            OutputHDF5WriteStructCoord, &
            WriteHDF5FluxVelocities, &
            DetermineNumVertices, &
            OutputHDF5WriteSnapShotAtts, &
            OutputHDF5Provenance, &
            OutputHDF5OpenGroup

contains

! ************************************************************************** !

subroutine OutputHDF5GetH5Filename(option, output_option, h5file, &
                                   var_list_type, filename)
  !
  ! Determine the propper hdf5 output file name and open it.
  !
  ! Return the file handle and 'first' flag indicating if this is the
  ! first time the file has been opened.
  !
  use Option_module
  use Utility_module, only : Equal

  implicit none

  type(option_type), intent(inout) :: option
  type(output_option_type), intent(in) :: output_option
  type(output_hdf5_type) :: h5file
  PetscInt, intent(in) :: var_list_type
  character(len=MAXSTRINGLENGTH) :: filename

  PetscBool :: first
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string,string2,string3

  filename = 'Uninitialized.h5'
  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') &
        int(option%time/output_option%periodic_snap_output_time_incr)
  end select

  if (output_option%print_single_h5_file) then
    first = h5file%first_write
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               trim(string2) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number, &
                output_option%times_per_h5_file) == 0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (Equal(mod((option%time-output_option%periodic_snap_output_time_incr)/ &
             output_option%periodic_snap_output_time_incr, &
             dble(output_option%times_per_h5_file)),0.d0)) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-' // trim(string) // trim(string2) // '.h5'
  endif
  h5file%first_write = first
  
end subroutine OutputHDF5GetH5Filename

! ************************************************************************** !

subroutine OutputHDF5OpenFile(option, h5file, filename, file_id)
  !
  ! Opens an HDF5 file, creating it if the first time
  !
  ! Author: Glenn Hammond
  ! Date: 10/18/19
  !
  use Option_module

  implicit none

  type(option_type) :: option
  type(output_hdf5_type) :: h5file
  character(len=MAXSTRINGLENGTH) :: filename
  integer(HID_T), intent(out) :: file_id

  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err

    ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  if (.not.h5file%first_write) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) h5file%first_write = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (h5file%first_write) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                      H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (h5file%first_write) then
    option%io_buffer = '--> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // trim(filename)
  endif
  call PrintMsg(option)

end subroutine OutputHDF5OpenFile

! ************************************************************************** !

subroutine OutputHDF5CloseFile(option, h5file, file_id)

  use Option_module

  implicit none

  type(option_type), intent(in) :: option
  type(output_hdf5_type) :: h5file
  integer(HID_T), intent(in) :: file_id
  
  integer :: hdf5_err
  PetscErrorCode :: ierr

  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)
  h5file%first_write = PETSC_FALSE

end subroutine OutputHDF5CloseFile

! ************************************************************************** !

function OutputHDF5FilenameID(output_option,option,var_list_type)
  ! 
  ! This subroutine creates an ID for HDF5 filename for:
  ! - Instantaneous, or
  ! - Temporally averaged variables.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/10/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  PetscInt :: var_list_type

  character(len=MAXWORDLENGTH) :: OutputHDF5FilenameID
  PetscInt :: file_number

  select case(var_list_type)
    case (INSTANTANEOUS_VARS)
      file_number = floor(real(output_option%plot_number)/ &
                               output_option%times_per_h5_file)
    case (AVERAGED_VARS)
      file_number = floor((option%time - &
                           output_option%periodic_snap_output_time_incr)/ &
                           output_option%periodic_snap_output_time_incr/ &
                           output_option%times_per_h5_file)
  end select

  if (file_number < 10) then
    write(OutputHDF5FilenameID,'("00",i1)') file_number
  else if (output_option%plot_number < 100) then
    write(OutputHDF5FilenameID,'("0",i2)') file_number  
  else if (output_option%plot_number < 1000) then
    write(OutputHDF5FilenameID,'(i3)') file_number  
  else if (output_option%plot_number < 10000) then
    write(OutputHDF5FilenameID,'(i4)') file_number
  else if (output_option%plot_number < 100000) then
    write(OutputHDF5FilenameID,'(i5)') file_number
  else
    option%io_buffer = 'Plot number exceeds current maximum of 10^5.'
    call PrintErrMsgToDev(option,'ask for a higher maximum')
  endif 
  
  OutputHDF5FilenameID = adjustl(OutputHDF5FilenameID)

end function OutputHDF5FilenameID

! ************************************************************************** !

subroutine OutputHDF5WriteStructCoord(realization_base,file_id)
  ! 
  ! Write structured grid coordinates to HDF5 file.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/18/19
  ! 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Patch_module

  use HDF5_module
  use HDF5_Aux_module
  
  implicit none

  class(realization_base_type) :: realization_base
  integer(HID_T) :: file_id

  integer(HID_T) :: grp_id
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscReal, allocatable :: array(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscMPIInt :: hdf5_err

  option => realization_base%option
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid

  ! create a group for the coordinates data set
  string = "Coordinates"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

  !GEH - Structured Grid Dependence - Begin
  ! write out coordinates in x, y, and z directions
  string = "X [m]"
  allocate(array(grid%structured_grid%nx+1))
  array(1) = discretization%origin_global(X_DIRECTION)
  do i=2,grid%structured_grid%nx+1
    array(i) = array(i-1) + grid%structured_grid%dx_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%nx+1, &
                            array,grp_id)
  deallocate(array)

  string = "Y [m]"
  allocate(array(grid%structured_grid%ny+1))
  array(1) = discretization%origin_global(Y_DIRECTION)
  do i=2,grid%structured_grid%ny+1
    array(i) = array(i-1) + grid%structured_grid%dy_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%ny+1, &
                            array,grp_id)
  deallocate(array)

  string = "Z [m]"
  allocate(array(grid%structured_grid%nz+1))
  array(1) = discretization%origin_global(Z_DIRECTION)
  do i=2,grid%structured_grid%nz+1
    array(i) = array(i-1) + grid%structured_grid%dz_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%nz+1, &
                            array,grp_id)
  deallocate(array)
  !GEH - Structured Grid Dependence - End

  call h5gclose_f(grp_id,hdf5_err)

end subroutine OutputHDF5WriteStructCoord

! ************************************************************************** !

function OutputHDF5GetGroupName_Time(option,output_option)
  ! 
  ! Returns a the group name, which is the formatted time (in the simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/18/19
  ! 
  use Option_module
  
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXSTRINGLENGTH) :: OutputHDF5GetGroupName_Time

  character(len=MAXSTRINGLENGTH) :: string
  
    ! create a group for the data set
  write(string,'(''Time:'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  
  OutputHDF5GetGroupName_Time = string
  
end function OutputHDF5GetGroupName_Time

! ************************************************************************** !

subroutine OutputHDF5OpenGroup(file_id,grp_id,string)
  ! 
  ! Write structured grid coordinates to HDF5 file.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/18/19
  ! 
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscMPIInt :: hdf5_err
  
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  
end subroutine OutputHDF5OpenGroup

! ************************************************************************** !

subroutine WriteHDF5Coordinates(name,option,length,array,file_id)
  ! 
  ! Writes structured coordinates to HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 
  use Option_module
  use Logging_module
  
  implicit none
  
  character(len=32) :: name
  type(option_type) :: option
  PetscInt :: length
  PetscReal :: array(:)
  integer(HID_T) :: file_id
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  PetscMPIInt :: rank
  integer :: hdf5_err
  PetscMPIInt :: hdf5_flag
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_output_coordinates_hdf5, &
                          ierr);CHKERRQ(ierr)

  ! write out grid structure
  rank = 1
  dims = 0
  ! x-direction
  dims(1) = length
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err) ! must be independent and only from p0
#endif
  if (option%myrank == option%io_rank) then
     call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
     ! this is due to a bug in hdf5-1.8.18 hwere H5S_ALL_F is an INTEGER.  It
     ! should be INTEGER(HID_T)
     memory_space_id = H5S_ALL_F
     call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
     call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_output_coordinates_hdf5, &
                        ierr);CHKERRQ(ierr)

end subroutine WriteHDF5Coordinates

! ************************************************************************** !

subroutine WriteHDF5CoordinatesUGrid(grid,option,file_id)
  ! 
  ! This subroutine writes unstructured coordinates to HDF5 file
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/31/12
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use HDF5_module, only : trick_hdf5
  use Variables_module
  use Logging_module
  
  implicit none

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err  

  PetscInt :: local_size
  PetscInt :: istart
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element
  PetscErrorCode :: ierr

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3
  
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(double_array(local_size*3))
  
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo
  
  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)


  call VecDestroy(global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_vertex_vec,ierr);CHKERRQ(ierr)


  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax
   
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size*NINE_INTEGER
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%nmax
  dims(1) = NINE_INTEGER
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = NINE_INTEGER
  
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(int_array(local_size*NINE_INTEGER))
  
  do i=1,local_size
    int_array((i-1)*9 + 1) = 0
    int_array((i-1)*9 + 2) = INT(vec_ptr((i-1)*8+1))
    int_array((i-1)*9 + 3) = INT(vec_ptr((i-1)*8+2))
    int_array((i-1)*9 + 4) = INT(vec_ptr((i-1)*8+3))
    int_array((i-1)*9 + 5) = INT(vec_ptr((i-1)*8+4))
    int_array((i-1)*9 + 6) = INT(vec_ptr((i-1)*8+5))
    int_array((i-1)*9 + 7) = INT(vec_ptr((i-1)*8+6))
    int_array((i-1)*9 + 8) = INT(vec_ptr((i-1)*8+7))
    int_array((i-1)*9 + 9) = INT(vec_ptr((i-1)*8+8))
    do j=2,9
      if (int_array((i-1)*9 + j)>0) int_array((i-1)*9 + 1)= int_array((i-1)*9 + 1) +1
    enddo
  enddo
  
  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteHDF5CoordinatesUGrid

! ************************************************************************** !

subroutine WriteHDF5CoordinatesUGridXDMF(realization_base,option,file_id)
  ! 
  ! This routine writes unstructured coordinates to HDF5 file in XDMF format
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/2012
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module
  use Logging_module
  
  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err  

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscInt :: temp_int, istart
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  Vec :: global_x_cell_vec,global_y_cell_vec,global_z_cell_vec
  Vec :: natural_x_cell_vec,natural_y_cell_vec,natural_z_cell_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element, ugdm_cell
  PetscErrorCode :: ierr

  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => realization_base%patch%grid

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)
  
  ! create a group for the coordinates data set
  string = "Domain"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on 
  !     Juqueen
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(double_array(local_size*3))
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)


  call VecDestroy(global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_vertex_vec,ierr);CHKERRQ(ierr)

  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax

  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on 
  !     Juqueen
  call MPI_Allreduce(vert_count,temp_int,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr)
  dims(1) = temp_int
  realization_base%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(vert_count, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (4) ! Tetrahedron
        int_array(vert_count) = TET_ID_XDMF
      case (5) ! Pyramid
        int_array(vert_count) = PYR_ID_XDMF
      case (6) ! Wedge
        int_array(vert_count) = WED_ID_XDMF
      case (8) ! Hexahedron
        int_array(vert_count) = HEX_ID_XDMF
    end select

    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*8+j))-1
      endif
    enddo
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

  ! Cell center X/Y/Z
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_z_cell_vec,ierr);CHKERRQ(ierr)

  call OutputGetCellCoordinates(grid, global_x_cell_vec,X_COORDINATE)
  call OutputGetCellCoordinates(grid, global_y_cell_vec,Y_COORDINATE)
  call OutputGetCellCoordinates(grid, global_z_cell_vec,Z_COORDINATE)


  call UGridCreateUGDM(grid%unstructured_grid,ugdm_cell,ONE_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_x_cell_vec,NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_y_cell_vec,NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_z_cell_vec,NATURAL,option)
                           
  call VecScatterBegin(ugdm_cell%scatter_gton,global_x_cell_vec, &
                       natural_x_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_x_cell_vec, &
                     natural_x_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_y_cell_vec, &
                       natural_y_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_y_cell_vec, &
                     natural_y_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_z_cell_vec, &
                       natural_z_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_z_cell_vec, &
                     natural_z_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecGetArrayF90(natural_x_cell_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_y_cell_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_z_cell_vec,vec_z_ptr,ierr);CHKERRQ(ierr)
  local_size = grid%unstructured_grid%nlmax

  ! XC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "XC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_x_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  ! YC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "YC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_y_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  ! ZC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "ZC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_z_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  
  call h5gclose_f(grp_id,hdf5_err)

  call VecRestoreArrayF90(natural_x_cell_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(natural_y_cell_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(natural_z_cell_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(global_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_cell_vec,ierr);CHKERRQ(ierr)

  call VecDestroy(natural_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_z_cell_vec,ierr);CHKERRQ(ierr)

  call UGridDMDestroy(ugdm_cell)

end subroutine WriteHDF5CoordinatesUGridXDMF

! ************************************************************************** !

subroutine DetermineNumVertices(realization_base,option)
  ! 
  ! Determine the number of vertices written out in the output HDF5 file
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/13/2015
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  PetscInt :: local_size,vert_count
  PetscInt :: i
  PetscInt :: temp_int

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  type(ugdm_type),pointer :: ugdm_element
  PetscErrorCode :: ierr

  grid => realization_base%patch%grid

  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER, &
                       option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax

  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call MPI_Allreduce(vert_count,temp_int,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr)
  realization_base%output_option%xmf_vert_len=temp_int

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine DetermineNumVertices

! ************************************************************************** !

subroutine WriteHDF5FluxVelocities(name,realization_base,iphase,direction, &
                                   file_id,h5file)
  ! 
  ! Print flux velocities to HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  use hdf5
  use HDF5_module, only : HDF5WriteStructuredDataSet, trick_hdf5

  implicit none

  character(len=32) :: name
  class(realization_base_type) :: realization_base
  PetscInt :: iphase
  PetscInt :: direction
  integer(HID_T) :: file_id
  type(output_hdf5_type) :: h5file

  PetscInt :: i, j, k
  PetscInt :: count, iconn
  PetscInt :: local_id, ghosted_id
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
    
  PetscReal, allocatable :: array(:)

  PetscBool, save :: trick_flux_vel_x = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_y = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_z = PETSC_FALSE

  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option  

  ! in a few cases (i.e. for small test problems), some processors may
  ! have no velocities to print.  This results in zero-length arrays
  ! in collective H5Dwrite().  To avoid, we switch to independent
  ! H5Dwrite() and don't write from the zero-length procs. 
!GEH - Structured Grid Dependence - Begin
  if (h5file%first_write) then
    trick_flux_vel_x = PETSC_FALSE
    trick_flux_vel_y = PETSC_FALSE
    trick_flux_vel_z = PETSC_FALSE
    
    nx_local = grid%structured_grid%nlx
    ny_local = grid%structured_grid%nly
    nz_local = grid%structured_grid%nlz
    if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (i == 0) trick_flux_vel_x = PETSC_TRUE
    if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (j == 0) trick_flux_vel_y = PETSC_TRUE
    if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (k == 0) trick_flux_vel_z = PETSC_TRUE
  endif

  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz

  select case(direction)
    case(X_DIRECTION)
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = PETSC_TRUE
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = PETSC_TRUE
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
        nz_local = grid%structured_grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = PETSC_TRUE
  end select  

  allocate(array(nx_local*ny_local*nz_local))
  call OutputCollectVelocityOrFlux(realization_base, iphase, direction, &
                                   PETSC_FALSE, array)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * output_option%tconv

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                                  option,nx_global,ny_global,nz_global, &
                                  nx_local,ny_local,nz_local, &
                                  grid%structured_grid%lxs, &
                                  grid%structured_grid%lys, &
                                  grid%structured_grid%lzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = PETSC_FALSE

end subroutine WriteHDF5FluxVelocities

! ************************************************************************** !

subroutine OutputHDF5Provenance(option, output_option, file_id)
  !
  ! write pflotran and petsc provenance information including a copy
  ! of the inputfile
  !

  use Option_module, only : option_type
  use Output_Aux_module, only : output_option_type
  use PFLOTRAN_Provenance_module, only : provenance_max_str_len

  implicit none

  type(option_type), intent(in) :: option
  type(output_option_type), intent(in) :: output_option
  integer(HID_T), intent(in) :: file_id

  character(len=32) :: filename, name
  integer(HID_T) :: prop_id, provenance_id, string_type
  PetscMPIInt :: hdf5_err
  PetscBool :: first
  integer(SIZE_T) :: size_t_int

  ! create the provenance group
  name = "Provenance"
  call h5gcreate_f(file_id, name, provenance_id, hdf5_err, &
                   OBJECT_NAMELEN_DEFAULT_F)

  ! create fixed length string datatype
  call h5tcopy_f(H5T_FORTRAN_S1, string_type, hdf5_err)
  size_t_int = provenance_max_str_len
  call h5tset_size_f(string_type, size_t_int, hdf5_err)

  call OutputHDF5Provenance_PFLOTRAN(option, provenance_id, string_type)
  call OutputHDF5Provenance_PETSc(provenance_id, string_type)

  ! close the provenance group
  call h5tclose_f(string_type, hdf5_err)
  call h5gclose_f(provenance_id, hdf5_err)

end subroutine OutputHDF5Provenance

! ************************************************************************** !

subroutine OutputHDF5Provenance_PFLOTRAN(option, provenance_id, string_type)
  !
  ! write the pflotran provenance data as attributes (small) or
  ! datasets (big details)
  !

  use Option_module, only : option_type
  use PFLOTRAN_Provenance_module

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type

  character(len=32) :: name
  integer(HID_T) :: pflotran_id
  PetscMPIInt :: hdf5_err

  ! Create the pflotran group under provenance
  name = "PFLOTRAN"
  call h5gcreate_f(provenance_id, name, pflotran_id, hdf5_err, &
                   OBJECT_NAMELEN_DEFAULT_F)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_date_time", &
                                    ONE_INTEGER, pflotran_compile_date_time)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_user", &
                                    ONE_INTEGER, pflotran_compile_user)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_hostname", &
                                    ONE_INTEGER, pflotran_compile_hostname)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, &
                                      "pflotran_status", &
                                      ONE_INTEGER, pflotran_status)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, &
                                      "pflotran_changeset", &
                                      ONE_INTEGER, pflotran_changeset)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_fflags", &
                                    detail_pflotran_fflags_len, &
                                    detail_pflotran_fflags)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_status", &
                                    detail_pflotran_status_len, &
                                    detail_pflotran_status)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_parent", &
                                    detail_pflotran_parent_len, &
                                    detail_pflotran_parent)

  ! FIXME(bja, 2013-11-25): break gcc when diffs are present  
  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_diff", &
                                    detail_pflotran_diff_len, &
                                    detail_pflotran_diff)

  call OutputHDF5Provenance_input(option, pflotran_id)

  ! close pflotran group
  call h5gclose_f(pflotran_id, hdf5_err)

end subroutine OutputHDF5Provenance_PFLOTRAN

! ************************************************************************** !

subroutine OutputHDF5Provenance_input(option, pflotran_id)
  !
  ! open the pflotran input file, figure out how long it is, read it
  ! into a buffer, then write the buffer as a pflotran provenance
  ! group dataset.
  !
  use Input_Aux_module, only : input_type, InputCreate, InputDestroy, &
       InputGetLineCount, InputReadToBuffer
  use Option_module, only : option_type
  use PFLOTRAN_Constants_module, only : IN_UNIT, MAXSTRINGLENGTH

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: pflotran_id

  integer(HID_T) :: input_string_type
  type(input_type), pointer :: input
  PetscInt :: i, input_line_count
  character(len=MAXSTRINGLENGTH), allocatable :: input_buffer(:)
  PetscMPIInt :: hdf5_err
  integer(SIZE_T) :: size_t_int

  input => InputCreate(IN_UNIT, option%input_filename, option)
  input_line_count = InputGetLineCount(input,option)
  allocate(input_buffer(input_line_count))
  call InputReadToBuffer(input, input_buffer, option)
  call h5tcopy_f(H5T_FORTRAN_S1, input_string_type, hdf5_err)
  size_t_int = MAXWORDLENGTH
  call h5tset_size_f(input_string_type, size_t_int, hdf5_err)
  call OutputHDF5DatasetStringArray(pflotran_id, input_string_type, &
                                    "pflotran_input_file", &
                                    input_line_count, input_buffer)
  call h5tclose_f(input_string_type, hdf5_err)
  deallocate(input_buffer)
  call InputDestroy(input)

end subroutine OutputHDF5Provenance_input

! ************************************************************************** !

subroutine OutputHDF5Provenance_PETSc(provenance_id, string_type)
  !
  ! write the petsc provenance data as attributes (small) or datasets
  ! (big details)
  !

  use PFLOTRAN_Provenance_module

  implicit none

  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type

  character(len=32) :: name
  integer(HID_T) :: petsc_id
  PetscMPIInt :: hdf5_err

  ! create the petsc group under provenance
  name = "PETSc"
  call h5gcreate_f(provenance_id, name, petsc_id, hdf5_err, &
                   OBJECT_NAMELEN_DEFAULT_F)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, "petsc_status", &
                                      ONE_INTEGER, petsc_status)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, &
                                      "petsc_changeset", &
                                      ONE_INTEGER, petsc_changeset)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_status", &
                                    detail_petsc_status_len, &
                                    detail_petsc_status)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_parent", &
                                    detail_petsc_parent_len, &
                                    detail_petsc_parent)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_config", &
                                    detail_petsc_config_len, &
                                    detail_petsc_config)

  ! close the petsc group
  call h5gclose_f(petsc_id, hdf5_err)

end subroutine OutputHDF5Provenance_PETSc

! ************************************************************************** !

subroutine OutputHDF5AttributeStringArray(parent_id, type, name, length, data)
  ! create the dataspaces and attributes consisting of an array of
  ! strings, then write the data and cleanup

  implicit none

  integer(HID_T), intent(in) ::  parent_id, type
  character(len=*), intent(in) :: name
  PetscInt, intent(in) :: length
  character(len=*), intent(in) :: data(length)

  integer(HID_T) :: dataspace_id, attribute_id
  integer(HSIZE_T), dimension(1:1) :: dims
  PetscMPIInt :: hdf5_err

  dims = length
  call h5screate_simple_f(1, dims, dataspace_id, hdf5_err)
  call h5acreate_f(parent_id, name, type, dataspace_id, attribute_id, hdf5_err)
  call h5awrite_f(attribute_id, type, data, dims, hdf5_err)
  call h5aclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5AttributeStringArray

! ************************************************************************** !

subroutine OutputHDF5DatasetStringArray(parent_id, type, name, length, data)
  ! create the dataspaces and dataset consisting of an array of
  ! strings, then write the data and cleanup

  implicit none

  integer(HID_T), intent(in) ::  parent_id, type
  character(len=*), intent(in) :: name
  PetscInt, intent(in) :: length
  character(len=*), intent(in) :: data(length)

  integer(HID_T) :: dataspace_id, attribute_id
  integer(HSIZE_T), dimension(1:1) :: dims
  PetscMPIInt :: hdf5_err

  dims = length
  call h5screate_simple_f(1, dims, dataspace_id, hdf5_err)
  call h5dcreate_f(parent_id, name, type, dataspace_id, attribute_id, hdf5_err)
  call h5dwrite_f(attribute_id, type, data, dims, hdf5_err)
  call h5dclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5DatasetStringArray

! ************************************************************************** !

subroutine OutputHDF5WriteSnapShotAtts(parent_id,option)
  ! 
  ! Writes attributes associated with a snapshot time in the output file.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/31/19
  ! 
  use Option_module

  implicit none

  integer(HID_T) :: parent_id
  type(option_type) :: option

  integer(HID_T) :: attribute_id
  integer(HID_T) :: dataspace_id
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  integer(HSIZE_T) :: dims(1)
  PetscMPIInt :: hdf5_err

  dims = 1
  call h5screate_simple_f(1,dims,dataspace_id,hdf5_err)
  string = 'Time (s)'
  call h5acreate_f(parent_id,string,H5T_NATIVE_DOUBLE,dataspace_id, &
                   attribute_id,hdf5_err)
  call h5awrite_f(attribute_id,H5T_NATIVE_DOUBLE,option%time,dims,hdf5_err)
  call h5aclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5WriteSnapShotAtts

end module Output_HDF5_Aux_module
