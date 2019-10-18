module Output_HDF5_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  use Output_HDF5_Aux_module
  
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  PetscMPIInt, private, parameter :: ON=1, OFF=0

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: hdf5_first
  
  public :: OutputHDF5Init, &
            OutputHDF5, &
            OutputHDF5UGridXDMF, &
            OutputHDF5UGridXDMFExplicit
            

contains

! ************************************************************************** !

subroutine OutputHDF5Init(num_steps)
  ! 
  ! Initializes module variables for HDF5 output
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

  use Option_module

  implicit none
  
  PetscInt :: num_steps
  
  if (num_steps == 0) then
    hdf5_first = PETSC_TRUE
  else
    hdf5_first = PETSC_FALSE
  endif

end subroutine OutputHDF5Init

! ************************************************************************** !

subroutine OutputHDF5(realization_base,var_list_type,h5file)
  ! 
  ! Print to HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module
  
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module
  use HDF5_Aux_module
  
  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type
  type(output_hdf5_type) :: h5file

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  
  Vec :: global_vec
  Vec :: global_vec_vx
  Vec :: global_vec_vy
  Vec :: global_vec_vz
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: h5_filename
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscBool :: include_gas_phase
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  grid => patch%grid
  output_option => realization_base%output_option

  call OutputHDF5GetH5Filename(option,output_option,h5file, &
                               var_list_type,h5_filename)
  call OutputHDF5OpenFile(option,h5file,h5_filename,file_id)
  if (h5file%first_write) then
    call OutputHDF5Provenance(option,output_option,file_id)
    call OutputHDF5WriteStructCoord(realization_base,file_id)
  endif
  
  group_name = OutputHDF5GetGroupName_Time(option,output_option)
  call OutputHDF5OpenGroup(file_id,grp_id,group_name)

  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)
  
  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vx)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vy)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vz)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        string = cur_variable%name
        call StringSwapChar(string," ","_")
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                             global_vec,grp_id, &
                                             H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                             global_vec,grp_id, &
                                             H5T_NATIVE_INTEGER)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      cur_variable => output_option%aveg_output_variable_list%first
      do ivar = 1,output_option%aveg_output_variable_list%nvars
        string = 'Aveg. ' // cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif

        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           field%avg_vars_vec(ivar),grp_id, &
                                           H5T_NATIVE_DOUBLE)

        cur_variable => cur_variable%next
      enddo

  end select

  include_gas_phase = PETSC_FALSE
  if (option%nphase > 1 .or. option%transport%nphase > 1) then
    include_gas_phase = PETSC_TRUE
  endif
  if (output_option%print_hdf5_vel_cent .and. &
      (var_list_type==INSTANTANEOUS_VARS)) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base, global_vec_vx, &
                                         global_vec_vy,global_vec_vz, &
                                         LIQUID_PHASE)

    string = "Liquid X-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vx,grp_id,H5T_NATIVE_DOUBLE)

    string = "Liquid Y-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vy,grp_id,H5T_NATIVE_DOUBLE)

    string = "Liquid Z-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vz,grp_id,H5T_NATIVE_DOUBLE)

    if (include_gas_phase) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             GAS_PHASE)
        string = "Gas X-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vx,grp_id, &
                                           H5T_NATIVE_DOUBLE)

        string = "Gas Y-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vy,grp_id, &
                                           H5T_NATIVE_DOUBLE)

        string = "Gas Z-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vz,grp_id, &
                                           H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_vel_face .and. &
     (var_list_type==INSTANTANEOUS_VARS)) then

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     X_DIRECTION,grp_id,h5file)
        if (include_gas_phase) then
          string = "Gas X-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       X_DIRECTION,grp_id,h5file)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     Y_DIRECTION,grp_id,h5file)
        if (include_gas_phase) then
          string = "Gas Y-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       Y_DIRECTION,grp_id,h5file)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     Z_DIRECTION,grp_id,h5file)
        if (include_gas_phase) then
          string = "Gas Z-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       Z_DIRECTION,grp_id,h5file)
        endif
    endif
   
  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(option, h5file, file_id)

end subroutine OutputHDF5

! ************************************************************************** !

subroutine OutputHDF5UGridXDMF(realization_base,var_list_type,h5file)
  ! 
  ! This routine writes unstructured grid data in HDF5 XDMF format.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/2012
  ! 

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use HDF5_Aux_module
  
  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type
  type(output_hdf5_type) :: h5file

  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: global_vec_vx,global_vec_vy,global_vec_vz
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr

  character(len=MAXSTRINGLENGTH) :: filename_header, h5_filename
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err  
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscInt :: vert_count
  Vec :: ivec
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      xmf_filename = OutputFilename(output_option,option,'xmf','')
    case (AVERAGED_VARS)
      xmf_filename = OutputFilename(output_option,option,'xmf','aveg')
  end select
  
  call OutputHDF5GetH5Filename(option,output_option,h5file, &
                               var_list_type,h5_filename)
  call OutputHDF5OpenFile(option,h5file,h5_filename,file_id)
  
  strings => StringSplit(h5_filename,'/')
  filename_header = strings(size(strings))
  deallocate(strings)

  grid => patch%grid

  if (h5file%first_write) then
    call WriteHDF5CoordinatesUGridXDMF(realization_base,option,file_id)
  endif

  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write xmf output file: ' // trim(h5_filename)
    call PrintMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeader(OUTPUT_UNIT, &
                         option%time/output_option%tconv, &
                         grid%nmax, &
                         realization_base%output_option%xmf_vert_len, &
                         grid%unstructured_grid%num_vertices_global,&
                         filename_header,PETSC_TRUE)
  endif

  group_name = OutputHDF5GetGroupName_Time(option,output_option)
  call OutputHDF5OpenGroup(file_id,grp_id,group_name)
  
  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vx)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vy)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vz)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        call DiscretizationGlobalToNatural(discretization,global_vec, &
                                           natural_vec,ONEDOF)
        string = cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename_header) // ":/" // trim(group_name) // &
                          "/" // trim(string)
        if (option%myrank == option%io_rank) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,CELL_CENTERED_OUTPUT_MESH)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      if (associated(output_option%aveg_output_variable_list%first)) then
        cur_variable => output_option%aveg_output_variable_list%first
        do ivar = 1,output_option%aveg_output_variable_list%nvars
          string = 'Aveg. ' // cur_variable%name
          if (len_trim(cur_variable%units) > 0) then
            word = cur_variable%units
            call HDF5MakeStringCompatible(word)
            string = trim(string) // ' [' // trim(word) // ']'
          endif

          call DiscretizationGlobalToNatural(discretization, &
                                             field%avg_vars_vec(ivar), &
                                             natural_vec,ONEDOF)
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_DOUBLE)
          att_datasetname = trim(filename_header) // ":/" // trim(group_name) // &
                            "/" // trim(string)
          if (option%myrank == option%io_rank) then
            call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                    att_datasetname,CELL_CENTERED_OUTPUT_MESH)
          endif
          cur_variable => cur_variable%next
        enddo
      endif

  end select

  !Output flowrates
  if (output_option%print_hdf5_mass_flowrate.or. &
     output_option%print_hdf5_energy_flowrate.or. &
     output_option%print_hdf5_aveg_mass_flowrate.or. &
     output_option%print_hdf5_aveg_energy_flowrate) then

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        call OutputGetFaceFlowrateUGrid(realization_base)
        if (output_option%print_hdf5_mass_flowrate.or.&
           output_option%print_hdf5_energy_flowrate) then
          call WriteHDF5FlowratesUGrid(realization_base,option,grp_id, &
                                       var_list_type)
        endif
      case (AVERAGED_VARS)
        if (output_option%print_hdf5_aveg_mass_flowrate.or.&
           output_option%print_hdf5_aveg_energy_flowrate) then
          call WriteHDF5FlowratesUGrid(realization_base,option,grp_id, &
                                       var_list_type)
        endif
    end select
  endif

  if (output_option%print_hdf5_vel_cent .and. &
      (var_list_type==INSTANTANEOUS_VARS)) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                         global_vec_vy,global_vec_vz, &
                                         LIQUID_PHASE)
    do i = 1, 3
      select case(i)
        case(1)
          word = 'X'
          ivec = global_vec_vx
        case(2)
          word = 'Y'
          ivec = global_vec_vy
        case(3)
          word = 'Z'
          ivec = global_vec_vz
      end select
      string = 'Liquid ' // trim(word) // '-Velocity [m_per_' // &
               trim(output_option%tunit) // ']'
      call DiscretizationGlobalToNatural(discretization,ivec, &
                                         natural_vec,ONEDOF)
      call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
      att_datasetname = trim(filename_header) // ":/" // &
                        trim(group_name) // "/" // trim(string)
      if (option%myrank == option%io_rank) then
      call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string,att_datasetname, &
                              CELL_CENTERED_OUTPUT_MESH)
      endif
    enddo

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             GAS_PHASE)
      do i = 1, 3
        select case(i)
          case(1)
            word = 'X'
            ivec = global_vec_vx
          case(2)
            word = 'Y'
            ivec = global_vec_vy
          case(3)
            word = 'Z'
            ivec = global_vec_vz
        end select
        string = 'Gas ' // trim(word) // '-Velocity [m_per_' // &
                 trim(output_option%tunit) // ']'
        call DiscretizationGlobalToNatural(discretization,ivec, &
                                         natural_vec,ONEDOF)
        call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                     H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (option%myrank == option%io_rank) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,CELL_CENTERED_OUTPUT_MESH)
        endif
      enddo
    endif
  endif

  ! Output velocity at cell-face
  if (output_option%print_hdf5_vel_face) then

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        call OutputGetFaceVelUGrid(realization_base)
        if (output_option%print_hdf5_vel_face) then
          call WriteHDF5FaceVelUGrid(realization_base,option,grp_id, &
                                     var_list_type)
        endif
      case (AVERAGED_VARS)
    end select
  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(option, h5file, file_id)

  if (option%myrank == option%io_rank) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

end subroutine OutputHDF5UGridXDMF

! ************************************************************************** !

subroutine OutputHDF5UGridXDMFExplicit(realization_base,var_list_type,h5file)
  ! 
  ! This subroutine prints the explicit
  ! unstructured grid information in xdmf format
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/17/2013
  ! 

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use HDF5_Aux_module
  
  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type
  type(output_hdf5_type) :: h5file

  integer(HID_T) :: file_id, new_file_id, file_id2
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id, new_grp_id, grp_id2
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id, new_prop_id
  PetscMPIInt :: rank
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  integer(HSIZE_T) :: dims(3), max_dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr

  character(len=MAXSTRINGLENGTH) :: filename_header, h5_filename
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: domain_filename_path, domain_filename_header
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err  
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscInt :: istart
  PetscInt :: vert_count
  PetscBool :: write_xdmf
  PetscBool :: include_cell_centers
  PetscInt :: num_vertices, num_cells
  PetscInt :: mesh_type
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      xmf_filename = OutputFilename(output_option,option,'xmf','')
    case (AVERAGED_VARS)
      xmf_filename = OutputFilename(output_option,option,'xmf','aveg')
  end select
  
  call OutputHDF5GetH5Filename(option,output_option,h5file, &
                               var_list_type,h5_filename)
  call OutputHDF5OpenFile(option,h5file,h5_filename,file_id)

  grid => patch%grid

  domain_filename_path = trim(option%global_prefix) // '-domain.h5'
  strings => StringSplit(h5_filename,'/')
  domain_filename_header = strings(size(strings))
  deallocate(strings)

  write_xdmf = PETSC_FALSE
  include_cell_centers = PETSC_FALSE
  mesh_type = grid%unstructured_grid%explicit_grid%output_mesh_type
  if (option%myrank == option%io_rank .and. &
      (output_option%print_explicit_primal_grid .or. &
       len_trim(grid%unstructured_grid%explicit_grid% &
                  domain_filename) > 0)) then
    if (.not.output_option%print_explicit_primal_grid) then
      ! have to open up domain file read the size of the vertex array
      domain_filename_path = & 
        grid%unstructured_grid%explicit_grid%domain_filename
      domain_filename_header = domain_filename_path
        ! initialize fortran hdf5 interface 
      option%io_buffer = 'Opening hdf5 file: ' // trim(domain_filename_path)
!      call PrintMsg(option)
      call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
      call HDF5OpenFileReadOnly(domain_filename_path,file_id2,prop_id,option)
      call h5pclose_f(prop_id,hdf5_err)
      string = 'Domain/Cells'      
      call h5dopen_f(file_id2,string,data_set_id,hdf5_err)
      if (hdf5_err /= 0) then
        option%io_buffer = 'HDF5 dataset "' // trim(string) // '" not found &
          &in file "' // trim(domain_filename_path) // '".'
        call PrintErrMsg(option)
      endif
      call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
      ! should be a rank=2 data space
      call h5sget_simple_extent_dims_f(file_space_id,dims, &
                                       max_dims,hdf5_err) 
      num_cells = int(dims(1))
      call h5sclose_f(file_space_id,hdf5_err)
      call h5dclose_f(data_set_id,hdf5_err)  
      string = 'Domain/Vertices'      
      call h5dopen_f(file_id2,string,data_set_id,hdf5_err)
      if (hdf5_err /= 0) then
        option%io_buffer = 'HDF5 dataset "' // trim(string) // '" not found &
          &in file "' // trim(domain_filename_path) // '".'
        call PrintErrMsg(option)
      endif
      call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
      ! should be a rank=2 data space
      call h5sget_simple_extent_dims_f(file_space_id,dims, &
                                       max_dims,hdf5_err) 
      num_vertices = int(dims(2))
      call h5sclose_f(file_space_id,hdf5_err)
      call h5dclose_f(data_set_id,hdf5_err)  
      call h5fclose_f(file_id2,hdf5_err)
      include_cell_centers = PETSC_TRUE
    else
      ! for primal grid output, num_cells is set in the call to  
      ! WriteHDF5CoordinatesUGridXDMFExplicit() below.  Therefore, this value
      ! for num_cells will be overwritten the first time called.
      num_cells = realization_base%output_option%xmf_vert_len
      num_vertices = grid%unstructured_grid%explicit_grid%num_vertices
    endif
    write_xdmf = PETSC_TRUE
  endif
  
  if (first .and. output_option%print_explicit_primal_grid) then
    call h5pcreate_f(H5P_FILE_ACCESS_F,new_prop_id,hdf5_err)
    call h5fcreate_f(domain_filename_path,H5F_ACC_TRUNC_F,new_file_id, &
                     hdf5_err,H5P_DEFAULT_F,new_prop_id)
    call h5pclose_f(new_prop_id,hdf5_err)
    ! create a group for the coordinates data set
    string = "Domain"
    call h5gcreate_f(new_file_id,string,new_grp_id,hdf5_err, &
                     OBJECT_NAMELEN_DEFAULT_F)
    call WriteHDF5CoordinatesUGridXDMFExplicit(realization_base,option, &
                                               new_grp_id)
    num_cells = realization_base%output_option%xmf_vert_len
    call h5gclose_f(new_grp_id,hdf5_err)
    call h5fclose_f(new_file_id,hdf5_err)    
  endif   
  
  if (write_xdmf) then
    option%io_buffer = '--> write xmf output file: ' // trim(xmf_filename)
    call PrintMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeader(OUTPUT_UNIT, &
                       option%time/output_option%tconv, &
                       grid%unstructured_grid%explicit_grid%num_elems, &
                       num_cells, &
                       num_vertices, &
                       domain_filename_header,include_cell_centers)
  endif

  group_name = OutputHDF5GetGroupName_Time(option,output_option)
  if (len_trim(output_option%plot_name) > 2) then
    group_name = trim(group_name) // ' ' // output_option%plot_name
  endif
  call OutputHDF5OpenGroup(file_id,grp_id,group_name)
  
  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        call DiscretizationGlobalToNatural(discretization,global_vec, &
                                           natural_vec,ONEDOF)
        string = cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (write_xdmf) then
          !call OutputXMFAttributeExplicit(OUTPUT_UNIT,grid%nmax,string, &
          !                                att_datasetname)
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,mesh_type)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      if (associated(output_option%aveg_output_variable_list%first)) then
        cur_variable => output_option%aveg_output_variable_list%first
        do ivar = 1,output_option%aveg_output_variable_list%nvars
          string = 'Aveg. ' // cur_variable%name
          if (len_trim(cur_variable%units) > 0) then
            word = cur_variable%units
            call HDF5MakeStringCompatible(word)
            string = trim(string) // ' [' // trim(word) // ']'
          endif

          call DiscretizationGlobalToNatural(discretization, &
                                             field%avg_vars_vec(ivar), &
                                             natural_vec,ONEDOF)
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_DOUBLE)
          att_datasetname = trim(filename_header) // ":/" // &
                            trim(group_name) // "/" // trim(string)
          if (write_xdmf) then
            call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                    att_datasetname,mesh_type)
          endif
          cur_variable => cur_variable%next
        enddo
      endif

  end select

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  
  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(option, h5file, file_id)

  if (write_xdmf) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

end subroutine OutputHDF5UGridXDMFExplicit

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

subroutine WriteHDF5CoordinatesUGridXDMFExplicit(realization_base,option, &
                                                 file_id)
  ! 
  ! Writes the coordinates of
  ! explicit grid to HDF5 file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/17/2013
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module
  
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
  PetscInt :: istart
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)

  PetscReal, pointer :: vec_ptr(:)
  Vec :: natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscErrorCode :: ierr

  PetscInt :: TRI_ID_XDMF = 4
  PetscInt :: QUAD_ID_XDMF = 5
  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => realization_base%patch%grid

  allocate(vec_x_ptr(grid%unstructured_grid%explicit_grid%num_vertices))
  allocate(vec_y_ptr(grid%unstructured_grid%explicit_grid%num_vertices))
  allocate(vec_z_ptr(grid%unstructured_grid%explicit_grid%num_vertices))

  do i = 1, grid%unstructured_grid%explicit_grid%num_vertices 
    vec_x_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%x
    vec_y_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%y
    vec_z_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%z
  enddo
 
  !local_size = grid%unstructured_grid%explicit_grid%num_cells_global
  local_size = grid%unstructured_grid%explicit_grid%num_vertices
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  !POc
  !dims(2) = grid%unstructured_grid%explicit_grid%num_cells_global
  dims(2) = grid%unstructured_grid%explicit_grid%num_vertices
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
  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)


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

  deallocate(vec_x_ptr)
  deallocate(vec_y_ptr)
  deallocate(vec_z_ptr)
   
  !
  !  Write elements
  !
  local_size = grid%unstructured_grid%explicit_grid%num_elems

  call VecCreate(PETSC_COMM_SELF,natural_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(natural_vec,local_size*EIGHT_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(natural_vec,EIGHT_INTEGER,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(natural_vec,ierr);CHKERRQ(ierr)

  call OutputGetCellVerticesExplicit(grid,natural_vec)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  vert_count=0

  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%unstructured_grid%explicit_grid%num_elems
  
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  realization_base%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
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
  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (3)
        int_array(vert_count) = TRI_ID_XDMF
      case (4) 
        if (grid%unstructured_grid%grid_type /= TWO_DIM_GRID) then   
        ! Tetrahedron
          int_array(vert_count) = TET_ID_XDMF
        else
          int_array(vert_count) = QUAD_ID_XDMF
        endif
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
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
    
end subroutine WriteHDF5CoordinatesUGridXDMFExplicit

! ************************************************************************** !

subroutine WriteHDF5FlowratesUGrid(realization_base,option,file_id, &
                                   var_list_type)
  ! 
  ! This routine writes (mass/energy) flowrate for unstructured grid.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/19/2013
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module
  
  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: var_list_type  

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

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  type(ugdm_type),pointer :: ugdm
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: idual
  PetscInt :: iconn
  PetscInt :: face_id
  PetscInt :: local_id_up,local_id_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscInt :: iface_up,iface_dn
  PetscInt :: dof
  PetscInt :: sum_connection
  PetscInt :: offset
  PetscInt :: cell_type
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: iface
  PetscInt :: ndof

  PetscReal, pointer :: flowrates(:,:,:)
  PetscReal, pointer :: vec_ptr1(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: double_array(:)
  PetscReal :: dtime
  PetscInt :: istart

  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate

  Vec :: global_flowrates_vec
  Vec :: natural_flowrates_vec

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: unit_string

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  field => realization_base%field

  select case(option%iflowmode)
    case (RICHARDS_MODE,RICHARDS_TS_MODE)
      ndof=1
    case (TH_MODE,TH_TS_MODE)
      ndof=1
      if (output_option%print_hdf5_mass_flowrate .and. &
          output_option%print_hdf5_energy_flowrate) ndof = 2
    case default
      option%io_buffer='FLOWRATE output not supported in this mode'
      call PrintErrMsg(option)
  end select

  call VecGetLocalSize(field%flowrate_inst,local_size,ierr);CHKERRQ(ierr)
  local_size = local_size/(option%nflowdof*MAX_FACE_PER_CELL + 1)

  allocate(double_array(local_size*(MAX_FACE_PER_CELL+1)))
  double_array = 0.d0

  offset = option%nflowdof*MAX_FACE_PER_CELL+1

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      call VecGetArrayF90(field%flowrate_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      mass_flowrate = output_option%print_hdf5_mass_flowrate
      energy_flowrate = output_option%print_hdf5_energy_flowrate
    case (AVERAGED_VARS)
      call VecGetArrayF90(field%flowrate_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%flowrate_aveg,vec_ptr2,ierr);CHKERRQ(ierr)
      mass_flowrate = output_option%print_hdf5_aveg_mass_flowrate
      energy_flowrate = output_option%print_hdf5_aveg_energy_flowrate
  end select


  do dof = 1,option%nflowdof

    if (dof==1 .and. (.not.mass_flowrate)) exit
    if (dof==2 .and. (.not.energy_flowrate)) exit

    select case(option%iflowmode)
      case(RICHARDS_MODE,RICHARDS_TS_MODE)
        string = "Mass_Flowrate [kg_per_s]" // CHAR(0)
      case(TH_MODE,TH_TS_MODE)
        if (dof==1) then
          string = "Mass_Flowrate [kg_per_s]" // CHAR(0)
        else
          string = "Energy_Flowrate [MJ_per_s]" // CHAR(0)
        endif
      case default
        option%io_buffer='FLOWRATE output not implemented in this mode.'
        call PrintErrMsg(option)
    end select

    if (var_list_type==AVERAGED_VARS) string = 'Aveg_' // trim(string) // &
                                               char(0)

    ! memory space which is a 1D vector
    rank_mpi = 1
    dims = 0
    dims(1) = local_size*(MAX_FACE_PER_CELL+1)
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
    ! file space which is a 2D block
    rank_mpi = 2
    dims = 0
    dims(2) = ugrid%nmax
    dims(1) = MAX_FACE_PER_CELL + 1
    call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

    call h5eset_auto_f(OFF,hdf5_err)
    call h5dopen_f(file_id,trim(string),data_set_id,hdf5_err)
    hdf5_flag = hdf5_err
    call h5eset_auto_f(ON,hdf5_err)
    if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
      call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
      call h5dcreate_f(file_id,trim(string),H5T_NATIVE_DOUBLE,file_space_id, &
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
    length(1) = MAX_FACE_PER_CELL + 1

    stride = 1
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                              hdf5_err,stride,stride)
    ! write the data
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
      
        do i=1,local_size
          ! Num. of faces for each cell (Note: Use vec_ptr1 not vec_ptr2)
          double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
            vec_ptr1((i-1)*offset+1)
          ! Flowrate values for each face
          do iface = 1,MAX_FACE_PER_CELL
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
            vec_ptr1((i-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface + 1)
          enddo
        enddo

      case (AVERAGED_VARS)

        do i=1,local_size
          ! Num. of faces for each cell (Note: Use vec_ptr1 not vec_ptr2)
          double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
            vec_ptr1((i-1)*offset+1)
          ! Divide the flowrate values by integration 'time'
          do iface = 1,MAX_FACE_PER_CELL
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
            vec_ptr2((i-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface + 1)/ &
            output_option%periodic_snap_output_time_incr
          enddo
        enddo
    end select

    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    call h5pclose_f(prop_id,hdf5_err)

    call h5dclose_f(data_set_id,hdf5_err)
    call h5sclose_f(file_space_id,hdf5_err)

  enddo

  ! Free up memory
  deallocate(double_array)

end subroutine WriteHDF5FlowratesUGrid

! ************************************************************************** !

subroutine WriteHDF5FaceVelUGrid(realization_base,option,file_id, &
                                 var_list_type)
  !
  ! This routine writes velocity at cell faces for unstructured grid.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 05/25/2014
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: var_list_type

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

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  type(ugdm_type),pointer :: ugdm
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: idual
  PetscInt :: iconn
  PetscInt :: face_id
  PetscInt :: local_id_up,local_id_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscInt :: iface_up,iface_dn
  PetscInt :: iphase
  PetscInt :: sum_connection
  PetscInt :: offset
  PetscInt :: cell_type
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: idir
  PetscInt :: istart
  PetscInt :: iface

  PetscReal, pointer :: flowrates(:,:,:)
  PetscReal, pointer :: vx_ptr(:)
  PetscReal, pointer :: vy_ptr(:)
  PetscReal, pointer :: vz_ptr(:)
  PetscReal, pointer :: vec_ptr1(:)
  PetscReal, pointer :: double_array(:)
  PetscReal :: dtime

  Vec :: global_flowrates_vec
  Vec :: natural_flowrates_vec

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: unit_string
  character(len=1) :: string_dir

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  field => realization_base%field

  if (option%nphase == 1 .and. option%transport%nphase > 1) then
    option%io_buffer = 'WriteHDF5FaceVelUGrid not supported for gas &
      &transport without flow in the gas phase.'
    call PrintErrMsg(option)
  endif
  call VecGetLocalSize(field%vx_face_inst,local_size,ierr);CHKERRQ(ierr)
  local_size = local_size/(option%nphase*MAX_FACE_PER_CELL + 1)

  allocate(double_array(local_size*(MAX_FACE_PER_CELL+1)))
  double_array = 0.d0

  offset = option%nphase*MAX_FACE_PER_CELL+1

  do idir = 1,3

    select case (idir)
      case (1)
        string_dir = 'X'
        call VecGetArrayF90(field%vx_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      case (2)
        string_dir = 'Y'
        call VecGetArrayF90(field%vy_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      case (3)
        string_dir = 'Z'
        call VecGetArrayF90(field%vz_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
    end select

    do iphase = 1,option%nphase

      select case (iphase)
        case (LIQUID_PHASE)
          string = "Liquid " // string_dir // "-Velocity at cell face [m_per_" &
            // trim(output_option%tunit) // "]"
        case (GAS_PHASE)
          string = "Gas " // string_dir // "-Velocity at cell face [m_per_" // &
            trim(output_option%tunit) // "]"
      end select

      ! memory space which is a 1D vector
      rank_mpi = 1
      dims = 0
      dims(1) = local_size*(MAX_FACE_PER_CELL+1)
      call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

      ! file space which is a 2D block
      rank_mpi = 2
      dims = 0
      dims(2) = ugrid%nmax
      dims(1) = MAX_FACE_PER_CELL + 1
      call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

      call h5eset_auto_f(OFF,hdf5_err)
      call h5dopen_f(file_id,trim(string),data_set_id,hdf5_err)
      hdf5_flag = hdf5_err
      call h5eset_auto_f(ON,hdf5_err)
      if (hdf5_flag < 0) then
        ! if the dataset does not exist, create it
        call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
        call h5dcreate_f(file_id,trim(string),H5T_NATIVE_DOUBLE,file_space_id, &
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
      length(1) = MAX_FACE_PER_CELL + 1

      stride = 1
      call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                                hdf5_err,stride,stride)
      ! write the data
      call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                              hdf5_err)
#endif

      select case (var_list_type)
        case (INSTANTANEOUS_VARS)

          do i=1,local_size
            ! Num. of faces for each cell
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
              vec_ptr1((i-1)*offset+1)
            ! Flowrate values for each face
            do iface = 1,MAX_FACE_PER_CELL
              double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
              vec_ptr1((i-1)*offset + (iphase-1)*MAX_FACE_PER_CELL + iface + 1)* &
              realization_base%output_option%tconv
            enddo
          enddo

        case (AVERAGED_VARS)

      end select

      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

      call h5pclose_f(prop_id,hdf5_err)

      call h5dclose_f(data_set_id,hdf5_err)
      call h5sclose_f(file_space_id,hdf5_err)

    enddo

    select case (idir)
      case (1)
        call VecRestoreArrayF90(field%vx_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
      case (2)
        call VecRestoreArrayF90(field%vy_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
      case (3)
        call VecRestoreArrayF90(field%vz_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
    end select

  enddo

  ! Free up memory
  deallocate(double_array)

end subroutine WriteHDF5FaceVelUGrid

end module Output_HDF5_module
