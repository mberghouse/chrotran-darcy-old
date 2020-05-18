!************************************************************
module kNNr_module

#include "petsc/finclude/petscsys.h"

  use kdtree2_module
  use PFLOTRAN_Constants_module

  use petscvec

  implicit none

 ! private

  type(kdtree2), pointer :: tree
  PetscReal, allocatable :: my_array(:,:)
  PetscInt :: n,d
  real, allocatable :: table_data(:,:)
  PetscReal :: eps = tiny(0.0d0)


contains
  
! ************************************************************************** !

subroutine read_my_h5_file()

  use hdf5
 ! use HDF5_Aux_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: h5_name = 'test.h5'
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: group_name = '/'
  character(len=MAXSTRINGLENGTH) :: dataset_name

  integer(HID_T) :: prop_id
  integer(HID_T) :: file_id
  integer(HID_T) :: parent_id
  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)

  integer :: ndims_h5

  real, dimension(:), allocatable :: dset_data
  
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

 
  call h5open_f(hdf5_err)
 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
 
  call h5fopen_f(h5_name,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
 
  if (hdf5_err /= 0) then
     print *, 'h5 file not found'
     stop
  endif

  call h5pclose_f(prop_id,hdf5_err)

  !hdf5groupopen
  string = adjustl(group_name)
  call h5gopen_f(file_id,trim(string),group_id,hdf5_err)
  if (hdf5_err < 0) then
    print *, 'HDF5 Group "' // trim(string) // '" not found.'
    stop
  endif

  dataset_name = 'Temp'
  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)
 
  if (hdf5_err < 0) then
    print *, "dataset not found"
    stop
  endif
 
  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
 

  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)

  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)

  allocate(table_data(7,dims_h5(1)))

  call h5dread_f(dataset_id,H5T_NATIVE_REAL, table_data(1,:), dims_h5, &
       hdf5_err)
  
  call h5dclose_f(dataset_id,hdf5_err)

  dataset_name = 'Env_CO3_2n'
  call get_h5_dataset(group_id,dims_h5,dataset_name,2)
  dataset_name = 'Env_O2'
  call get_h5_dataset(group_id,dims_h5,dataset_name,3)
  dataset_name = 'Env_Fe_2p'
  call get_h5_dataset(group_id,dims_h5,dataset_name,4)
  dataset_name = 'Env_H2'
  call get_h5_dataset(group_id,dims_h5,dataset_name,5)
  dataset_name = 'Dose Rate d0'
  call get_h5_dataset(group_id,dims_h5,dataset_name,6)
  dataset_name = 'UO2 Surface Flux'
  call get_h5_dataset(group_id,dims_h5,dataset_name,7)

  deallocate(dims_h5)
  deallocate(max_dims_h5)
  
  call h5gclose_f(group_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)


  table_data(2,:) = log10(table_data(2,:))
  table_data(3,:) = log10(table_data(3,:))
  table_data(4,:) = log10(table_data(4,:))
  table_data(5,:) = log10(table_data(5,:))
  table_data(6,:) = log10(table_data(6,:))


end subroutine read_my_h5_file

! ************************************************************************** !

subroutine get_h5_dataset(group_id,dims_h5,dataset_name,i)

  use hdf5

  implicit none

  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T),allocatable :: dims_h5(:)

  PetscInt :: i

  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: hdf5_err

  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)
 
  if (hdf5_err < 0) then
    print *, "dataset not found"
    stop
  endif
 
  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  call h5dread_f(dataset_id,H5T_NATIVE_REAL, table_data(i,:), dims_h5, &
       hdf5_err)
  
  call h5dclose_f(dataset_id,hdf5_err)

end subroutine get_h5_dataset

! ************************************************************************** !

subroutine knnr_init()

  implicit none

  PetscInt :: i_n,i_d
  PetscInt :: data_array_shape(2)

  call read_my_h5_file()

  data_array_shape = shape(table_data)
  n = data_array_shape(1)-1
  d = data_array_shape(2) ! Quantities of Interest (QoI) is not part of the search query of the search query.

  allocate(my_array(n,d))

  do i_d = 1, n
    my_array(i_d,:) = table_data(i_d,:)
  end do
 
  tree => kdtree2_create(my_array,sort=.false.,rearrange=.false.) 

end subroutine knnr_init

! ************************************************************************** !

subroutine knnr_query (burnup,sTme,current_temp_C,decay_time,conc,fuelDisRate)

  implicit none

  PetscReal, intent(in) :: current_temp_C
  PetscReal, intent(in) :: decay_time
  PetscReal, intent(in) :: conc(:)
  PetscReal, intent(in) :: burnup
  PetscReal, intent(in) :: sTme

  PetscReal, intent(out) :: fuelDisRate

  ! features
  PetscReal, dimension(6) :: f
  PetscReal :: yTme
  PetscReal :: f1, f2, f3, f4, f5
  PetscReal :: AOF, rad0a, rad0

  PetscInt :: nn

  PetscReal :: qoi_ave
  PetscReal, parameter :: UO2_molar_mass = 270.0d0 !g/mol
      
  type(kdtree2_result),allocatable :: results(:)

  !Testing parameters
  integer   :: rind

  real(kdkind) :: rv
 
  yTme = sTme/60.0d0/60.0d0/24.0d0/DAYS_PER_YEAR  

  ! calculate dose rate at the fuel surface (rad0)
  AOF = yTme + decay_time

  f2 = log(AOF)
  f1 = f2**2.0d0
  f3 = 1.0d0/f2
  f4 = f2/AOF
  f5 = exp(burnup/25.26892627636246d0)

  rad0a = -206.0634818750711d0   - 0.7631591788870090d0*f1 &
        + 20.97112373957833d0*f2 + 678.8463343193430d0*f3 &
        - 506.7149017370657d0*f4 + 0.1555448893425319d0*f5
  rad0 = max(exp(rad0a),5.0d-3)

  f(1) = current_temp_C + 273.15d0
  f(2) = log10(conc(1)) ! Env_CO3_2n
  f(3) = log10(conc(2)) ! Env_O2
  f(4) = log10(conc(3)) ! Env_Fe_2p
  f(5) = log10(conc(4)) ! Env_H2
  f(6) = log10(rad0)    ! Dose Rate

  !Testing Purposes
!  call random_number(rv)
!  rind = floor(rv*tree%n)+1
!  f = tree%the_data(:,rind)

  !Nearest Neighbor 
  nn = 7

  allocate(results(nn))

  call kdtree2_n_nearest(tp=tree,qv=f,nn=nn,results=results)

  call inverse_distance(results,nn,qoi_ave)

  fuelDisRate = (qoi_ave) * UO2_molar_mass !convert units

!  print *, 'mol/m2/yr', fuelDisRate/270.0
!  print *, 'known value', table_data(rind,d+1)

end subroutine knnr_query

! ************************************************************************** !

subroutine knnr_close()

  implicit none
  
  deallocate(table_data)
  deallocate(my_array)
   
end subroutine knnr_close

! ************************************************************************** !

subroutine inverse_distance(results,nn,qoi_ave)

  implicit none

  PetscReal :: qoi_i, qoi_sum, qoi_ave, qoi_weights, weight, dis
  type(kdtree2_result) :: myresult
  type(kdtree2_result),allocatable :: results(:)

  PetscInt :: i_d,nn
  
  do i_d = 1,nn
    myresult = results(i_d)
    qoi_i = table_data(myresult%idx,d+1)
    dis = myresult%dis
     
    if ( abs(dis) <= eps ) then
      qoi_weights = 1.0
      qoi_sum = qoi_i

      exit
    elseif ( isinfinite(abs(1/dis)) ) then
      qoi_weights = 1.0
      qoi_sum = qoi_i
         
      exit
    else 

      weight = 1 / dis
      qoi_sum = qoi_sum + qoi_i * weight
      qoi_weights = qoi_weights+weight

    endif

  end do
 
!   print *,'qoi_sum=',qoi_sum
!   print *, 'qoi_weights=',qoi_weights
  qoi_ave = qoi_sum/qoi_weights
   
end subroutine inverse_distance

! ************************************************************************** !

function isinfinite(value1)

  implicit none

  PetscBool :: isinfinite
  PetscReal :: value1
  PetscReal :: infinity

  isinfinite = PETSC_FALSE
   
  
  infinity = huge(0.0d0)
 
  if ( value1 >= infinity ) then
    isinfinite = PETSC_TRUE
  endif
 
end function isinfinite

end module kNNr_module







 

