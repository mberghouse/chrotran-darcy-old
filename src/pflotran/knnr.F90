!************************************************************
module kNNr_module

#include "petsc/finclude/petscsys.h"

  use kdtree2_module

  use petscvec
 

  use PFLOTRAN_Constants_module
  use Utility_module

  implicit none

  type(kdtree2), pointer :: tree
  PetscReal, dimension(:,:), allocatable :: my_array
  PetscInt :: n,d
  real, allocatable :: table_data(:,:)
  PetscReal :: eps = tiny(0.0d0)


contains
  
  subroutine read_my_csv_data ( )
!*****************************************************************************80
!
!! reads a CSV table data file, with one header line, and real valued entries after that.
!
!*****************************************************************************80

  implicit none

  

  character ( len = 80 ) :: csv_file_name = 'rosie_test_data.csv'

  integer   ( kind = 4 ) csv_file_status
  integer   ( kind = 4 ) csv_file_unit
  integer   ( kind = 4 ) csv_record_status
  integer   ( kind = 8 ) i,i_v,vl,al
  integer   ( kind = 8 ) line_num
  character ( len = 200 ) record
  integer   ( kind = 4 ) value_count, this_value_count
  character ( len = 30 ), dimension(10) :: my_values ! hard wired for now

  vl = 30 !value length in my_values array
  al = 10 !length of array my_values


  call csv_file_line_count ( csv_file_name, line_num )


  call csv_file_open_read ( csv_file_name, csv_file_unit )

  read ( IUNIT_TEMP, '(a)', iostat = csv_file_status ) record

  call csv_value_count ( record, csv_record_status, value_count )


  allocate(table_data(line_num-1,value_count))



  do i = 1, line_num - 1

    read ( IUNIT_TEMP, '(a)', iostat = csv_file_status ) record


    ! TODO: Need to add some error checking that throws a fit if any lines have less than the expected number of records
    call csv_values_extract ( record, csv_record_status, this_value_count, my_values, vl, al)

!    if (this_value_count .ne. value_count) then
!        write (*,*) 'Poorly formed input line'
!    endif

    do i_v = 1, value_count
      read(my_values(i_v), *) table_data(i,i_v)
    end do
  end do

  call csv_file_close_read ( csv_file_name )



  return
end subroutine read_my_csv_data

subroutine csv_value_count ( csv_record, csv_record_status, value_count )

!*****************************************************************************80
!
!! CSV_COUNT counts the number of values in a CSV record.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  character csv_char
  character csv_char_old
  integer ( kind = 4 ) csv_len
  integer ( kind = 4 ) csv_loc
  character ( len = * ) csv_record
  integer ( kind = 4 ) csv_record_status
  character :: TAB = achar ( 9 )
  integer ( kind = 4 ) value_count
  integer ( kind = 4 ) word_length
!
!  No values so far.
!
  value_count = 0
!
!  We begin in "unquoted" status.
!
  csv_record_status = 0
!
!  How many characters in the record?
!
  csv_len = len_trim ( csv_record )
!
!  Count number of characters in each word.
!
  word_length = 0
!
!  Consider each character.
!
  csv_char_old = ','

  do csv_loc = 1, csv_len

    csv_char = csv_record(csv_loc:csv_loc)
!
!  Each comma divides one value from another.
!
    if ( csv_char_old == ',' ) then

      value_count = value_count + 1
      word_length = 0
!
!  For quotes, try using CSV_RECORD_STATUS to count the number of
!  quoted characters.
!
    else if ( csv_char == '"' ) then

      if ( 0 < csv_record_status ) then
        csv_record_status = 0
      else
        csv_record_status = csv_record_status + 1
      end if
!
!  Ignore blanks
!
    else if ( csv_char == ' ' .or. csv_char == TAB ) then
!
!  Add character to length of word.
!
    else

      word_length = word_length + 1

      if ( value_count == 0 ) then
        value_count = 1
      end if

    end if

    csv_char_old = csv_char

  end do

  return
end subroutine csv_value_count

subroutine csv_values_extract (csv_record, csv_record_status, value_count,array_of_values,vl,al)

!*****************************************************************************80
!
!! extracts the values from a record
!
!  Licensing:
!
!
!
!  Modified:
!
!   02/11/2020
!
!  Author:
!
!    Bert Debusschere
!
!  Parameters:
!
  implicit none

  character csv_char
  character csv_char_old
  integer ( kind = 4 ) csv_len
  integer ( kind = 4 ) csv_loc
  character ( len = * ) csv_record
  integer ( kind = 4 ) csv_record_status
  character :: TAB = achar ( 9 )
  integer ( kind = 4 ) value_count
  integer ( kind = 8 ) word_length, vl
  integer (kind = 8) al
  character (len = vl ) :: array_of_values(al)

  character (len = 30) word  ! one word / value in record

  ! Start with assumption there is at least 1 value
  value_count = 1

!
!  We begin in "unquoted" status.
!
  csv_record_status = 0
!
!  How many characters in the record?
!
  csv_len = len_trim ( csv_record )
!
!  Count number of characters in each word.
!
  word_length = 0
  word = ''

!
!  Consider each character.
!
  ! csv_char_old = ','

  do csv_loc = 1, csv_len

    csv_char = csv_record(csv_loc:csv_loc)
    ! write(*, *) value_count, word, csv_char
!
!  Each comma divides one value from another.
!
    if ( csv_char == ',' ) then ! we met the end of a word / value

      ! write(*,*) 'Assigning the word ',trim(word), ' to the array'
      array_of_values(value_count) = trim(word)
      value_count = value_count + 1
      word_length = 0
      word = ''
!
!  For quotes, try using CSV_RECORD_STATUS to count the number of
!  quoted characters.
!
    else if ( csv_char == '"' ) then

      if ( 0 < csv_record_status ) then
        csv_record_status = 0
      else
        csv_record_status = csv_record_status + 1
      end if
!
!  Ignore blanks
!
    else if ( csv_char == ' ' .or. csv_char == TAB ) then
!
!  Add character to length of word.
!
    else

      word_length = word_length + 1
      word = trim(word) // csv_char

    end if

    csv_char_old = csv_char

  end do

  array_of_values(value_count) = trim(word) ! Add the last word to the array

  return
end subroutine csv_values_extract


subroutine csv_file_close_read ( csv_file_name )

!*****************************************************************************80
!
!! CSV_FILE_CLOSE_READ closes a CSV file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!

  implicit none

  character ( len = * ) csv_file_name

  close (IUNIT_TEMP )

  return
end subroutine csv_file_close_read



subroutine csv_file_line_count ( csv_file_name, line_num )

!*****************************************************************************80
!
!! CSV_FILE_LINE_COUNT counts the number of lines in a CSV file.
!
!  Discussion:
!
!    This routine does not try to distinguish the possible header line,
!    blank lines, or cases where a single CSV record extends over multiple
!    lines.  It simply counts the number of lines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!  use Option_module
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 1023 ) line
  integer ( kind = 8 ) line_num
!  type(option_type) :: option

  line_num = -1

  open (IUNIT_TEMP, file = csv_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
     print *,'CSV file ' // trim(csv_file_name) // 'not found.'
!     option%io_buffer = 'CSV file ' !// trim(csv_file_name) // ' not found.'
!     call PrintErrMsg(option,'csv file wrong')

    stop
  end if

  line_num = 0

  do

    read ( IUNIT_TEMP, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

  end do

  close ( IUNIT_TEMP )

  return
end subroutine csv_file_line_count



subroutine csv_file_open_read ( csv_file_name, csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_OPEN_READ opens a CSV file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!
  !

  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_status
  integer ( kind = 4 ) csv_file_unit


  open (IUNIT_TEMP, file = csv_file_name, status = 'old', &
    iostat = csv_file_status )

  
!  if ( csv_file_status /= 0 ) then
!     option%io_buffer = 'CSV file ' !// trim(csv_file_name) // 'not found.'
!     call PrintErrMsg(option)

!    stop
!  end if

  return
end subroutine csv_file_open_read

!***********************************
!***********************************

subroutine read_my_h5_file()

  use hdf5
  use HDF5_Aux_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: h5_name = 'test.h5'
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: group_name = '/'
  character(len=MAXSTRINGLENGTH) :: dataset_name1 = 'Temp'
  character(len=MAXSTRINGLENGTH) :: dataset_name2 = 'Env_CO3_2n'
  character(len=MAXSTRINGLENGTH) :: dataset_name3 = 'Env_O2'
  character(len=MAXSTRINGLENGTH) :: dataset_name4 = 'Env_Fe_2p'
  character(len=MAXSTRINGLENGTH) :: dataset_name5 = 'Env_H2'
  character(len=MAXSTRINGLENGTH) :: dataset_name6 = 'Dose Rate d0'
  character(len=MAXSTRINGLENGTH) :: dataset_name7 = 'UO2 Surface Flux'

  integer(HID_T) :: prop_id
  integer(HID_T) :: file_id
  integer(HID_T) :: parent_id
  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)

  integer :: ndims_h5
!  integer :: 
  real, dimension(:), allocatable :: dset_data
  
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

 
  call h5open_f(hdf5_err)
 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
 
  !hdf5openfilereadonlt
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

  !!!!!!!!!!!!!!!!!!! Repeat
  call h5dopen_f(group_id,dataset_name1,dataset_id,hdf5_err)
 
  if (hdf5_err < 0) then
    print *, "dataset not found"
    stop
  endif
 
  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
 
  ! 366 ???
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)

  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  !
 
  allocate(dset_data(dims_h5(1)))
  allocate(table_data(7,dims_h5(1)))
  ! 501
  
  !dset_data
  call h5dread_f(dataset_id,H5T_NATIVE_REAL, table_data(1,:), dims_h5, &
       hdf5_err)!, memory_space_id, file_space_id,prop_id)

  !  print *,log10(table_data(1,:))
 ! print *,dset_data(1)
  
  call h5dclose_f(dataset_id,hdf5_err)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call get_h5_dataset(group_id,dims_h5,dataset_name3,2)
  call get_h5_dataset(group_id,dims_h5,dataset_name3,3)
  call get_h5_dataset(group_id,dims_h5,dataset_name4,4)
  call get_h5_dataset(group_id,dims_h5,dataset_name5,5)
  call get_h5_dataset(group_id,dims_h5,dataset_name6,6)
  call get_h5_dataset(group_id,dims_h5,dataset_name7,7)
  
  call h5gclose_f(group_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)


  table_data(2,:) = log10(table_data(2,:))
  table_data(3,:) = log10(table_data(3,:))
  table_data(4,:) = log10(table_data(4,:))
  table_data(5,:) = log10(table_data(5,:))
  table_data(6,:) = log10(table_data(6,:))

  !testing purposes
!  deallocate(table_data)
end subroutine read_my_h5_file


subroutine get_h5_dataset(group_id,dims_h5,dataset_name,i)

  use hdf5
  use HDF5_Aux_module

  implicit none

  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

   integer :: ndims_h5


  integer(HSIZE_T),allocatable :: dims_h5(:), max_dims_h5(:)

  PetscInt :: i

  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: hdf5_err


    !!!!!!!!!!!!!!!!!!! Repeat
  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)
 
  if (hdf5_err < 0) then
    print *, "dataset not found"
    stop
  endif
 
  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

!    call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)

!  allocate(dims_h5(ndims_h5))
!  allocate(max_dims_h5(ndims_h5))
  
!  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
 
  
  !dset_data
  call h5dread_f(dataset_id,H5T_NATIVE_REAL, table_data(i,:), dims_h5, &
       hdf5_err)!, memory_space_id, file_space_id,prop_id)

  !  print *,log10(table_data(1,:))
  !print *, tabl
  
  call h5dclose_f(dataset_id,hdf5_err)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine get_h5_dataset

!***********************************
!***********************************

subroutine knnr_init()



  ! this is how you declare a tree in your main program

  PetscInt :: i_n,i_d


  PetscInt, dimension(2) :: data_array_shape

  call read_my_h5_file()
!  call read_my_csv_data ()
 

!  data_array_shape = shape(table_data)
!  n = data_array_shape(1)
!  d = data_array_shape(2)-1 ! Quantities of Interest (QoI) is not part of the search query of the search query.

!  allocate(my_array(d,n))
!  print *, 'n= ', n
!  print *, 'd= ', d
  
!  do i_d = 1, d
!    do i_n = 1, n
!       if (i_d == 1) then
!          my_array(i_d,i_n) = table_data(i_n,i_d)
!       else
!          my_array(i_d,i_n) = log10(table_data(i_n,i_d))
!       endif     
!    end do
! end do

  data_array_shape = shape(table_data)
  n = data_array_shape(1)-1
  d = data_array_shape(2) ! Quantities of Interest (QoI) is not part of the search query of the search query.

  allocate(my_array(n,d))
  my_array(1,:) = table_data(1,:)
  my_array(2,:) = table_data(2,:)
  my_array(3,:) = table_data(3,:)
  my_array(4,:) = table_data(4,:)
  my_array(5,:) = table_data(5,:)
  my_array(6,:) = table_data(6,:)
  
!  print *, my_array
 !my_array
 tree => kdtree2_create(my_array,sort=.false.,rearrange=.false.) 


  end subroutine knnr_init
  

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
      
  type(kdtree2_result),allocatable :: results(:)


  !Testing parameters
  integer   :: rind

  real(kdkind) :: rv

  
  yTme = sTme/60.0d0/60.0d0/24.0d0/365.0d0   !DAYS_PER_YEAR

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
  f(6) = log10(rad0)



  !Testing Purposes
!  call random_number(rv)
!  rind = floor(rv*tree%n)+1
!  f = tree%the_data(:,rind)

  
  nn=7

  allocate(results(nn))


  call kdtree2_n_nearest(tp=tree,qv=f,nn=nn,results=results)

  call inverse_distance(results,nn,qoi_ave)

  fuelDisRate = (qoi_ave)*270.0 !convert units

!  print *, 'mol/m2/yr', fuelDisRate/270.0
!  print *, 'known value', table_data(rind,d+1)

end subroutine knnr_query

subroutine knnr_close()
   deallocate(table_data)
   deallocate(my_array)
 end subroutine knnr_close

subroutine inverse_distance(results,nn,qoi_ave)

   PetscReal :: qoi_i, qoi_sum, qoi_ave, qoi_weights, weight, dis
   type(kdtree2_result)::myresult
   type(kdtree2_result),allocatable :: results(:), resultsb(:)

   PetscInt :: i_d,nn

   
   do i_d = 1,nn
      myresult = results(i_d)
      qoi_i = table_data(myresult%idx,d+1)
      dis = myresult%dis
     
      if (abs(dis) <= eps) then
         qoi_weights = 1.0
         qoi_sum = qoi_i

         exit
      elseif (isinfinite(abs(1/dis))) then
         qoi_weights = 1.0
         qoi_sum = qoi_i
         
         exit
      else 

         weight = 1 / dis

         qoi_sum = qoi_sum +qoi_i *weight
      

         qoi_weights = qoi_weights+weight
      endif

   end do
 
!   print *,'qoi_sum=',qoi_sum
!   print *, 'qoi_weights=',qoi_weights
   qoi_ave = qoi_sum/qoi_weights
   

 end subroutine inverse_distance

  function isinfinite(value1)

    implicit none

    PetscBool :: isinfinite
    PetscReal :: value1
    PetscReal :: infinity


   isinfinite = .false.
   
  
   infinity = huge(0.0d0)
 

   if (value1 >= infinity) then
      isinfinite =.true.
   endif

   

 end function isinfinite

end module kNNr_module







 

