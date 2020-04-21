module table_data_mod
  real ( kind = 8 ), allocatable :: table_data(:,:)
end module table_data_mod

!*****************************************************************************80


!************************************************************
module kNNr_module

#include "petsc/finclude/petscsys.h"
  use table_data_mod
  use kdtree2_module
!  use time_kdtree
  use petscvec
 
  use petscsnes

  implicit none

  type(kdtree2), pointer :: tree, tree2, tree3
  real(kdkind), dimension(:,:), allocatable :: my_array
  integer :: n, d

contains
  
  subroutine read_my_csv_data ( )
!*****************************************************************************80
!
!! reads a CSV table data file, with one header line, and real valued entries after that.
!
!*****************************************************************************80
  use table_data_mod
  implicit none

  character ( len = 80 ) :: csv_file_name = 'test_data.csv'

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

  read ( csv_file_unit, '(a)', iostat = csv_file_status ) record

  call csv_value_count ( record, csv_record_status, value_count )


  allocate(table_data(line_num-1,value_count))



  do i = 1, line_num - 1

    read ( csv_file_unit, '(a)', iostat = csv_file_status ) record


    ! TODO: Need to add some error checking that throws a fit if any lines have less than the expected number of records
    call csv_values_extract ( record, csv_record_status, this_value_count, my_values, vl, al)

!    if (this_value_count .ne. value_count) then
!        write (*,*) 'Poorly formed input line'
!    endif

    do i_v = 1, value_count
      read(my_values(i_v), *) table_data(i,i_v)
    end do
  end do

  call csv_file_close_read ( csv_file_name, csv_file_unit )



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

  ! write(*,*) 'Testing array of strings:'
  !
  !
  ! word = '235'
  ! array_of_values(1) = 'abc'
  ! write(*,*) array_of_values(1)

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

  ! write(*,*) 'Assigning the word ',trim(word), ' to the array'
  array_of_values(value_count) = trim(word) ! Add the last word to the array

  return
end subroutine csv_values_extract


subroutine csv_file_close_read ( csv_file_name, csv_file_unit )

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
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_unit

  close ( unit = csv_file_unit )

  return
end subroutine csv_file_close_read

subroutine csv_file_close_write ( csv_file_name, csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_CLOSE_WRITE closes a CSV file for writing.
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
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_unit

  close ( unit = csv_file_unit )

  return
end subroutine csv_file_close_write

subroutine csv_file_header_write ( csv_file_name, csv_file_unit, header )

!*****************************************************************************80
!
!! CSV_FILE_HEADER_WRITE writes a header to a CSV file.
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
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
!    Input, character ( len = * ) HEADER, the header.
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_unit
  character ( len = * ) header


  return
end subroutine csv_file_header_write

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
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 1023 ) line
  integer ( kind = 8 ) line_num

  line_num = -1

  call get_unit ( input_unit )

  open ( unit = input_unit, file = csv_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CSV_FILE_LINE_COUNT - Fatal error!'
!    write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
    stop
  end if

  line_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

  end do

  close ( unit = input_unit )

  return
end subroutine csv_file_line_count

subroutine csv_file_record_write ( csv_file_name, csv_file_unit, record )

!*****************************************************************************80
!
!! CSV_FILE_RECORD_WRITE writes a record to a CSV file.
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
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
!    Input, character ( len = * ) RECORD, the record.
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_unit
  character ( len = * ) record


  return
end subroutine csv_file_record_write

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
!    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_status
  integer ( kind = 4 ) csv_file_unit

  call get_unit ( csv_file_unit )

  open ( unit = csv_file_unit, file = csv_file_name, status = 'old', &
    iostat = csv_file_status )

  if ( csv_file_status /= 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CSV_FILE_OPEN_READ - Fatal error!'
!    write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
    csv_file_unit = - 1
    stop
  end if

  return
end subroutine csv_file_open_read

subroutine csv_file_open_write ( csv_file_name, csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_OPEN_WRITE opens a CSV file for writing.
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
!    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_status
  integer ( kind = 4 ) csv_file_unit

  call get_unit ( csv_file_unit )

  open ( unit = csv_file_unit, file = csv_file_name, status = 'replace', &
    iostat = csv_file_status )



  return
end subroutine csv_file_open_write





subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end subroutine get_unit






!***********************************
!***********************************
!***********************************
!***********************************

subroutine knnr_init()


 
 
  real(kdkind), allocatable :: query_vec(:)


  ! this is how you declare a tree in your main program

  integer :: i_n,i_d,k




  integer, dimension(2) :: data_array_shape

  print *, 'reading csv data'

  call read_my_csv_data ()

  print *, 'done reading csv data'

  

  data_array_shape = shape(table_data)
  n = data_array_shape(1)
  d = data_array_shape(2)-1 ! Quantities of Interest (QoI) is not part of the search query of the search query.

  allocate(my_array(d,n))

  print *, 'populating array'

  do i_d = 1, d
    do i_n = 1, n
       if (i_d == 1) then
          my_array(i_d,i_n) = table_data(i_n,i_d)
       else
          my_array(i_d,i_n) = log(table_data(i_n,i_d))
       endif     
    end do
 end do

! do i_n =1,n
!    my_array(1,i_n) = table_data (i_n,1)
 !end do


 print *, 'beginning tree'

 tree => kdtree2_create(my_array,sort=.false.,rearrange=.false.)  ! this is how you create a tree.

  print *, 'end tree'

  !!!!SAVE TREE!!!!

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
    PetscReal, dimension(6) :: f !!!change!!
    PetscReal :: yTme
    PetscReal :: f1, f2, f3, f4, f5
    PetscReal :: AOF, rad0a, rad0

   integer   :: nnbrute, rind, nn, i_d
  real      :: t0, t1, sps, avgnum, maxdeviation
  real(kdkind) :: rv, qoi_i, qoi_sum, qoi_ave, qoi_int

  real(kdkind), allocatable :: query_vec(:)

  integer, parameter  :: nnn = 3
  integer   :: nnarray(nnn)
  data nnarray / 1, 5, 10 /
      
  type(kdtree2_result),allocatable :: results(:), resultsb(:)
  type(kdtree2_result)::myresult

  !create query vector! ???
  print *, 'inside success'

  
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

  print *, f
  
  nn=7

  allocate(results(nn))
  print *, 'calling nearest'
  allocate(query_vec(6))
  query_vec = f

  call kdtree2_n_nearest(tp=tree,qv=f,nn=nn,results=results)

  print *, 'out'

    qoi_sum = 0.d0 

  do i_d = 1,nn
     myresult = results(i_d)

     qoi_i = table_data(myresult%idx,d+1)
     print *, qoi_i
    qoi_sum = qoi_sum + qoi_i
  end do

  print *, 'almost there'
  fuelDisRate = (qoi_sum/float(nn))*270.0

end subroutine knnr_query

subroutine knnr_close()
   deallocate(table_data)
   deallocate(my_array)
end subroutine knnr_close

end module kNNr_module







 

