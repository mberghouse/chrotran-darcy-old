module Subsurface_Reset_module

  ! The purpose of this module is provide objects and function that allow
  ! the user to overwrite portions of the solution vector at points in time.
 
#include "petsc/finclude/petscsys.h"
  use petscsys 
  use PFLOTRAN_Constants_module
  use Dataset_Base_class

  implicit none

  private

  type, public :: subsurface_reset_type
    PetscReal :: time
    character(len=MAXWORDLENGTH), pointer :: material_names(:)
    character(len=MAXWORDLENGTH) :: dataset_name
    type(dataset_base_type), pointer :: dataset
    type(subsurface_reset_type), pointer :: next
  end type subsurface_reset_type
  
  public :: SubsurfaceResetCreate, &
            SubsurfaceResetRead, &
            SubsurfaceResetAddToList, &
            SubsurfaceResetDestroy

contains

! ************************************************************************** !

function SubsurfaceResetCreate()
  ! 
  ! Creates a subsurface reset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/18/17
  ! 
  
  implicit none

  type(subsurface_reset_type), pointer :: SubsurfaceResetCreate
  
  type(subsurface_reset_type), pointer :: subsurface_reset
  
  allocate(subsurface_reset)
  subsurface_reset%time = UNINITIALIZED_DOUBLE
  nullify(subsurface_reset%material_names)
  subsurface_reset%dataset_name = ''
  nullify(subsurface_reset%dataset)
  nullify(subsurface_reset%next)
  SubsurfaceResetCreate => subsurface_reset

end function SubsurfaceResetCreate

! ************************************************************************** !

subroutine SubsurfaceResetRead(subsurface_reset,input,option)
  ! 
  ! Reads in contents of a subsurface reset card
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/18/17
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Units_module

  implicit none
  
  type(subsurface_reset_type) :: subsurface_reset
  type(input_type), pointer :: input
  type(option_type) :: option
  
  PetscInt :: icount
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: materials(100)

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','FLUID_PROPERTY')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('TIME') 
        call InputReadDouble(input,option,subsurface_reset%time)
        call InputErrorMsg(input,option,'time', &
                           'SOLUTION_RESET')
        call InputReadAndConvertUnits(input,subsurface_reset%time, &
                                      'sec','SOLUTION_RESET,time',option)
      case('MATERIALS')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          icount = icount + 1
          if (icount > 100) then
            option%io_buffer = 'Material name buffer in SubsurfaceResetRead &
              &must be increased.'
            call printErrMsg(option)
          endif
          call InputReadWord(input,option,materials(icount),PETSC_TRUE)
          call InputErrorMsg(input,option,'material','SOLUTION_RESET')
        enddo
        allocate(subsurface_reset%material_names(icount))
        subsurface_reset%material_names = materials(1:icount)
      case('DATASET')
        call InputReadWord(input,option,subsurface_reset%dataset_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'dataset','SOLUTION_RESET')
      case default
        call InputKeywordUnrecognized(keyword,'SOLUTION_RESET',option)
    end select
    
  enddo  

end subroutine SubsurfaceResetRead

! ************************************************************************** !

subroutine SubsurfaceResetAddToList(subsurface_reset,list)
  ! 
  ! Adds a subsurface reset object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/18/17
  ! 

  implicit none
  
  type(subsurface_reset_type), pointer :: subsurface_reset
  type(subsurface_reset_type), pointer :: list

  type(subsurface_reset_type), pointer :: cur_subsurface_reset
  
  if (associated(list)) then
    cur_subsurface_reset => list
    ! loop to end of list
    do
      if (.not.associated(cur_subsurface_reset%next)) exit
      cur_subsurface_reset => cur_subsurface_reset%next
    enddo
    cur_subsurface_reset%next => subsurface_reset
  else
    list => subsurface_reset
  endif
  
end subroutine SubsurfaceResetAddToList

! ************************************************************************** !

recursive subroutine SubsurfaceResetDestroy(subsurface_reset)
  ! 
  ! Destroys a subsurface reset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/18/17
  ! 

  implicit none
  
  type(subsurface_reset_type), pointer :: subsurface_reset
  
  if (.not.associated(subsurface_reset)) return
  
  call SubsurfaceResetDestroy(subsurface_reset%next)

  deallocate(subsurface_reset)
  nullify(subsurface_reset)
  
end subroutine SubsurfaceResetDestroy

end module Subsurface_Reset_module
