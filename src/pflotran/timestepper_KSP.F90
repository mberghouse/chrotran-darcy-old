module Timestepper_KSP_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Solver_module
  use Timestepper_Base_class
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  type, public, extends(timestepper_base_type) :: timestepper_ksp_type

    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations
    PetscInt :: cumulative_wasted_linear_iterations
  
  contains
    
    procedure, public :: ReadInput => TimestepperKSPRead
    procedure, public :: Init => TimestepperKSPInit
    procedure, public :: StepDT => TimestepperKSPStepDT
    procedure, public :: UpdateDT => TimestepperKSPUpdateDT
    procedure, public :: CheckpointBinary => TimestepperKSPCheckpointBinary
    procedure, public :: CheckpointHDF5 => TimestepperKSPCheckpointHDF5
    procedure, public :: RestartBinary => TimestepperKSPRestartBinary
    procedure, public :: RestartHDF5 => TimestepperKSPRestartHDF5
    procedure, public :: Reset => TimestepperKSPReset
    procedure, public :: PrintInfo => TimestepperKSPPrintInfo
    procedure, public :: InputRecord => TimestepperKSPInputRecord
    procedure, public :: FinalizeRun => TimestepperKSPFinalizeRun
    procedure, public :: Strip => TimestepperKSPStrip
    procedure, public :: Destroy => TimestepperKSPDestroy
    
  end type timestepper_ksp_type
  
  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: stepper_ksp_header_type
    PetscInt :: cumulative_linear_iterations
  end type stepper_ksp_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: stepper_ksp_header_type
      implicit none
#include "petsc/finclude/petscbag.h"
      PetscBag :: bag
      class(stepper_ksp_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public :: TimestepperKSPCreate, TimestepperKSPPrintInfo, &
            TimestepperKSPInit

contains

! ************************************************************************** !

function TimestepperKSPCreate()
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  implicit none
  
  class(timestepper_ksp_type), pointer :: TimestepperKSPCreate
  
  class(timestepper_ksp_type), pointer :: this
  
  allocate(this)
  call this%Init()

  this%solver => SolverCreate()
  
  TimestepperKSPCreate => this
  
end function TimestepperKSPCreate

! ************************************************************************** !

subroutine TimestepperKSPInit(this)
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  implicit none
  
  class(timestepper_ksp_type) :: this

  call TimestepperBaseInit(this)

  this%num_linear_iterations = 0
  this%cumulative_linear_iterations = 0
  this%cumulative_wasted_linear_iterations = 0

  nullify(this%solver)

end subroutine TimestepperKSPInit

! ************************************************************************** !

subroutine TimestepperKSPRead(this,input,option)
  ! 
  ! Reads parameters associated with time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none

  class(timestepper_ksp_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER_KSP')
    call StringToUpper(keyword)   

    select case(trim(keyword))
      case default
        call TimestepperBaseProcessKeyword(this,input,option,keyword)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine TimestepperKSPRead

! ************************************************************************** !

subroutine TimestepperKSPUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  use PM_Base_class
  
  implicit none

  class(timestepper_ksp_type) :: this
  class(pm_base_type) :: process_model

  PetscReal, pointer :: dummy_array(:)
  PetscBool :: update_time_step

  update_time_step = PETSC_TRUE

  if (this%time_step_cut_flag) then
    this%num_constant_time_steps = 1
  else if (this%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    this%num_constant_time_steps = &
      this%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (this%num_constant_time_steps > &
      this%constant_time_step_threshold) then
    this%num_constant_time_steps = 0
  else if (this%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif

  if (update_time_step) then

    call process_model%UpdateTimestep(this%dt, &
                                      this%dt_min, &
                                      this%dt_max, &
                                      0,0,dummy_array, &  ! for SNES
                                      this%time_step_max_growth_factor)

  endif

end subroutine TimestepperKSPUpdateDT

! ************************************************************************** !

subroutine TimestepperKSPStepDT(this,process_model,stop_flag)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 
  use PM_Base_class
  use Option_module
  
  implicit none

  class(timestepper_ksp_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag
  
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  PetscInt :: icut

  solver => this%solver
  option => process_model%option
  
  call process_model%InitializeTimestep()

  do

    call process_model%PreSolve()

!    call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

     call process_model%Solve(option%time,ierr)

!    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

    if (.not. process_model%AcceptSolution()) then
      icut = icut + 1
      this%time_step_cut_flag = PETSC_TRUE
      ! add logic for shutdown, moving it from Timestepper_SNES to _Base
      this%target_time = this%target_time - this%dt
      this%dt = this%time_step_reduction_factor * this%dt
      this%target_time = this%target_time + this%dt
      option%dt = this%dt
      call process_model%TimeCut()
    else
      exit
    endif
  enddo

  this%steps = this%steps + 1
#if 0
!  this%cumulative_linear_iterations = &
!    this%cumulative_linear_iterations + sum_linear_iterations
!  this%cumulative_wasted_linear_iterations = &
!    this%cumulative_wasted_linear_iterations + sum_wasted_linear_iterations
!  this%cumulative_time_step_cuts = &
!    this%cumulative_time_step_cuts + icut

!  this%num_linear_iterations = num_linear_iterations

  ! print screen output
  call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
    write(*, '(/," Step ",i6," Time= ",1pe12.5, &
         & " Dt= ",1pe12.5, " [",a,"]",/, &
         & " linear = ",i5," [",i10,"]"," cuts = ",i2, " [",i4,"]")') &
         this%steps, &
         this%target_time/tconv, &
         this%dt/tconv, &
         trim(tunit), &
         sum_linear_iterations, &
         this%cumulative_linear_iterations,icut, &
         this%cumulative_time_step_cuts
  endif
  
  if (option%print_file_flag) then
    write(option%fid_out, '(/," Step ",i6," Time= ",1pe12.5, &
         & " Dt= ",1pe12.5, " [",a,"]",/, &
         & " linear = ",i5," [",i10,"]"," cuts = ",i2, " [",i4,"]")') &
         this%steps, &
         this%target_time/tconv, &
         this%dt/tconv, &
         trim(tunit), &
         sum_linear_iterations, &
         this%cumulative_linear_iterations,icut, &
         this%cumulative_time_step_cuts
  endif
#endif

  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
end subroutine TimestepperKSPStepDT

! ************************************************************************** !

subroutine TimestepperKSPCheckpointBinary(this,viewer,option)
  ! 
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 
  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_ksp_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

end subroutine TimestepperKSPCheckpointBinary

! ************************************************************************** !

subroutine TimestepperKSPRestartBinary(this,viewer,option)
  ! 
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_ksp_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

end subroutine TimestepperKSPRestartBinary

! ************************************************************************** !

subroutine TimestepperKSPCheckpointHDF5(this, h5_chk_grp_id, option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 12/02/19
  !
  use Option_module
  use hdf5
  use Checkpoint_module, only : CheckPointWriteIntDatasetHDF5
  use Checkpoint_module, only : CheckPointWriteRealDatasetHDF5

  implicit none

  class(timestepper_ksp_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  type(option_type) :: option

end subroutine TimestepperKSPCheckpointHDF5

! ************************************************************************** !

subroutine TimestepperKSPRestartHDF5(this, h5_chk_grp_id, option)
  !
  ! Restarts parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 12/02/19
  !
  use Option_module
  use hdf5
  use HDF5_Aux_module
  use Checkpoint_module, only : CheckPointReadIntDatasetHDF5
  use Checkpoint_module, only : CheckPointReadRealDatasetHDF5

  implicit none

  class(timestepper_ksp_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  type(option_type) :: option

end subroutine TimestepperKSPRestartHDF5

! ************************************************************************** !

subroutine TimestepperKSPReset(this)
  ! 
  ! Zeros timestepper object members.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  implicit none

  class(timestepper_ksp_type) :: this

  this%cumulative_linear_iterations = 0

  call TimestepperBaseReset(this)

end subroutine TimestepperKSPReset

! ************************************************************************** !

subroutine TimestepperKSPPrintInfo(this,option)
  ! 
  ! Prints information about time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  use Option_module
  use String_module
  
  implicit none
  
  class(timestepper_ksp_type) :: this
  type(option_type) :: option

  ! if adding new code, see TimestepperBEPrintInfo for an example of using
  ! multiple strings

  call TimestepperBasePrintInfo(this,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)

end subroutine TimestepperKSPPrintInfo

! ************************************************************************** !

subroutine TimestepperKSPInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 
  
  implicit none
  
  class(timestepper_ksp_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'initial timestep size: '
  write(word,*) this%dt_init
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperKSPInputRecord

! ************************************************************************** !

recursive subroutine TimestepperKSPFinalizeRun(this,option)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  use Option_module
  
  implicit none
  
  class(timestepper_ksp_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (OptionPrintToScreen(option)) then
    write(*,'(/,a," TS KSP steps = ",i6," cuts = ",i6)') &
            trim(this%name), &
            this%steps, &
            this%cumulative_time_step_cuts
    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,'(a)') trim(this%name) // ' TS KSP solver time = ' // &
      trim(adjustl(string)) // ' seconds'
  endif
  
end subroutine TimestepperKSPFinalizeRun

! ************************************************************************** !

subroutine TimestepperKSPStrip(this)
  ! 
  ! Deallocates members of a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  implicit none
  
  class(timestepper_ksp_type) :: this

  call TimestepperBaseStrip(this)
  call SolverDestroy(this%solver)
  
end subroutine TimestepperKSPStrip

! ************************************************************************** !

subroutine TimestepperKSPDestroy(this)
  ! 
  ! Deallocates a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/19
  ! 

  implicit none
  
  class(timestepper_ksp_type) :: this
  
  call TimestepperKSPStrip(this)
    
end subroutine TimestepperKSPDestroy

end module Timestepper_KSP_class
