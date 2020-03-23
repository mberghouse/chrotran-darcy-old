module PM_THS_class

#include "petsc/finclude/petscsnes.h"

  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  
  type, public, extends(pm_subsurface_flow_type) :: pm_ths_type

    PetscReal, pointer :: residual_abs_inf_tol(:) !MAN I think this should be of dimension n eqns
    PetscReal, pointer :: residual_scaled_inf_tol(:)
    PetscReal, pointer :: abs_update_inf_tol(:,:) !MAN Eqns x phase states?
    PetscReal, pointer :: rel_update_inf_tol(:,:) !MAN Eqnx x phase states?
    
    
  contains
    procedure, public :: ReadSimulationBlock => PMTHSRead
    procedure, public :: InitializeRun => PMTHSInitializeRun
    procedure, public :: UpdateTimestep => PMTHSUpdateTimestep
    procedure, public :: Residual => PMTHSResidual
    procedure, public :: Jacobian => PMTHSJacobian
    procedure, public :: CheckUpdatePre => PMTHSCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHSCheckUpdatePost
    procedure, public :: TimeCut => PMTHSTimeCut
    procedure, public :: UpdateSolution => PMTHSUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHSUpdateAuxVars
    procedure, public :: MaxChange => PMTHSMaxChange
    procedure, public :: InputRecort => PMTHSInputRecord
    procedure, public :: CheckpointBinary => PMTHSCheckpointBinary
    procedure, public :: RestartBinary => PMTHSRestartBinary
    procedure, public :: Destroy => PMTHSDestroy
  
  end type pm_ths_type

  public :: PMTHSCreate, &
            PMTHSSetFlowMode
  
contains

! ************************************************************************** !

function PMTHSCreate()
  !
  ! Creates THS process model shell
  !
  ! Author: Michael Nole
  ! Date: 03/05/20
  !
  
  implicit none
  
  class(pm_ths_type), pointer :: PMTHSCreate,this
  
  allocate(this)
  
  call PMSubsurfaceFlowInit(this)
  this%name = 'THS Flow'
  this%header = 'THS Flow'
  
  ! MAN: need to allocate tolerances after we know how many equations are being used.
!   this%residual_abs_inf_tol = residua_labs_inf_tol
!   this%residual_scaled_inf_tol = residual_scaled_inf_tol
!   this%abs_update_inf_tol = abs_update_inf_tol
!   this%rel_update_inf_tol = rel_update_inf_tol

  
  PMTHSCreate => this


end function PMTHSCreate

! ************************************************************************** !

subroutine PMTHSSetFlowMode(option)
  !
  ! Initializes default flow mode parameters for THS mode
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  
  option%iflowmode = THS_MODE
  option%nphase = 1
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3

  ! MAN: this should vary depending on inputs
  option%nflowdof = 3
  option%nflowspec = 2
  option%use_isothermal = PETSC_FALSE
  
end subroutine PMTHSSetFlowMode
  
! ************ Local subroutines ************ !
! ************************************************************************** !

subroutine PMTHSRead(this,input)
  !
  ! Reads THS options
  ! Author: Michael Nole
  ! Date: 03/05/20
  !
  
  Use THS_module
  Use THS_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  class(pm_ths_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  
  option => this%option
  error_string = 'THS Options'
  
  input%ierr = 0
  call InputPushBlock(input,option)
  
  do
  
    call InputReadPflotranString(input,option)
    
    if(InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found, &
                                        error_string,option)    
    if (found) cycle
    
    select case(trim(keyword))
    
    end select
    
  enddo
  
  call InputPopBlock(input,option)
  
  
end subroutine PMTHSRead
  
! ************************************************************************** !

recursive subroutine PMTHSInitializeRun(this)
  !
  ! Initializes timestepping
  !
  ! Author: Michael Nole
  ! Date: 04/05/20
  !
  
  use Realization_Base_class
  
  implicit none
  
  class(pm_ths_type) :: this
  
!   do i = 1, nvars !MAN need this to change depending on number of potential primary variables
!     call RealizationGetVariable(this%realization, &
!                                 this%realization%field%max_change_vecs(i), &
!                                 this%max_change_ivar(i), &
!                                 this%max_change_isubvar(i))
!   enddo
  
  call PMSubsurfaceFlowInitializeRun(this)
  

end subroutine PMTHSInitializeRun


! ************************************************************************** !

subroutine PMTHSUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations, tfac,&
                               time_step_max_growth_factor)
  !
  ! Updates timestep considering timestepping restrictions. 
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  implicit none
  
  class(pm_ths_type) :: this
  PetscReal :: dt, dt_min, dt_max
  PetscInt :: iacceleration, num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  if (initialized(this%cfl_governor)) then
  
  
  endif
                               
end subroutine PMTHSUpdateTimestep

! ************************************************************************** !

subroutine PMTHSResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
!   use THS_Module, only: THSResidual

  implicit none
  
  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: xx,r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
!   call THSResidual(snes,xx,r,this%realization,ierr)
  
end subroutine PMTHSResidual

! ************************************************************************** !

subroutine PMTHSJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
!   use THS_module, only: THSJacobian

  implicit none
  
  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A,B
  PetscErrorCode :: ierr
  
!   call THSJacobian(snes,xx,A,B,this%realization,ierr)
  
end subroutine PMTHSJacobian

! ************************************************************************** 

subroutine PMTHSCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Performs a check on the update before solving
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  implicit none
  
  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: X,dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  

end subroutine PMTHSCheckUpdatePre

! ************************************************************************** 

subroutine PMTHSCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                X1_changed,ierr)
  !
  ! Performs a check on the update after solving
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  implicit none
  
  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: X0,dX,X1
  PetscBool :: dX_changed,X1_changed
  PetscErrorCode :: ierr
  
  

end subroutine PMTHSCheckUpdatePost

! ************************************************************************** 

subroutine PMTHSCheckConvergence(this,snes,it,xnorm,unorm,fnorm,&
                                 reason,ierr)
  !
  ! Checks convergence using THS-specific criteria
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  implicit none
  
  class(pm_ths_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm,unorm,fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr


end subroutine PMTHSCheckConvergence

! ************************************************************************** 

subroutine PMTHSTimeCut(this)
  !
  ! Cuts time time
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
!   use THS_module, only : THSTimeCut
  
  implicit none

  class(pm_ths_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
!   call THSTimeCut(this%realization)

end subroutine PMTHSTimeCut

! ************************************************************************** 

subroutine PMTHSUpdateSolution(this)
  ! 
  ! Updates solution
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
 
!  use THS_module, only : THSUpdateSolution, &
!                         THSMapBCAuxVarsToGlobal
 
  implicit none
 
  class(pm_ths_type) :: this
 
  call PMSubsurfaceFlowUpdateSolution(this)
!  call THSUpdateSolution(this%realization)
!  call THSMapBCAuxVarsToGlobal(this%realization)
 

end subroutine PMTHSUpdateSolution

! ************************************************************************** 

subroutine PMTHSUpdateAuxVars(this)
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !

!   use THS_module, only: THSUpdateAuxVars

  implicit none
  
  class(pm_ths_type) :: this
  
!   call THSUpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE)
  
end subroutine PMTHSUpdateAuxVars

! ************************************************************************** 

subroutine PMTHSMaxChange(this)
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !

  implicit none
 
  class(pm_ths_type) :: this
 

end subroutine PMTHSMaxChange

! ************************************************************************** 

subroutine PMTHSInputRecord(this)
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !

  implicit none
 
  class(pm_ths_type) :: this
 
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id
 
  id = INPUT_RECORD_UNIT
 
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'ths'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_RESIDUAL: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMTHSInputRecord

! **************************************************************************

subroutine PMTHSCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with the THS PM
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_ths_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMTHSCheckpointBinary

! **************************************************************************

subroutine PMTHSRestartBinary(this,viewer)
  !
  ! Restarts data associatd with the THS PM
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_ths_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMTHSRestartBinary

! **************************************************************************

subroutine PMTHSDestroy(this)
  !
  ! Destroys a THS process model
  !
  ! Author: Michael Nole
  ! Date: 03/06/20
  !
  
!   use THS_module, only :: THSDestroy

  implicit none
  
  class(pm_ths_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif
  
!   deallocate(this%max_change_ivar)
!   nullify(this%max_change_ivar)
!   deallocate(this%max_change_isubvar)
!   nullify(this%max_change_isubvar)
  
!   call THSDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHSDestroy

! **************************************************************************

end module PM_THS_class
