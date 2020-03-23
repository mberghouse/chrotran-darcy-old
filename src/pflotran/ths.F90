module THS_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use THS_Aux_module
  use THS_Common_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module
  
  public :: THSSetup, &
            THSInitializeTimestep, &
            THSUpdateSolution, &
            THSUpdateAuxVars, &
            THSUpdateFixedAccum, &
            THSComputeMassBalance, &
            THSResidual, &
            THSJacobian, &
            THSGetTecplotHeader, &
            THSSetPlotVariables, &
            THSMapBCAuxVarsToGlobal, &
            THSDestroy
            
contains

! ************************************************************************** !

subroutine THSSetup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Michael Nole
  ! Date: 03/19/20
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Material_Aux_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(material_parameter_type), pointer :: material_parameter
  PetscBool :: error_found
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%THS => THSAuxCreate(option)
  
  ths_analytical_derivatives = .not.option%flow%numerical_derivatives
  
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE
  

end subroutine THSSetup

! ************************************************************************** !

subroutine THSInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Upwind_Direction_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  ths_newton_iteration_number = -999
  update_upwind_direction = PETSC_FALSE
  
  call THSUpdateFixedAccum(realization)

end subroutine THSInitializeTimestep

! ************************************************************************** !

subroutine THSUpdateSolution(realization)
  !
  ! Updates data in module after successful timestep
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization


end subroutine THSUpdateSolution

! ************************************************************************** !

subroutine THSUpdateAuxVars(realization)
  !
  ! Updates auxiliary variables associated with THS
  !
  ! Author: Michael Nole
  ! Date: 03/22/20
  !
  
  use Realization_Subsurface_class
  use Patch_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
end subroutine THSUpdateAuxVars

! ************************************************************************** !

subroutine THSUpdateFixedAccum(realization)
  !
  ! Updates fixed portion of the accumulation term
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

end subroutine THSUpdateFixedAccum

! ************************************************************************** !

subroutine THSComputeMassBalance(realization,mass_balance)
  !
  ! Initializes mass mass balance
  !
  ! Author: Michael Nole
  ! Date: 03/22/20
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)
                            
  

end subroutine THSComputeMassBalance

! ************************************************************************** !

subroutine THSResidual(snes,xx,r,realization,ierr)
  !
  ! Computes the residual
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Option_module
  
  implicit none
  
  SNES :: snes
  Vec :: xx,r
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  

end subroutine THSResidual

! ************************************************************************** !

subroutine THSJacobian(snes,xx,A,B,realization,ierr)
  !
  ! Computes the Jacobian
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none
  
  SNES :: snes
  Vec :: xx
  Mat :: A,B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  

end subroutine THSJacobian

! ************************************************************************** !

subroutine THSGetTecplotHeader(realization,icolumn)
  !
  ! Returns THS contribution to tecplot file header
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn

end subroutine THSGetTecplotHeader

! ************************************************************************** !

subroutine THSSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

end subroutine THSSetPlotVariables

! ************************************************************************** !

subroutine THSMapBCAuxVarsToGlobal(realization)
  !
  ! Maps variables in THS auxvar to global equivalent
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization


end subroutine THSMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine THSDestroy(realization)
  !
  ! Deallocates variables associated with THS
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

end subroutine THSDestroy


end module THS_module
