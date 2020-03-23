module THS_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  
  implicit none
  
  private
  
  PetscBool, public :: ths_analytical_derivatives = PETSC_FALSE
  PetscInt, public :: ths_newton_iteration_number = 0

  type, public :: ths_auxvar_type
    PetscReal, pointer :: pres(:)
    PetscReal, pointer :: sat(:)
    PetscReal, pointer :: den(:)
    PetscReal, pointer :: den_kg(:)
    PetscReal :: temp
    PetscReal, pointer :: xmol(:,:)
    PetscReal, pointer :: H(:)
    PetscReal, pointer :: U(:)
    PetscReal, pointer :: kr(:)
    PetscReal, pointer :: mobility(:)
    PetscReal :: effective_porosity
    PetscReal :: pert
    type(ths_derivative_auxvar_type), pointer :: d
  end type ths_auxvar_type
  
  type, public :: ths_derivative_auxvar_type
  
  end type ths_derivative_auxvar_type
  
  type, public :: ths_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscBool :: check_post_converged
  end type ths_parameter_type
  
  type, public :: ths_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux
    PetscInt :: num_aux_bc
    PetscInt :: num_aux_ss
    type(ths_parameter_type), pointer :: ths_parameter
    type(ths_auxvar_type), pointer :: auxvars(:,:)
    type(ths_auxvar_type), pointer :: auxvars_bc(:)
    type(ths_auxvar_type), pointer :: auxvars_ss(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type ths_type
  
  interface THSAuxVarDestroy
    module procedure THSAuxVarSingleDestroy
    module procedure THSAuxVarArray1Destroy
    module procedure THSAuxVarArray2Destroy
  end interface THSAuxVarDestroy
  
  interface THSOutputAuxVars
    module procedure THSOutputAuxVars1
    module procedure THSOutputAuxVars2
  end interface THSOutputAuxVars
  
  public :: THSAuxCreate, &
            THSAuxVarInit, &
            THSAuxVarCompute, &
            THSAuxVarCopy, &
            THSAuxDestroy, &
            THSAuxVarStrip, &
            THSAuxVarPerturb, &
            THSPrintAuxVars, &
            THSOutputAuxVars, &
            THSAuxVarDestroy
  
  
contains

! ************************************************************************** !

function THSAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(ths_type), pointer :: THSAuxCreate
  type(ths_type), pointer :: aux
  
  allocate(aux)
  
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)
  
  THSAuxCreate => aux

end function THSAuxCreate
  
! ************************************************************************** !

subroutine THSAuxVarInit(auxvar,allocate_derivative,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: auxvar
  PetscBool :: allocate_derivative
  type(option_type) :: option
  
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%pert = 0.d0
  
  allocate(auxvar%pres(option%nphase + FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%H(option%nphase))
  auxvar%H = 0.d0
  allocate(auxvar%U(option%nphase))
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  allocate(auxvar%kr(option%nphase))
  auxvar%kr = 0.d0
  
  if (allocate_derivative) then
  
  else
    nullify(auxvar%d)
  endif

end subroutine THSAuxVarInit

! ************************************************************************** !

subroutine THSAuxVarCompute(x,ths_auxvar,global_auxvar,material_auxvar, &
                            characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  ! 
  
  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  type(ths_auxvar_type) :: ths_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: natural_id
  
  
  

end subroutine THSAuxVarCompute

! ************************************************************************** !

subroutine THSAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: auxvar,auxvar2
  type(option_type) :: option
  
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmol = auxvar%xmol
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%kr = auxvar%kr
  auxvar%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert


end subroutine THSAuxVarCopy

! ************************************************************************** !

subroutine THSAuxDestroy(aux)
  !
  ! Deallocates a THS auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  use Utility_module, only : DeallocateArray

  implicit none
  
  type(ths_type), pointer :: aux
  
  if(.not.associated(aux)) return
  
  call THSAuxVarDestroy(aux%auxvars)
  call THSAuxVarDestroy(aux%auxvars_bc)
  call THSAuxVarDestroy(aux%auxvars_ss)
  
  call MatrixZeroingDestroy(aux%matrix_zeroing)
  
  if (associated(aux%ths_parameter)) then
    call DeallocateArray(aux%ths_parameter%diffusion_coefficient)
    deallocate(aux%ths_parameter)
  endif
  
  nullify(aux%ths_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine THSAuxDestroy

! ************************************************************************** !

subroutine THSAuxVarPerturb(ths_auxvar,global_auxvar,material_auxvar,&
                            characteristic_curves, natural_id, option)
  !
  ! Calculates auxiliary variables for a perturbed system
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  type(ths_auxvar_type) :: ths_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id
  type(option_type) :: option

end subroutine THSAuxVarPerturb

! ************************************************************************** !

subroutine THSPrintAuxVars(ths_auxvar,global_auxvar,material_auxvar, &
                           natural_id, string, option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Michael Nole
  ! Date: 03/23/20
  !
  
  use Global_Aux_module
  use Material_Aux_class
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: ths_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
                           
end subroutine THSPrintAuxVars

! ************* Local Subroutines ************* !
! ************************************************************************** !

subroutine THSOutputAuxVars1(ths_auxvar, global_auxvar, material_auxvar, &
                             natural_id, string, append, option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Michael Nole
  ! Date: 03/19/20
  !
  
  use Global_Aux_module
  use Material_Aux_class
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: ths_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: append
  type(option_type) :: option
  
                             
end subroutine THSOutputAuxVars1

! ************************************************************************** !

subroutine THSOutputAuxVars2(ths_auxvars, global_auxvars, option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Michael Nole
  ! Date: 03/19/20
  !
  
  use Global_Aux_module
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: ths_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option
  

end subroutine THSOutputAuxVars2
                             
! ************************************************************************** !

subroutine THSAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  implicit none
  
  type(ths_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call THSAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  
  nullify(auxvar)

end subroutine THSAuxVarSingleDestroy

! ************************************************************************** !

subroutine THSAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates an auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  implicit none
  
  type(ths_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call THSAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THSAuxVarArray1Destroy

! ************************************************************************** !

subroutine THSAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates an auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  implicit none
  
  type(ths_auxvar_type), pointer :: auxvars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1,size(auxvars,1)
        call THSAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THSAuxVarArray2Destroy

! ************************************************************************** !

subroutine THSAuxVarStrip(auxvar)
  !
  ! Author: Michael Nole
  ! Date: 03/10/20
  !
  
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(ths_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%den_kg)
  call DeallocateArray(auxvar%xmol)
  call DeallocateArray(auxvar%H)
  call DeallocateArray(auxvar%U)
  call DeallocateArray(auxvar%mobility)
  
  if (associated(auxvar%d)) then
    deallocate(auxvar%d)
    nullify(auxvar%d)
  endif

end subroutine THSAuxVarStrip
  
! ************************************************************************** !

end module THS_Aux_module
