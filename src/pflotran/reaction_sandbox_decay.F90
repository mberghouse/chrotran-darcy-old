module Reaction_Sandbox_Decay_class

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_decay_type
    ! Aqueous species
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
  contains
    procedure, public :: ReadInput => DecayRead
    procedure, public :: Setup => DecaySetup
    procedure, public :: Evaluate => DecayReact
  end type reaction_sandbox_decay_type

  public :: DecayCreate

contains

! ************************************************************************** !

function DecayCreate()
  ! 
  ! Allocates decay reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  implicit none
  
  class(reaction_sandbox_decay_type), pointer :: DecayCreate

  allocate(DecayCreate)
  DecayCreate%species_name = ''
  DecayCreate%species_id = UNINITIALIZED_INTEGER
  DecayCreate%rate_constant = UNINITIALIZED_DOUBLE
  nullify(DecayCreate%next)  
      
end function DecayCreate

! ************************************************************************** !

subroutine DecayRead(this,input,option)
  ! 
  ! Reads input deck for decay reaction parameters (if any)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/22/18
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_decay_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,DECAY')
    call StringToUpper(word)

    select case(trim(word))
      case('SPECIES_NAME')
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'species_name', &
                           'CHEMISTRY,REACTION_SANDBOX,DECAY')
      case('RATE_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,DECAY')
        ! Read the units
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          ! If units do not exist, assume default units of 1/s which are the
          ! standard internal PFLOTRAN units for this rate constant.
          input%err_buf = 'REACTION_SANDBOX,DECAY,RATE CONSTANT UNITS'
          call InputDefaultMsg(input,option)
        else
          ! If units exist, convert to internal units of 1/s
          internal_units = 'unitless/sec'
          this%rate_constant = this%rate_constant * &
            UnitsConvertToInternal(word,internal_units,option)
        endif
      case default
        call InputKeywordUnrecognized(word, &
                       'CHEMISTRY,REACTION_SANDBOX,DECAY',option)
    end select
  enddo

end subroutine DecayRead


! ************************************************************************** !

subroutine DecaySetup(this,reaction,option)
  ! 
  ! Sets up the decay reaction with hardwired parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_decay_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)
  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'Uninitialized rate constant for species ' // &
      trim(this%species_name) // ' in decay reaction sandbox.'
  endif
      
end subroutine DecaySetup

! ************************************************************************** !

subroutine DecayReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/15
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_decay_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: kg_water               ! kg water
  
  PetscReal :: rate
  PetscReal :: derivative 
  
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  kg_water = porosity*liquid_saturation*volume*global_auxvar%den_kg(iphase)
  
  ! mol/sec = 1/sec * mol/kg water * kg water
  rate = -1.d0 * this%rate_constant * rt_auxvar%pri_molal(this%species_id) * &
         kg_water 
  ! kg water/sec = 1/sec * kg water
  derivative = -1.d0 * this%rate_constant * kg_water
  
  Residual(this%species_id) = Residual(this%species_id) - rate

  if (compute_derivative) then
    Jacobian(this%species_id,this%species_id) = &
      Jacobian(this%species_id,this%species_id) - derivative
  endif
  
end subroutine DecayReact

end module Reaction_Sandbox_Decay_class
