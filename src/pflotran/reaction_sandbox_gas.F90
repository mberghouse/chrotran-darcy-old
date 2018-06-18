module Reaction_Sandbox_Gas_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_gas_type
  contains
    procedure, public :: ReadInput => GasRead
    procedure, public :: Setup => GasSetup
    procedure, public :: Evaluate => GasReact
    procedure, public :: Destroy => GasDestroy
  end type reaction_sandbox_gas_type

  public :: GasCreate

contains

! ************************************************************************** !

function GasCreate()
  ! 
  ! Allocates gas reaction object.
  ! 
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  ! 

  implicit none
  
  class(reaction_sandbox_gas_type), pointer :: GasCreate

  allocate(GasCreate)
  nullify(GasCreate%next)  
      
end function GasCreate

! ************************************************************************** !

subroutine GasRead(this,input,option)
  ! 
  ! Reads input deck for gas reaction parameters (if any)
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_gas_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
  enddo
  
end subroutine GasRead

! ************************************************************************** !

subroutine GasSetup(this,reaction,option)
  ! 
  ! Sets up the gas reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_gas_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

end subroutine GasSetup

! ************************************************************************** !

subroutine GasReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_gas_type) :: this  
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
  PetscReal :: L_water

  ! Key gas variables:
  ! rt_auxvar%%gas_pp(igas): partial pressure of gas in bars [10^5 Pa]
  !   convert to gas concentration by multiplying by 
  !     1.d5 / (IDEAL_GAS_CONSTANT * (global_auxvar%temp+273.15d0)*1.d3
  ! rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase): block in Jacobian 
  !   representing chemistry in the grid cell.  The block is square
  !   measuring (number of aqueous components + number of immobile 
  !   components) on a side. Aqueous components are numbered in order 
  !   specified in PRIMARY_SPECIES block followed by IMMOBILE_SPECIES
  !     units: kg water / L gas (mixture)
  !     icomp,jcomp = derivative of icomp total concentation [mol/L gas] 
  !                   with respect to jcomp aqueous free ion concentration 
  !                   [mol/kg water]
  !     iphase = 2

  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives 
  !   or a numerical approximation unless always running with 
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of 
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error 
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in Gas'
  !   call printErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.  
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations 
  !                                 [mol/L water] for iphase = 1
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3] 
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3] 
  !     global_auxvar%sat(iphase) - saturation 
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]                     
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.
  
  ! Unit of the residual must be in moles/second
  ! global_auxvar%sat(iphase) = saturation of cell
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! always subtract contribution from residual
!  Residual(this%species_id) = Residual(this%species_id) - &
!    this%rate_constant * &  ! 1/sec
!    L_water * & ! L water
    ! rt_auxvar%total(this%species_id,iphase) = species total component 
    !   concentration
!    rt_auxvar%total(this%species_id,iphase) ! mol/L water
    
  
  
  if (compute_derivative) then

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
!    Jacobian(this%species_id,this%species_id) = &
!    Jacobian(this%species_id,this%species_id) + &
!      this%rate_constant * & ! 1/sec
!      L_water * & ! L water
                  ! kg water/L water
      ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) = 
      !   derivative of total component concentration with respect to the
      !   free ion concentration of the same species.
!      rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) 

  endif
  
end subroutine GasReact

! ************************************************************************** !

subroutine GasDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  implicit none
  
  class(reaction_sandbox_gas_type) :: this  

end subroutine GasDestroy

end module Reaction_Sandbox_Gas_class
