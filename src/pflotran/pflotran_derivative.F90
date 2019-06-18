program pflotran_derivative
  
  use Option_module
  use General_Derivative_module
  use EOS_module
    
  implicit none

#include "petsc/finclude/petscsys.h"

  class(option_type), pointer :: option
  
  option => OptionCreate()
  call OptionInitMPI(option)  
  call OptionInitPetsc(option)  
  call EOSInit()
  call GeneralDerivativeDriver(option)
  call OptionFinalize(option)
  
end program pflotran_derivative
