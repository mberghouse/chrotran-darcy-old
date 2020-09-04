module Sensibility_Analysis_module

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Richards_Common_module
  use Richards_Aux_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Material_Aux_class
  
  implicit none

  private
  
  PetscInt, parameter :: PERMEABILITY = 1
  PetscInt, parameter :: POROSITY = 2
  
  ! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: RichardsPermeabilitySensibility, &
            RichardsPorositySensibility

contains

! ************************************************************************** !

subroutine RichardsPermeabilitySensibility(realization,ierr)

  use Realization_Subsurface_class
  
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  call RichardsSensibility(realization,PERMEABILITY,ierr)

end subroutine

! ************************************************************************** !

subroutine RichardsPorositySensibility(realization,ierr)

  use Realization_Subsurface_class
  
  Mat :: J
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  call RichardsSensibility(realization,POROSITY,ierr)

end subroutine

! ************************************************************************** !

subroutine RichardsSensibility(realization,ivar,ierr)
  ! 
  ! Computes derivative of the residual according to the permeability at 
  ! each grid cell
  ! Most of the subroutines are taken from RichardsJacobian
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 
  ! Note:
  ! Inline surface flow not considered

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module
  use Debug_module
  use Discretization_module

  implicit none

  type(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: J_mat_type
  PetscViewer :: viewer
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  call PetscLogEventBegin(logging%event_r_jacobian,ierr);CHKERRQ(ierr)
  
  option => realization%option
  
  !prepare J matrix
  J_mat_type = MATBAIJ
  call DiscretizationCreateJacobian(realization%discretization, &
                                    NFLOWDOF, J_mat_type, J, option)

  call MatSetOptionsPrefix(J,"sensibility_",ierr);CHKERRQ(ierr)
  
  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  !call RichardsSensibilityInternalConn(J,realization,ivar,ierr)
  !call RichardsSensibilityBoundaryConn(J,realization,ivar,ierr)
  !call RichardsSensibilitySourceSink(J,realization,ivar,ierr)
  !update here when porosity ok
  !call RichardsSensibilityAccumulation(J,realization,ivar,ierr)

  select case(ivar)
    case(PERMEABILITY)
      call DebugWriteFilename(realization%debug,string,'K_sensibility','', &
                              richards_ts_count,richards_ts_cut_count, &
                              richards_ni_count)
      call DebugCreateViewer(realization%debug,string,option,viewer)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call DebugViewerDestroy(realization%debug,viewer)
    case(POROSITY)
      call DebugWriteFilename(realization%debug,string, &
                              'Porosity_sensibility', '', &
                              richards_ts_count,richards_ts_cut_count, &
                              richards_ni_count)
    case default
      call PrintErrMsg(option, "Wrong value of ivar in RichardsSensibility")
  end select
  
  !destroy J
  call MatDestroy(J,ierr);CHKERRQ(ierr)
    
end subroutine

! ************************************************************************** !

subroutine RichardsSensibilityInternalConn(A,realization,ivar,ierr)
  ! 
  ! Computes the interior flux terms of the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 
       
  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  !use Coupler_module
  !use Field_module
  !use Debug_module
  use Material_Aux_class
  !use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscErrorCode :: ierr

  PetscInt :: icc_up,icc_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  !PetscInt :: region_id_up, region_id_dn
  PetscInt :: istart_up, istart_dn, istart

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: unit_z(3) = [0.d0,0.d0,1.d0]

!  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  !type(field_type), pointer :: field
  !type(region_type), pointer :: region
  !type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  !field => realization%field
  !material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
   
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (associated(patch%aux%Richards%matrix_buffer)) then
      call MatrixBufferZero(patch%aux%Richards%matrix_buffer)
    else
      patch%aux%Richards%matrix_buffer => MatrixBufferCreate()
      call MatrixBufferInit(A,patch%aux%Richards%matrix_buffer,grid)
    endif
  endif
#endif

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (abs(dot_product(cur_connection_set%dist(1:3,iconn),unit_z)) < &
            0.99d0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)
      
      call RichardsFluxSensibility(rich_auxvars(ghosted_id_up), &
                 global_auxvars(ghosted_id_up), &
                 material_auxvars(ghosted_id_up), &
                 rich_auxvars(ghosted_id_dn), &
                 global_auxvars(ghosted_id_dn), &
                 material_auxvars(ghosted_id_dn), &
                 cur_connection_set%area(iconn), &
                 cur_connection_set%dist(-1:3,iconn),&
                 option,&
                 patch%characteristic_curves_array(icc_up)%ptr, &
                 patch%characteristic_curves_array(icc_dn)%ptr, &
                 Jup,Jdn,ivar)
      
      if (local_id_up > 0) then

#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_up,ghosted_id_up,Jup(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_up,ghosted_id_dn,Jdn(1,1))
        else
#endif
          istart_up = (ghosted_id_up-1)*option%nflowdof + 1
          istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1

          call MatSetValuesLocal(A,1,istart_up-1,1,istart_up-1, &
                                        Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValuesLocal(A,1,istart_up-1,1,istart_dn-1, &
                                        Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif

      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_dn,ghosted_id_dn,Jdn(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_dn,ghosted_id_up,Jup(1,1))
        else
#endif
          istart_up = (ghosted_id_up-1)*option%nflowdof + 1
          istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1

          call MatSetValuesLocal(A,1,istart_dn-1,1,istart_dn-1, &
                                        Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValuesLocal(A,1,istart_dn-1,1,istart_up-1, &
                                        Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

end subroutine RichardsSensibilityInternalConn

! ************************************************************************** !

subroutine RichardsFluxSensibility(rich_auxvar_up,global_auxvar_up, &
                                  material_auxvar_up, & 
                                  rich_auxvar_dn,global_auxvar_dn, &
                                  material_auxvar_dn, &
                                  area, dist, &
                                  option, &
                                  characteristic_curves_up, &
                                  characteristic_curves_dn, &
                                  Jup,Jdn, ivar)
  ! 
  ! Computes the sensibility of the internal flux terms
  ! For K or poro, just change the value in the numerical derivative with a 
  ! select statement
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 
  use Option_module 
  use Characteristic_Curves_module
  use Material_Aux_class
  use Connection_module
  
  implicit none
  
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy, area, dist(-1:3)
  type(characteristic_curves_type) :: characteristic_curves_up
  type(characteristic_curves_type) :: characteristic_curves_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ivar

  PetscInt :: ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_up, rich_auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert_up, material_auxvar_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  v_darcy = 0.D0  

  call GlobalAuxVarInit(global_auxvar_pert_up,option)
  call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
  call MaterialAuxVarInit(material_auxvar_pert_up,option)
  call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
  call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
  call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
  call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
  call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
  call MaterialAuxVarCopy(material_auxvar_up,material_auxvar_pert_up,option)
  call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn,option)
  x_up(1) = global_auxvar_up%pres(1)
  x_dn(1) = global_auxvar_dn%pres(1)
  call RichardsFlux(rich_auxvar_up,global_auxvar_up,material_auxvar_up, &
                    rich_auxvar_dn,global_auxvar_dn,material_auxvar_dn, &
                    area, dist, &
                    option,v_darcy,res)
  ideriv = 1
  
  
  select case (ivar)
    case (PERMEABILITY)
      ! TODO
   
!    pert_up = x_up(ideriv)*perturbation_tolerance
  pert_up = max(dabs(x_up(ideriv)*perturbation_tolerance),0.1d0)
  if (x_up(ideriv) < option%reference_pressure) pert_up = -1.d0*pert_up
!    pert_dn = x_dn(ideriv)*perturbation_tolerance
  pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
  if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
  x_pert_up = x_up
  x_pert_dn = x_dn
  x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
  x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn

    case (POROSITY)
    
    case default
      
  end select
  
  call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                             global_auxvar_pert_up, &
                             material_auxvar_pert_up, &
                             characteristic_curves_up, &
                             -999, &
                             PETSC_TRUE,option)
  call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                             global_auxvar_pert_dn, &
                             material_auxvar_pert_dn, &
                             characteristic_curves_dn, &
                             -999, &
                             PETSC_TRUE,option)
  call RichardsFlux(rich_auxvar_pert_up,global_auxvar_pert_up, &
                    material_auxvar_pert_up, &
                    rich_auxvar_dn,global_auxvar_dn, &
                    material_auxvar_dn, &
                    area, dist, &
                    option,v_darcy,res_pert_up)
  call RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                    material_auxvar_up, &
                    rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                    material_auxvar_pert_dn, &
                    area, dist, &
                    option,v_darcy,res_pert_dn)
  
  J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
  J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
  Jup = J_pert_up
  Jdn = J_pert_dn
  call GlobalAuxVarStrip(global_auxvar_pert_up)
  call GlobalAuxVarStrip(global_auxvar_pert_dn)    
  call MaterialAuxVarStrip(material_auxvar_pert_up)
  call MaterialAuxVarStrip(material_auxvar_pert_dn)    

end subroutine RichardsFluxSensibility

! ************************************************************************** !

subroutine RichardsSensibilityBoundaryConn(A,realization,ivar,ierr)
  ! 
  ! Computes the boundary flux terms of the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  !use Debug_module
  use Material_Aux_class
  use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscErrorCode :: ierr

  PetscInt :: icc_up,icc_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, region_id, i
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: istart_up, istart_dn, istart
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(region_type), pointer :: region
  type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icc_dn = patch%cc_id(ghosted_id) 

      call RichardsBCFluxSensibility(boundary_condition%flow_condition%itype, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     rich_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     rich_auxvars(ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     option, &
                     patch%characteristic_curves_array(icc_dn)%ptr, &
                     Jdn, ivar)
      Jdn = -Jdn

#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jdn(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jdn, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine RichardsSensibilityBoundaryConn

! ************************************************************************** !

subroutine RichardsBCFluxSensibility(ibndtype,auxvars, &
                                    rich_auxvar_up,global_auxvar_up, &
                                    rich_auxvar_dn,global_auxvar_dn, &
                                    material_auxvar_dn, &
                                    area,dist,option, &
                                    characteristic_curves_dn, &
                                    Jdn, ivar)
  ! 
  ! Computes numerically the derivatives of the boundary flux
  ! terms for the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 
  use Option_module
  use Characteristic_Curves_module
  use EOS_Water_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
  PetscReal :: area
  ! dist(-1) = fraction_upwind
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(-1:3)
  type(characteristic_curves_type) :: characteristic_curves_dn
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ivar
  
  PetscReal :: v_darcy
  PetscInt :: pressure_bc_type

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_dn, rich_auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  class(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  PetscErrorCode :: ierr

  v_darcy = 0.d0
  Jdn = 0.d0 

  call GlobalAuxVarInit(global_auxvar_pert_up,option)
  call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
  allocate(material_auxvar_pert_up,material_auxvar_pert_dn)
  call MaterialAuxVarInit(material_auxvar_pert_up,option)  
  call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
  call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
  call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
  call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
  call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
  call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_up, &
                          option)
  call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn, &
                          option)
  x_up(1) = global_auxvar_up%pres(1)
  x_dn(1) = global_auxvar_dn%pres(1)
  ideriv = 1
  if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
    x_up(ideriv) = x_dn(ideriv)
  endif
  call RichardsBCFlux(ibndtype,auxvars, &
                      rich_auxvar_up,global_auxvar_up, &
                      rich_auxvar_dn,global_auxvar_dn, &
                      material_auxvar_dn, &
                      area,dist,option,v_darcy,res)
  
  select case(ivar)
    case (PERMEABILITY)
      ! TODO
    case (POROSITY)
    
    case default
    
  end select

  if (pressure_bc_type == ZERO_GRADIENT_BC) then
    x_pert_up = x_up
  endif
  ideriv = 1
!    pert_dn = x_dn(ideriv)*perturbation_tolerance    
  pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
  if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
  x_pert_dn = x_dn
  x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
  x_pert_up = x_up
  if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
    x_pert_up(ideriv) = x_pert_dn(ideriv)
  endif
  
  
  call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                             global_auxvar_pert_dn, &
                             material_auxvar_pert_dn, &
                             characteristic_curves_dn, &
                             -999, &
                             PETSC_TRUE,option)
  call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                             global_auxvar_pert_up, &
                             material_auxvar_pert_up, &
                             characteristic_curves_dn, &
                             -999, &
                             PETSC_TRUE,option)
  call RichardsBCFlux(ibndtype,auxvars, &
                      rich_auxvar_pert_up,global_auxvar_pert_up, &
                      rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                      material_auxvar_pert_dn, &
                      area,dist,option,v_darcy,res_pert_dn)
  J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
  Jdn = J_pert_dn
  call GlobalAuxVarStrip(global_auxvar_pert_up)
  call GlobalAuxVarStrip(global_auxvar_pert_dn)   
  call MaterialAuxVarStrip(material_auxvar_pert_up)
  call MaterialAuxVarStrip(material_auxvar_pert_dn)
  deallocate(material_auxvar_pert_up,material_auxvar_pert_dn)

end subroutine RichardsBCFluxSensibility

! ************************************************************************** !

subroutine RichardsSensibilitySourceSink(A,realization,ivar,ierr)
  ! 
  ! Computes the accumulation and source/sink terms of
  ! the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscErrorCode :: ierr

  PetscReal :: qsrc
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: flow_pc
  PetscViewer :: viewer
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: ukvr, Dq, dphi, v_darcy
  Vec, parameter :: null_vec = tVec(0)
  character(len=MAXSTRINGLENGTH) :: string

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    if (source_sink%flow_condition%itype(1)/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%itype(1)/=HET_MASS_RATE_SS .and. &
       source_sink%flow_condition%itype(1)/=WELL_SS) &
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS,SCALED_MASS_RATE_SS,HET_MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O* &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_VOL_RATE_SS)
          Jup(1,1) = -source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* &
                    rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(WELL_SS) ! production well, SK 12/19/13
          ! if node pessure is lower than the given extraction pressure, 
          ! shut it down
          !  well parameter explanation
          !   1. well status. 1 injection; -1 production; 0 shut in!
          !   2. well factor [m^3],  the effective permeability [m^2/s]
          !   3. bottomhole pressure:  [Pa]
          !   4. max pressure: [Pa]
          !   5. min pressure: [Pa]
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)
    
          ! production well (well status = -1)
          if (dabs(well_status + 1.D0) < 1.D-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor 
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.D0) then ! outflow only
                ukvr = rich_auxvars(ghosted_id)%kvr
                if (ukvr < 1.e-20) ukvr = 0.D0
                v_darcy = 0.D0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  Jup(1,1) = -1.d0*(-Dq*rich_auxvars(ghosted_id)%dkvr_dp*dphi* &
                             global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*1.d0*global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*dphi*rich_auxvars(ghosted_id)%dden_dp)
                endif
              endif
            endif
          endif 
      end select
#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jup(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup,ADD_VALUES, &
                               ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
    enddo
    source_sink => source_sink%next
  enddo

  !call RichardsSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
  !                       global_auxvars,rich_auxvars,option)

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_srcsink','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (patch%aux%Richards%inactive_cells_exist) then
      call MatrixBufferZeroRows(patch%aux%Richards%matrix_buffer, &
                                patch%aux%Richards%matrix_zeroing% &
                                  n_inactive_rows, &
                                patch%aux%Richards%matrix_zeroing% &
                                  inactive_rows_local_ghosted)
    endif
    call MatrixBufferSetValues(A,patch%aux%Richards%matrix_buffer)
  endif
#endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

! zero out isothermal and inactive cells
#ifdef BUFFER_MATRIX
  if (.not.option%use_matrix_buffer) then
#endif
    if (patch%aux%Richards%inactive_cells_exist) then
      qsrc = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%Richards%matrix_zeroing% &
                              n_inactive_rows, &
                            patch%aux%Richards%matrix_zeroing% &
                              inactive_rows_local_ghosted, &
                            qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                            ierr);CHKERRQ(ierr)
    endif
#ifdef BUFFER_MATRIX
  endif
#endif

end subroutine RichardsSensibilitySourceSink

! ************************************************************************** !

subroutine RichardsSensibilityAccumulation(A,realization,ivar,ierr)
  ! 
  ! Computes the accumulation terms of the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Debug_module
  use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscErrorCode :: ierr

  PetscInt :: local_id, ghosted_id, region_id
  PetscInt :: istart

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(region_type), pointer :: region
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  !type(inlinesurface_auxvar_type), pointer :: inlinesurface_auxvars(:)
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  if (.not.option%steady_state) then

    ! Accumulation terms ------------------------------------
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RichardsAccumSensibility(rich_auxvars(ghosted_id), &
           global_auxvars(ghosted_id), &
           material_auxvars(ghosted_id), &
           option, &
           patch%characteristic_curves_array( &
           patch%cc_id(ghosted_id))%ptr, &
           Jup, ivar)

#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
             ghosted_id,Jup(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup, &
             ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
    enddo

#if 0
    if (option%inline_surface_flow) then
      do region_id = 1, region%num_cells
        local_id = region%cell_ids(region_id)
        ghosted_id = grid%nL2G(local_id)         
        if (patch%imat(ghosted_id) <= 0) cycle
        call InlineSurfaceAccumulationJac(inlinesurface_auxvars(region_id), &
             material_auxvars(ghosted_id),option,Jup)
        istart = (ghosted_id-1)*option%nflowdof + 1
        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup, &
             ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    endif
#endif

  endif

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_accum','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine RichardsSensibilityAccumulation

! ************************************************************************** !

subroutine RichardsAccumSensibility(rich_auxvar,global_auxvar, &
                                   material_auxvar, &
                                   option, &
                                   characteristic_curves, &
                                   J,ivar)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the sensibility
  ! 
  ! Author: Moise Rousseau
  ! Date: 09/03/2020
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialAuxVarInit, &
                                 MaterialAuxVarCopy, &
                                 MaterialAuxVarStrip, &
                                 MaterialCompressSoil
  
  implicit none

  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(characteristic_curves_type) :: characteristic_curves
  PetscReal :: J(option%nflowdof,option%nflowdof)
  PetscInt :: ivar
     
  PetscInt :: ispec 
  PetscReal :: vol_over_dt
  PetscReal :: tempreal
  PetscReal :: derivative

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  vol_over_dt = material_auxvar%volume/option%flow_dt
  
  select case(ivar)
    case (PERMEABILITY)
      ! TODO it's null
    case (POROSITY)
      ! I may compute this analytically
      ! accumulation term units = dkmol/dp
      J(1,1) = (material_auxvar%dporosity_dp*global_auxvar%sat(1)* &
                global_auxvar%den(1) + &
                (global_auxvar%sat(1)*rich_auxvar%dden_dp + &
                 rich_auxvar%dsat_dp*global_auxvar%den(1)) * &
                material_auxvar%porosity) * &
                vol_over_dt
    case default
    
  end select
  
  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert,option)  
    call MaterialAuxVarInit(material_auxvar_pert,option)  
    call RichardsAuxVarCopy(rich_auxvar,rich_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)
    x(1) = global_auxvar%pres(1)
    call RichardsAccumulation(rich_auxvar,global_auxvar,material_auxvar, &
                              option,res)
    ideriv = 1
    pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
    x_pert = x
    if (x_pert(ideriv) < option%reference_pressure) pert = -1.d0*pert
    x_pert(ideriv) = x_pert(ideriv) + pert
    
    call RichardsAuxVarCompute(x_pert(1),rich_auxvar_pert,global_auxvar_pert, &
                               material_auxvar_pert, &
                               characteristic_curves, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsAccumulation(rich_auxvar_pert,global_auxvar_pert, &
                              material_auxvar_pert, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)  
    call MaterialAuxVarStrip(material_auxvar_pert)  
  endif
   
end subroutine RichardsAccumSensibility

! ************************************************************************** !
    
end module
