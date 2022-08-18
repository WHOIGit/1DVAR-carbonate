!----------------------------------------------------------------------------
!     CVS:$Id: adjoint_driver.F90,v 1.9 2005/03/07 15:02:35 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

program adjoint_driver
!---------------------------------------------------------------------------
! Based on driver.F90, this is the code to run the adjoint of the driver
! forward model
!---------------------------------------------------------------------------

  use adj_common, only : dxmin,df1,epsg,impres,imode,niter,nsim
  use adj_common, only : adjoint_init
  use common_mod, only : homedir,datdir,bio0,bioparams0,bioparams_default
  use common_mod, only : ecopar_fname, iotempflag
  use const, only : c0,c1,c10
  use cost, only : readdata
  use eco_common, only : read_eco_params,write_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : NumStateVar,nparams_opt,nparams_bio
  use eco_params, only : get_opt_param_names
  use forcing, only : readforcing
  use grid, only : grid_init
  use io, only : io_init,io_final,driver_init,write_matrix,results_un
  use model_mod, only : bio_to_opt,opt_to_bio
  use admodel_mod, only : simul
  use model_mod, only : model
  use m1qn3_mod, only : euclid,ctonbe,ctcabe,m1qn3
  use hessian_mod, only : hessian_init, hessian_inv
  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! bioparams           optimizable biological parameters (various units)
  double precision, dimension(nparams_opt) :: optparams,adoptparams
  double precision, dimension(nparams_bio) :: bioparams
  ! fcost               value of cost function
  double precision :: fcost


  ! Variables needed for use with n1qn3 optimizer

  ! izs,rzs             working (dummy) arrays
  integer, dimension(1) :: izs
  double precision, dimension(1) :: rzs
  double precision, dimension(1) :: dzs

  ! stdout              unit number for stdout
  integer, parameter :: stdout=6

  ! iz                  working array for n1qn3
  integer, dimension(5) :: iz

  ! nrz                 see n1qn3 documentation
  integer, parameter :: nrz=4*nparams_opt+7*(2*nparams_opt+1)

  ! rz
  double precision, dimension(nrz) :: rz

  integer :: indic=4
  
  integer :: omode
  
  logical :: reverse=.false.

  character(len=20), dimension(nparams_opt) :: opt_param_names

  !-------------------------------------------------------------------------
  ! Perform various initializations
  !-------------------------------------------------------------------------
  !! [HK] save error/output from derivs_mod.f90 
  !! open(unit = 100, file = 'derivs_mod_int_output')

  call driver_init
  call adjoint_init
  call read_eco_params(bioparams_default,bioparams0)
! write out initial and default values for record keeping
  call write_eco_params(bioparams_default,trim(ecopar_fname)//'.default')
  call write_eco_params(bioparams0,trim(ecopar_fname)//'.initial')
  call bio_to_opt(bioparams0,optparams)
  call grid_init
  call readforcing
  call readdata
  call read_bio
  call io_init(.true.)
  call hessian_init

  !-------------------------------------------------------------------------
  ! Run model in forward/adjoint mode once to set up gradient
  !-------------------------------------------------------------------------
  call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)

  ! write out initial results
  write(results_un,'(102(1x,1PG19.12,:))') 0,fcost,optparams

  !-------------------------------------------------------------------------
  ! Run optimization
  !-------------------------------------------------------------------------
  call m1qn3(simul,euclid,ctonbe,ctcabe,nparams_opt,optparams,fcost, &
       adoptparams,dxmin,df1,epsg,impres,stdout,imode,omode,niter, &
       nsim, iz,rz,nrz,reverse,indic,izs,rzs,dzs)
  
  ! rerun model one last time to generate output
  iotempflag =.true.
  call model(optparams,fcost)

  call io_final
  call opt_to_bio(optparams,bioparams)
  call write_eco_params(bioparams,trim(ecopar_fname)//'.new')
  
  call get_opt_param_names(opt_param_names)
  call write_matrix(trim(homedir)//'/'//'invhessian.out',hessian_inv, &
       opt_param_names)


end program adjoint_driver

