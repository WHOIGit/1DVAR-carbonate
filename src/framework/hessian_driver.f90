!----------------------------------------------------------------------------
!     CVS:$Id: hessian_driver.F90,v 1.7 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

program hessian_driver
!---------------------------------------------------------------------------
! Based on driver.F90, this is the code to compute the hessian using the
! adjoint.  
!---------------------------------------------------------------------------
  
  use admodel_mod, only : simul
  use common_mod, only : bio0,bioparams0,bioparams_default,homedir
  use const, only : c0,c2
  use cost, only : readdata
  use eco_common, only : read_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : nparams_opt,nparams_bio
  use eco_params, only : get_opt_param_names
  use forcing, only : readforcing
  use grid, only : grid_init
  use io, only : driver_init, write_matrix,io_init
  use model_mod, only : bio_to_opt
  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! optparams           optimizable biological parameters (normalized)
  ! adoptparams         ajoint variable to optparams
  ! optparams0          initial parameters
  ! adoptparams0        initial gradient
  ! adoptparams_m       gradient at negative perturbation
  ! adoptparams_p       gradient at positive perturbation
  double precision, dimension(nparams_opt) ::  &
       optparams,adoptparams,optparams0, &
       adoptparams0,adoptparams_m,adoptparams_p

  ! fcost               value of cost function
  ! fcost0              value of cost function at starting point
  ! delta               perturbation
  double precision :: fcost,fcost0,delta

  ! hessian             hessian matrix
  double precision, dimension(nparams_opt,nparams_opt) :: hessian

  ! Variables needed for use with n1qn3 optimizer
  ! izs,rzs             working (dummy) arrays
  integer, dimension(1) :: izs
  double precision, dimension(1) :: rzs
  double precision, dimension(1) :: dzs

  ! iz                  working array for n1qn3
  !integer, dimension(5) :: iz

  ! counters
  integer :: ipar

  integer :: indic=4

  character(len=20), dimension(nparams_opt) :: opt_param_names

  !-------------------------------------------------------------------------
  ! Perform various initializations
  !-------------------------------------------------------------------------
  
  call driver_init
  call read_eco_params(bioparams_default,bioparams0)
  call bio_to_opt(bioparams0,optparams)
  call grid_init
  call readforcing
  call readdata
  call read_bio
  call io_init(.false.)

  !-------------------------------------------------------------------------
  ! Run model in forward/adjoint mode once to set up gradient
  !-------------------------------------------------------------------------
  adoptparams = c0
  indic=4
  call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)
  
  !-------------------------------------------------------------------------
  ! Compute hessian using numerical approximation
  !-------------------------------------------------------------------------

  optparams0=optparams
  fcost0=fcost
  adoptparams0=adoptparams
  delta=1.0e-12

  do ipar=1,nparams_opt
     optparams=optparams0   

     ! optparams0 should be 1, but just in case...
     ! note that this fails if optparams0 = 0

     ! do negative perturbation
     adoptparams_m = c0
     fcost=c0
     optparams(ipar)=optparams0(ipar)-delta
     indic=4
     call simul(indic,nparams_opt,optparams,fcost,adoptparams_m,izs,rzs,dzs)
     
     ! do positive perturbation
     adoptparams_p = c0
     optparams=optparams0
     fcost=c0
     optparams(ipar)=optparams0(ipar)+delta
     indic=4
     call simul(indic,nparams_opt,optparams,fcost,adoptparams_p,izs,rzs,dzs)
     
     ! two sided finite difference approximation for hessian
     hessian(:,ipar) = (adoptparams_p - adoptparams_m)/ &
          (c2*delta)
  end do



  call get_opt_param_names(opt_param_names)
  call write_matrix(trim(homedir)//'/'//'hessian.out',hessian, &
       opt_param_names)


end program hessian_driver



