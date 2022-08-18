!----------------------------------------------------------------------------
!     CVS:$Id: param_srch.F90,v 1.6 2005/01/19 19:20:45 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

program param_srch
  !--------------------------------------------------------------------------
  ! This code does a search in the region near the starting point, stopping
  ! when a drop in cost function is seen.
  !--------------------------------------------------------------------------
 
  ! TODO: clean these up
  use adj_common, only : dxmin,df1,epsg,impres,mode,niter,nsim
  use adj_common, only : adjoint_init
  use admodel_mod, only : simul
  use common_mod, only : homedir,datdir,bio0,bioparams0,bioparams_default
  use common_mod, only : ecopar_fname
  use const, only : c1,c10
  use cost, only : readdata
  use eco_common, only : read_eco_params,write_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : NumStateVar,nparams_opt,nparams_bio
  use forcing, only : readforcing
  use grid, only : grid_init
  use grid, only : NumDepth
  use io, only : driver_init, io_init
  use kinds_mod, only : dbl_kind,int_kind
  !use light, only : light_init
  use model_mod, only : bio_to_opt,opt_to_bio
  use model_mod, only : model
  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! bioparams           optimizable biological parameters (various units)
  real(kind=dbl_kind), dimension(nparams_opt) :: &
       optparams,adoptparams,optparams0,adoptparams0,lnsrch_optparams
  real(kind=dbl_kind), dimension(nparams_bio) :: bioparams

  ! fcost               value of cost function
  real(kind=dbl_kind) :: fcost,fcost0,delta,delta0,mult


  ! Variables needed for use with n1qn3 optimizer

  ! izs,rzs             working (dummy) arrays
  integer(kind=int_kind), dimension(1) :: izs
  real(kind=dbl_kind), dimension(1) :: rzs
  real(kind=dbl_kind), dimension(1) :: dzs

  ! iz                  working array for n1qn3
  integer(kind=int_kind), dimension(5) :: iz

  ! nrz                 see n1qn3 documentation
  integer(kind=int_kind), parameter :: nrz=4*nparams_opt + 7*(2*nparams_opt +1)
  
  ! rz
  real(kind=dbl_kind), dimension(nrz) :: rz

  ! normg
  real(kind=dbl_kind) :: normg

  ! array of random numbers
  real(kind=dbl_kind), dimension(nparams_opt) :: randnums

  integer(kind=int_kind) :: i, indic=4

!  integer lnsrch_dir   ! search direction for line search
!  logical lnsrch_success, lnsrch_done, lnsrch_found, lnsrch_pass2
!  real lnsrch_fcost

  !-------------------------------------------------------------------------
  ! Perform various initializations
  !-------------------------------------------------------------------------
  call driver_init
  call adjoint_init
  call read_eco_params(bioparams_default,bioparams0)
  call bio_to_opt(bioparams0,optparams)
  call grid_init
  call readforcing
  call readdata
  call read_bio
  call io_init(.true.)

  !-------------------------------------------------------------------------
  ! Run model in forward mode once to set up gradient
  !-------------------------------------------------------------------------
  call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)
  
  !-------------------------------------------------------------------------
  ! Run tests
  !-------------------------------------------------------------------------

  optparams0=optparams
  ! fcost0 is the target, epsg (from n1qn3 code) is the search criteria
  ! make more restrictive, to cut down # searches.
  ! fcost0=fcost * (c1-epsg)
  fcost0=fcost * (c1-1.0e-4)
  adoptparams0=adoptparams
  delta0 = 1e-6
  delta=delta0
  mult=2.0
  normg = sqrt(sum(adoptparams*adoptparams))
  write(6,*) 'param_srch: initial cost ',fcost0
  write (6,*) 'param_srch:: normg ',normg
 
  write(6,*) 'param_srch: trying random search'
  mult=1.02

  do i=1,2   ! number of attempts - if first fails, second is likely also
     delta = 1e-4
     call random_seed()
     do while (delta < 1.0_dbl_kind)
        ! note: delta is scaled by optparams0 (the starting location)
        call random_number(randnums)
        !     write(6,*) randnums
        randnums = (randnums * -2.0_dbl_kind) + 1.0_dbl_kind
        optparams = optparams0+delta*randnums*optparams0

        ! attempt to reduce roundoff problem
        call opt_to_bio(optparams,bioparams)
        call bio_to_opt(bioparams,optparams)
        
        indic=4
        call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)
        write (6,*) 'param_srch:: d, J, dJ',delta,fcost,fcost-fcost0
        if (fcost < fcost0) then
           write (6,*) 'param_srch: found lower cost function (random)'
           call opt_to_bio(optparams,bioparams)
           call write_eco_params(bioparams,trim(ecopar_fname)//'.new')
           stop
        end if
        delta = delta * mult
     end do
  end do

  
  write(6,*) 'param_srch: no lower cost function found'
  stop 1
    
  end program param_srch



