!----------------------------------------------------------------------------
!     CVS:$Id: adjoint_test.F90,v 1.8 2005/01/31 21:59:43 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

program adjoint_test
!---------------------------------------------------------------------------
! Based on driver.F90, this is the code to run the adjoint of the driver
! forward model
!---------------------------------------------------------------------------
  
  use adj_common, only : dxmin,df1,epsg,impres,imode,niter,nsim
  use adj_common, only : adjoint_init
  use admodel_mod, only : simul
  use common_mod, only : homedir,datdir,bio0,bioparams0,bioparams_default
  use const, only : c10
  use cost, only : readdata
  use eco_common, only : read_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : NumStateVar,nparams_opt,nparams_bio
  use eco_params, only : get_opt_param_names !Luo_g
  use forcing, only : readforcing
  use grid, only : grid_init
  use io, only : driver_init,io_init
  use model_mod, only : bio_to_opt
  use model_mod, only : model
  use m1qn3_mod, only : euclid,ctonbe,ctcabe,m1qn3

  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! optparams           optimizable biological parameters (normalized)
  ! adoptparams         adjoint variable to optparams
  ! optparams0          initial parameter values
  ! adoptparams0        initial gradient
  double precision, dimension(nparams_opt) :: &
       optparams,adoptparams,optparams0,adoptparams0
  ! fcost               value of cost function
  ! fcost0              value of cost function at starting point
  ! delta               perturbation
  double precision :: fcost,fcost0,delta

  ! Variables needed for use with n1qn3 optimizer
  ! izs,rzs             working (dummy) arrays
  integer, dimension(1) :: izs
  double precision, dimension(1) :: rzs
  double precision, dimension(1) :: dzs

  ! iz                  working array for n1qn3
  !integer, dimension(5) :: iz

  ! nrz                 see n1qn3 documentation
  integer, parameter :: nrz=4*nparams_opt + 7*(2*nparams_opt +1)
  ! rz
  !double precision, dimension(nrz) :: rz

  ! normg               norm of the gradient
  double precision :: normg

  integer :: indic=4
  
  character(len=*), parameter :: & !Luo_gradient
      fname_gradient_out='gradient.csv' 
     
  integer, parameter :: g_un=50 !Luo_gradient
  
  character(len=20), dimension(nparams_opt) :: opt_param_names !Luo_gradient
  
  integer :: iopt !Luo_gradient

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
  call io_init(.false.)

  !-------------------------------------------------------------------------
  ! Run model in forward mode once to set up gradient
  !-------------------------------------------------------------------------
  call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)
  
  !-------------------------------------------------------------------------
  ! Run tests
  !-------------------------------------------------------------------------
  open(unit=g_un,file=trim(homedir)//'/'//trim(fname_gradient_out), &
            status='replace')  !Luo_g
  call get_opt_param_names(opt_param_names) !Luo_g
  do iopt=1, nparams_opt !Luo_g
      write(g_un,'(a20,",",g15.3,",")') opt_param_names(iopt), adoptparams(iopt)
  end do
  close(g_un) !Luo_g
  
  optparams0=optparams
  fcost0=fcost
  adoptparams0=adoptparams
  open(unit=g_un,file=trim(homedir)//'/adjtest_each_opt.dat', &
            status='replace')
  do iopt = 1, nparams_opt
      delta=1d-3
      normg = sqrt(adoptparams0(iopt)*adoptparams0(iopt))
      write (6,*) opt_param_names(iopt), 'adjoint_test:: normg ',normg
      write (g_un,*) iopt, ' ==============================='
      write (g_un,*) opt_param_names(iopt), 'adjoint_test:: normg ',normg
      do while (delta > 1d-9)
         optparams = optparams0
         optparams(iopt) = optparams0(iopt)-delta*adoptparams0(iopt)/normg
         indic=4
         call simul(indic,nparams_opt,optparams,fcost,adoptparams,izs,rzs,dzs)
         ! print out diagnostic based on deltaJ = dJ/dp deltap
         write (6,*) opt_param_names(iopt), 'adjoint_test:: d, dJ:d, J',delta, &
              (fcost-fcost0)/(-delta*(adoptparams0(iopt)*adoptparams0(iopt))/normg),fcost
         write (g_un,*) opt_param_names(iopt), 'adjoint_test:: d, dJ:d, J',delta, &
              (fcost-fcost0)/(-delta*(adoptparams0(iopt)*adoptparams0(iopt))/normg),fcost
         delta = delta * 0.1
      end do
  end do
  close(g_un)

end program adjoint_test



