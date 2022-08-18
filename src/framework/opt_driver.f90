!----------------------------------------------------------------------------
!     CVS:$Id: opt_driver V1.0 05/09/2008 01:00 LUO
!     CVS:$Name:  $
!----------------------------------------------------------------------------
program opt_driver
!----------------------------------------------------------------
!
!     Marjorie Friedrichs 03/03/03
!
!     Subroutines
!       readdata.f:      reads in data for computation of cost fn
!       readforcing.f:   reads in forcing files
!       readinit.f:      reads in initial conditions
!                                 ecosystem parameters
!       model.f:         calls physical and biological subroutines
!          biosub.f:     general biological processes
!              derivs.f: specific ecosystem model equations
!          horizadv.f:   applies horizontal advection
!          tridag.f:     tridiagonal solver for diffusion
!          sink.f:       applies detrital sinking
!       output.f:        outputs model equivalents of data
!       calccost.f:      outputs cost function
!
!-----------------------------------------------------------------
! 

!---------------------------------------------------------------------------
! Revised April 2004, January 2005 - Jeff Dusenberry
!
!---------------------------------------------------------------------------

  use common_mod, only : homedir,datdir,bio0,bioparams0,bioparams_default
  use cost, only : readdata
  use eco_common, only : read_eco_params,write_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : NumStateVar,nparams_opt,nparams_bio
  use forcing, only : readforcing,setup_forcing
  use grid, only : grid_init
  use io, only : io_init
  use io, only : io_init,io_final,driver_init,write_matrix,results_un
  use model_mod, only : bio_to_opt,opt_to_bio
  use model_mod, only : model,CALFUN
  use newuoa_mod, only: NEWUOA
  use common_mod, only : ecopar_fname
  use opt_common, only : opt_init,NPT,IPRINT,MAXFUN,RHOBEG,RHOEND,W

  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! optparams           optimizable biological parameters (various units)
  double precision, dimension(nparams_opt) :: optparams
  double precision, dimension(nparams_bio) :: bioparams
  ! fcost               value of cost function
  double precision :: fcost

  !-------------------------------------------------------------------------
  ! Perform various initializations
  !-------------------------------------------------------------------------
  call driver_init
  call read_eco_params(bioparams_default,bioparams0)
  ! write out initial and default values for record keeping
  call write_eco_params(bioparams_default,trim(ecopar_fname)//'.default')
  call write_eco_params(bioparams0,trim(ecopar_fname)//'.initial')
  call bio_to_opt(bioparams0,optparams)
  call grid_init
  call readforcing
  call readdata
  call read_bio
  call opt_init
  call io_init(.true.)

  !-------------------------------------------------------------------------
  ! Run model in forward mode
  !-------------------------------------------------------------------------
  call CALFUN(nparams_opt,optparams,fcost)
  ! write out initial results
  write(results_un,'(102(1x,1PG19.12,:))') 0,fcost,optparams
  
  call NEWUOA (nparams_opt,NPT,optparams,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
  
  ! rerun model one last time to generate output
  call model(optparams,fcost)

  call io_final
  call opt_to_bio(optparams,bioparams)
  call write_eco_params(bioparams,trim(ecopar_fname)//'.new')
  
end program opt_driver
