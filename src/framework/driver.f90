!----------------------------------------------------------------------------
!     CVS:$Id: driver.F90,v 1.14 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
program driver
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

  use common_mod, only : homedir,datdir,bio0,bioparams0,bioparams_default,iotempflag
  use cost, only : readdata
  use eco_common, only : read_eco_params
  use eco_derivs, only : read_bio
  use eco_params, only : NumStateVar,nparams_opt,nparams_bio
  use forcing, only : readforcing,setup_forcing
  use grid, only : grid_init
  use io, only : io_init
  use io, only : io_final,driver_init
  use model_mod, only : bio_to_opt
  use model_mod, only : model

  implicit none

  !-------------------------------------------------------------------------
  ! Local variables
  !-------------------------------------------------------------------------

  ! optparams           optimizable biological parameters (various units)
  double precision, dimension(nparams_opt) :: optparams
  ! fcost               value of cost function
  double precision :: fcost

  !-------------------------------------------------------------------------
  ! Perform various initializations
  !-------------------------------------------------------------------------
  iotempflag = .true.
  call driver_init
  call read_eco_params(bioparams_default,bioparams0)
  call bio_to_opt(bioparams0,optparams)
  call grid_init
  call readforcing
  call readdata
  call read_bio
  call io_init(.false.)

  !-------------------------------------------------------------------------
  ! Run model in forward mode
  !-------------------------------------------------------------------------
  call model(optparams,fcost)

  !-------------------------------------------------------------------------
  ! Clean up
  !-------------------------------------------------------------------------
  call io_final

  
end program driver


