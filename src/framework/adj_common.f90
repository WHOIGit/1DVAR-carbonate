!----------------------------------------------------------------------------
!     CVS:$Id: adj_common.F90,v 1.3 2005/01/31 16:25:49 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module adj_common
  implicit none

  ! save everything
  save

  ! dxmin               resolution in x for l-infinity norm
  ! df1                 expected decrease in cost function first iteration
  ! epsg                stopping criterion
  double precision :: dxmin,df1,epsg
  
  ! mode                specified input and output modes of m1qn3
  ! niter               maximum/actual number of iterations
  ! nsim                maximum/actual number of simulations
  ! impres              controls level of debugging info printed
  integer :: niter,nsim,impres
  
  integer :: imode(2)


contains

  subroutine adjoint_init
    
    use common_mod, only : adjpar_fname
    implicit none

    integer, parameter :: iun=46

    imode(1) = 0 ! DIS mode
    imode(2) = 0 ! cold start

    namelist /adjoint_parms_nml/ &
         df1, &
         dxmin, &
         epsg, &
         impres, &
         niter, &
         nsim

    !-----------------------------------------------------------------------
    !   default namelist settings
    !-----------------------------------------------------------------------
    df1=1.0
    dxmin=1e-8
    epsg=1e-9
    impres=5
    niter=1500
    nsim=2500

    !-----------------------------------------------------------------------
    !   read in parameter values
    !-----------------------------------------------------------------------
    open(unit=iun, file=trim(adjpar_fname), status='old')
    read(unit=iun, nml=adjoint_parms_nml)
    close(unit=iun)

  end subroutine adjoint_init

end module adj_common
