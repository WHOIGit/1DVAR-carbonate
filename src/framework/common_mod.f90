!----------------------------------------------------------------------------
!     CVS:$Id: common_mod.F90,v 1.18 2005/01/31 21:59:43 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module common_mod
  use const, only : c0
  use eco_params, only : NumStateVar,NumDiagVar,nparams_bio
  implicit none

  ! save everything
  save

  !-----------------------------------------------------------------------
  ! parameters
  !-----------------------------------------------------------------------

  integer, parameter :: &
       nsite_max=19    , & ! maximum number of sites
       nz_max=25      , & ! maximum number of depth levels
       nt_max=12000        ! maximum number of time steps

  integer :: &
       nsite              ! number of sites
  ! homedir
  ! datdir
  character(len=80) :: homedir,datdir

  ! input file names
  !
  ! ecopar_fname       filename for ecosystem parameters
  ! adjpar_fname       filename for optimizer parameters
  character(len=132) :: ecopar_fname,adjpar_fname,optpar_fname
  
  ! aeoflux            Aeolian flux forcing
  double precision, dimension(NumStateVar,nt_max) :: aeoflux = c0

  ! bio0               state variable initial conditions
  double precision, dimension(nz_max+1,NumStateVar) :: bio0

  ! bbcnsv             bottom boundary conditions for each state variable
  !                      use < -1 to extrapolate
  double precision, dimension(NumStateVar,nt_max) :: bbcnsv
  
  integer, dimension(NumStateVar) :: flag_bbcnsv


  !double precision, dimension(nz_max,numstatevar,nt_max) :: horiz_adv
  double precision, dimension(:,:,:), allocatable :: horiz_adv
  logical, dimension(numstatevar) :: lhoriz_adv

  ! bioparams0         initial biological parameter values
  ! bioparams_default  default biological parameters, used for normalization
  double precision, dimension(nparams_bio) :: bioparams0,bioparams_default

  ! ioflag            if .true. write to costs.out, params.out
  logical :: ioflag = .true.
  ! iotempflag            if .true. write to iotemp.dat
  logical :: iotempflag = .false.
      
  integer, parameter :: iotemp_un=38

  character(len=12), dimension(nsite_max) :: dat_prefix_glob

  character(len=12) :: dat_prefix_loc


contains
  
  ! set up any site-specific common variables
  subroutine setup_common(isite)
    implicit none

    integer :: isite

    dat_prefix_loc = dat_prefix_glob(isite)
        
    if (iotempflag) then
         open(unit=iotemp_un, &
            file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'// &
            'physics.dat', &
                status='replace')
    end if

  end subroutine setup_common


end module common_mod
