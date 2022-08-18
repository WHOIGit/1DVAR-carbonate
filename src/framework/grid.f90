!----------------------------------------------------------------------------
!     CVS:$Id: grid.F90,v 1.7 2005/01/31 21:59:43 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module grid
  use common_mod, only : nz_max,nt_max,nsite
  use common_mod, only : datdir,dat_prefix_glob
  implicit none
  
  ! save everything
  save

  ! everything is private unless specified otherwise
  private

  public :: &
       grid_init      , &
       setup_grid


  ! TODO: find way to keep this private
  integer, allocatable, public :: &
       nz_glob(:)        , & ! number of depth levels, all sites
       nt_glob(:)            ! number of time steps, all sites
  
  integer, public :: &
       nz             , & ! number of depth levels, current site
       nt                 ! number of time steps, current site
  
  ! TODO: find way to keep this private
  double precision, allocatable, public :: &
       delt_glob(:)          ! timestep (seconds), all sites
  
  double precision, public :: &
       delt               ! timestep (seconds), current site
       
  integer, allocatable, public :: &
       ntsout_glob(:)        ! number of output intervals, all sites
  
  integer, public :: &
       ntsout                  ! number of output intervals, current site

  ! TODO: find way to keep this private
  double precision, allocatable, target, public :: &
       zmid_glob(:,:)      , & ! depth at midpoint of model levels (+1 extra) (m)
       zifc_glob(:,:)      , & ! depth at interface (starting at surface) (m)
       dzv_glob(:,:)       , & ! distance from midpoint to midpoint (m)
       rdzv_glob(:,:)          ! reciprocal of dzv (1/m)

  ! use aliases for these fields (as opposed to an "array of pointers")
  ! to keep it simpler - JAD.

  ! aliases (per site) for above fields

  double precision, dimension(:), pointer, public :: &
       zmid=>null()   , & ! depth at midpoint of model levels (+1 extra) (m)
       zifc=>null()   , & ! depth at interface (starting at surface) (m)
       dzv=>null()    , & ! distance from midpoint to midpoint (m)
       rdzv=>null()       ! reciprocal of dzv (1/m)

  double precision, allocatable, target :: &
       dzt_glob(:,:)       , & ! level thickness for each level (meters)
       rdzt_glob(:,:)          ! reciprocal of dzt (1/m)

  ! aliases (per site) for above fields

  double precision, dimension(:), pointer, public :: &
       dzt=>null()    , & ! level thickness for each level (meters)
       rdzt=>null()       ! reciprocal of dzt (1/m)

contains
  

  subroutine grid_init
    use const, only : c0,c1,c2,p5
    implicit none

    ! counters, unit numbers
    integer :: iz,iun,isite

    integer :: &
         nz_site     , & ! number of depths at a particular site
         nt_site         ! number of timesteps at a particular site

    double precision, dimension(nz_max) :: &
         dzt_site        ! level thickness (m) at given site
    
    double precision :: &
         delt_site       ! timestep at given site (seconds)
         
    integer :: &
         ntsout_site   ! intervals of output - how many steps will be averaged for one output


    namelist /grid_io_nml/ &
         nz_site, &
         dzt_site, &
         nt_site, &
         delt_site, &
         ntsout_site

    allocate ( &
        nz_glob(nsite), nt_glob(nsite), &
        delt_glob(nsite), ntsout_glob(nsite) )
    allocate (&
       zmid_glob(nz_max+1,nsite)      , & ! depth at midpoint of model levels (+1 extra) (m)
       zifc_glob(nz_max+1,nsite)      , & ! depth at interface (starting at surface) (m)
       dzv_glob(nz_max+1,nsite)      , & ! distance from midpoint to midpoint (m)
       rdzv_glob(nz_max+1,nsite), &
       dzt_glob(nz_max,nsite), &
       rdzt_glob(nz_max,nsite) )         ! reciprocal of dzv (1/m) 

    siteloop: do isite=1,nsite
    
       !----------------------------------------------------------------------
       ! default namelist settings
       !----------------------------------------------------------------------
       nz_site = nz_max
       dzt_site(:) = 10.0   ! 10m default
       nt_site = 10660
       delt_site = 3600.0   ! 3600s default
       ntsout_site = 4 ! default: output at every time step

       ! read in site-specific grid and timestep
       iun = get_lun()
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_grid.in', &
            status='old')

       read(unit=iun, nml=grid_io_nml)
       close(iun)

       ! sanity check
       if (nz_site > nz_max) then
          write(6,*) 'ERROR: nz_max not large enough.'
          write(6,*) 'change value of nz_max and recompile.'
          write(6,*) 'nz_max = ',nz_max,' but need nz = ',nz
          stop
       end if

       if (nt_site > nt_max) then
          write(6,*) 'ERROR: nt_max not large enough.'
          write(6,*) 'change value of nt_max and recompile.'
          write(6,*) 'nt_max = ',nt_max,' but need nt = ',nt_site
          stop
       end if

       ! set global values
       nz_glob(isite) = nz_site
       dzt_glob(1:nz_site,isite) = dzt_site(1:nz_site)
       nt_glob(isite) = nt_site
       delt_glob(isite) = delt_site
       ntsout_glob(isite) = ntsout_site
       
       ! set interface (top of level) and midpoint values
       zifc_glob(1,isite) = c0
       zmid_glob(1,isite) = zifc_glob(1,isite) + dzt_glob(1,isite)/c2
       do iz=2,nz_site
          zifc_glob(iz,isite) = zifc_glob(iz-1,isite) + dzt_glob(iz-1,isite) 
          zmid_glob(iz,isite) = zifc_glob(iz,isite) + dzt_glob(iz,isite)/c2
       end do
       iz=nz_site+1
       zifc_glob(iz,isite) = zifc_glob(iz-1,isite) + dzt_glob(iz-1,isite)
! fix bug June 2007?       zmid_glob(iz,isite) = zifc_glob(iz,isite) + dzt_glob(iz-1,isite)/c2
	zmid_glob(iz,isite) = zifc_glob(iz,isite)

       ! velocity defined at interface
       dzv_glob(1,isite) = dzt_glob(1,isite)   ! extend across boundary
       dzv_glob(2:nz_site,isite) = p5*(dzt_glob(1:nz_site-1,isite) + &
            dzt_glob(2:nz_site,isite))
       dzv_glob(nz_site+1,isite) = dzt_glob(nz_site,isite)

       ! reciprocals of dzt, dzv
       rdzv_glob(1:nz_site+1,isite) = c1/dzv_glob(1:nz_site+1,isite)
       rdzt_glob(1:nz_site,isite) = c1/dzt_glob(1:nz_site,isite)

    end do siteloop

    !print *,'dzt_glob',dzt_glob
    !print *,'dzv_glob',dzv_glob


  end subroutine grid_init



  ! assign aliases for use on per-site basis
  subroutine setup_grid(isite)
    implicit none

    integer :: isite

    nz = nz_glob(isite)
    nt = nt_glob(isite)
    delt = delt_glob(isite)
    ntsout = ntsout_glob(isite)


    dzt  => dzt_glob (1:nz  ,isite)
    zifc => zifc_glob(1:nz+1,isite)
    zmid => zmid_glob(1:nz+1,isite)
    dzv  => dzv_glob (1:nz+1,isite)
    rdzv => rdzv_glob(1:nz+1,isite)
    rdzt => rdzt_glob(1:nz  ,isite)

  end subroutine setup_grid




  ! TODO: maybe move this to make more general?
  ! get an available unit number
  function get_lun()  result (lun)
    implicit none

    integer :: lun
    logical :: exist, opened
    integer :: iostat

    do lun = 11,99
       inquire (unit=lun,exist=exist,opened=opened,iostat=iostat)
       if (exist .and. (.not. opened) .and. (iostat .eq. 0)) return
    end do

    ! consider writing error message here.
    lun = -1

 end function get_lun




end module grid
