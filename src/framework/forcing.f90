
!----------------------------------------------------------------------------
!     CVS:$Id: forcing.F90,v 1.16 2005/02/16 15:27:40 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

module forcing
  use common_mod, only : nt_max,nz_max
  implicit none


  ! all fields private, unless specified otherwise
  private

  ! save everything
  save

! zref_kv         ref. depth below which Kv values set to background (m)
! kv_bgd          background diffusivity () 
  double precision, parameter :: zref_kv=250.0,kv_bgd=1.0e-5

!!  ! nz_forc         number of depth levels in forcing file
!!  integer, parameter :: nz_forc=8
  
  ! filenames for forcing input files
  character(len=*), parameter :: &
       fname_par_suffix='PAR', &
       fname_mld_suffix='mld', &
       fname_w_suffix='w', &
       fname_kv_suffix='Kv', &
       fname_t_suffix='T'

  !-------------------------------------------------------------------------
  ! "global" arrays of forcing fields
  !-------------------------------------------------------------------------

  ! TODO: remove public attribute


  ! qi_glob            surface irradiance ()
  ! dml_glob           mixed layer depth (m)
  ! cdays_glob         time points in model space (days)
  double precision, allocatable, target, public ::  &
       qi_glob(:,:)       , & !
       dml_glob(:,:)      , & !
       cdays_glob(:,:)

  ! TODO: why was cdays_glob public???
!  double precision, dimension(nt_max,nsite_max), target :: cdays_glob

  ! forcing fields (profiles)
  double precision, allocatable, target ::  &
       Tdat_glob(:,:,:)     , & ! temperature
       HANi_glob(:,:,:)     , & ! Nitrate horizontal advection term
       HAAm_glob(:,:,:)     , & ! Ammonium horizontal advection term
       HAFe_glob(:,:,:)         ! Iron horizontal advection term


! JAD - Note:  rkz changed from tracer levels to level interfaces 10/4/2006

  double precision, allocatable, target :: &
       wvel_glob(:,:,:)     , & ! vertical velocity
       rkz_glob(:,:,:)          ! diffusivity

  !-------------------------------------------------------------------------
  ! site-specific aliases of forcing fields
  !-------------------------------------------------------------------------

  ! dml                mixed layer depth (m)
  ! cdays              time points in model space (days)
  ! qi                 surface irradiance ()
!  double precision, dimension(nt_max), public :: dml,cdays,qi

  double precision, dimension(:), pointer, public ::  &
       dml           , & ! mixed layer depth (m)
       cdays         , & ! time points in model space (days)
       qi                ! surface irradiance ()

  ! forcing fields (profiles)

  double precision, dimension(:,:), pointer, public ::  &
       Tdat          , & ! temperature
       rkz           , & ! diffusivity
       wvel              ! vertical velocity

  public readforcing
  public setup_forcing

contains

  !-------------------------------------------------------------------------
  ! Setup the site specific forcing fields 
  !-------------------------------------------------------------------------
  subroutine setup_forcing(isite)
    use grid, only : nt,nz
    implicit none

    ! NOTE:  this routine assumes that the grid information in grids_mod
    ! is correct for the current site.  Need to call setup_grid first.

    ! isite           location index
    integer :: isite
    
    qi    => qi_glob   (1:nt,isite)
    dml   => dml_glob  (1:nt,isite)
    cdays => cdays_glob(1:nt,isite)

    Tdat  => Tdat_glob(1:nz  ,1:nt,isite)
    rkz   => rkz_glob (1:nz+1,1:nt,isite)
    wvel  => wvel_glob(1:nz+1,1:nt,isite)

  end subroutine setup_forcing
  

  !-------------------------------------------------------------------------
  !     This subroutine reads in the forcing files: PAR
  !     vertical velocity, mixed layer depth, vertical diffusivity,
  !     and temperature.
  !-------------------------------------------------------------------------
  subroutine readforcing
    use common_mod, only : datdir, dat_prefix_glob,nsite,nt_max,nz_max
    use const, only : SecPerDay,c0,c2
    use grid, only : zifc_glob,nz_glob,zmid_glob,nt_glob,delt_glob,zmid_glob
    implicit none

    !-----------------------------------------------------------------------
    ! Local Variables
    !-----------------------------------------------------------------------
    
    ! work variables
    double precision, dimension(nz_max+1) ::  swvel
    double precision :: day,adml

    ! counters and indices
    integer :: it,it2,iz,idx,isite

    double precision :: nct
    
    
    ! t_forc          temperature on forcing grid
    double precision, dimension(nz_max,nt_max) :: t_forc
    
    double precision, dimension(nz_max+1,nt_max) :: kv_forc

    ! unit number
    integer, parameter :: iun=19

    ! Steps per day
    integer, dimension(nsite) :: StepPerDay !Luo
    
    allocate ( &
       qi_glob(nt_max,nsite), & 
       dml_glob(nt_max,nsite), & 
       cdays_glob(nt_max,nsite) )
    allocate (  &
       Tdat_glob(nz_max,nt_max,nsite)     , & ! temperature
       HANi_glob(nz_max,nt_max,nsite)     , & ! Nitrate horizontal advection term
       HAAm_glob(nz_max,nt_max,nsite)     , & ! Ammonium horizontal advection term
       HAFe_glob(nz_max,nt_max,nsite)  )       ! Iron horizontal advection term
    allocate ( &
       wvel_glob(nz_max+1,nt_max,nsite)    , & ! vertical velocity
       rkz_glob(nz_max+1,nt_max,nsite) )         ! diffusivity
       
    ! read in forcing fields for each site
    siteloop: do isite=1,nsite
       StepPerDay = SecPerDay/delt_glob(isite)  !Luo
       ! read in PAR
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_par_suffix)// &
            '.dat', &
            status='old')
       it = 1
       read(iun,*)cdays_glob(it,isite),qi_glob(it,isite) 
       do it=2,nt_glob(isite)
          do it2=1,delt_glob(isite)/1800
              read(iun,*)cdays_glob(it,isite),qi_glob(it,isite)
          enddo
       end do
       close(iun)

       ! read in mixed layer depth
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_mld_suffix)// &
            '.dat', &
            status='old')


       ! NOTE: this approach for reading in forcing files assumes all
       ! forcing files have a 1 day time step.

       nct=c0
       ! read in first time point
       read(iun,*)day,adml
       do it=1,nt_glob(isite)
!          dml_glob(it,isite) = min(60.d0,adml)
          dml_glob(it,isite) = adml
          do while (nct .ge. secperday)
             nct = nct - secperday
             read(iun,*)day,adml
          end do
	  nct=nct+delt_glob(isite)
       enddo
       close(iun)

       ! read in vertical velocity
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_w_suffix)//'.dat', &
            status='old')
       nct=c0
       wvel_glob(:,:,isite) = c0
       read(iun,*)(swvel(iz),iz=1,nz_glob(isite)+1)
       do it=1,nt_glob(isite)
          do iz=1,nz_glob(isite)+1
             wvel_glob(iz,it,isite)=-swvel(iz)
          enddo
          do while (nct .ge. secperday)
             nct = nct - secperday
             read(iun,*)(swvel(iz),iz=1,nz_glob(isite)+1)
          end do
	  nct=nct+delt_glob(isite)
       end do
       close(iun)

       ! read in vertical diffusivity
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_kv_suffix)//'.dat', &
            status='old')
       do it=1,nt_glob(isite)/StepPerDay(isite)  !Luo
          read(iun,*)day,(kv_forc(iz,it),iz=1,nz_glob(isite)+1)
       end do
       close(iun)

       ! put Kv data onto hourly timestep grid
       do iz = 1,nz_glob(isite)+1
          do it=1,nt_glob(isite)/StepPerDay(isite)   !Luo
             do idx=1,StepPerDay(isite)   !Luo
		        if (zmid_glob(iz,isite) .gt. zref_kv) then
                   rkz_glob(iz,(it-1)*StepPerDay(isite)+idx,isite)=kv_bgd
		        !elseif (iz.gt.nz_glob(isite)) then
        	    !   rkz_glob(iz,(it-1)*StepPerDay(isite)+idx,isite)=1.d-3+(kv_forc(iz,it)+ &
                !           ((idx-1)/StepPerDay(isite))*(kv_forc(iz,it+1)-kv_forc(iz,it)))  !Luo
		        else
	               rkz_glob(iz,(it-1)*StepPerDay(isite)+idx,isite)=1.d-6+(kv_forc(iz,it)+ &
                           ((idx-1)/StepPerDay(isite))*(kv_forc(iz,it+1)-kv_forc(iz,it)))  !Luo
		        endif
             end do
          end do
       end do

       ! read in temperature
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_t_suffix)//'.dat', &
            status='old')

       do it=1,nt_glob(isite)/StepPerDay(isite)  !Luo
          do iz=1,nz_glob(isite)
             read(iun,*)day,t_forc(iz,it)
          enddo
       enddo
       close(iun)

       ! put temperature data onto model grid
       do iz = 1,nz_glob(isite)
          do it=1,nt_glob(isite)/StepPerDay(isite)  !Luo
             do idx=1,StepPerDay(isite) 
               Tdat_glob(iz,(it-1)*StepPerDay(isite)+idx,isite)=t_forc(iz,it)  !Luo
             end do
          end do
       end do

    end do siteloop

  end subroutine readforcing

end module forcing
