!----------------------------------------------------------------------------
!     CVS:$Id: cost.F90,v 1.19 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module cost
  use common_mod, only : nsite,nz_max,nt_max
  use const, only : c1
  use kinds_mod, only : dbl_kind,int_kind
  implicit none

  ! keep everything private, unless specified otherwise
  private
  ! save everything
  save

  public calccost
  public output
  public readdata
  public setup_cost

  ! ndat_max         maximum number data points for cost data
  integer(kind=int_kind), parameter, public :: ndat_max=500

  ! indices for cost function variables
  ! xcnit            nitrate (mmolN/m3)
  ! xczoo            zooplankton (mmolC/m3)
  ! xcchl            phytoplankton/chlorophyll (mgChl/m3)
  ! xcppr            primary productivity (mgC/m3/day)
  ! xcbac           Bacterial Biomass (mmolC/m3) !Luo_MB
  ! xcbpr            Bacterial Production (mmolC/m3/day) !Luo_MB
  ! xcdoc           semi-labile DOC Concentration (mmolC/m3) !Luo_MB
  integer(kind=int_kind), parameter, public :: xcnit=1,xczoo=2,xcchl=3,xcppr=4, &
                                                                                xcbac=5,xcbpr=6,xcdoc=7 !Luo_MB

  ! ncnsv            number of cost function variables
  ! ncnsv_sed        number of export-type cost function variables
  !                  NOTE:  ncnsv_sed has limited usage, changing
  !                         from value of 1 will require some recoding
  integer(kind=int_kind), parameter, public :: ncnsv=7 !Luo_MB
  integer(kind=int_kind), parameter :: ncnsv_sed=1

  ! cost_un          unit number for cost output values
  integer(kind=int_kind), parameter, public :: cost_un=40

  ! weights for cost function - must be in order of xc index params above
  ! units are inverse of units listed above.
  real(kind=dbl_kind), allocatable:: w_cost_glob(:,:)
  real(kind=dbl_kind), allocatable :: w_cost_sed_glob(:)

  real(kind=dbl_kind), dimension(ncnsv), public :: w_cost
  real(kind=dbl_kind), public :: w_cost_sed

  ! input files are formed as (note double use of prefix):
  ! datdir/dat_prefix/dat_prefix_fname_cost_suffix.dat
  ! for corresponding output files:
  ! homedir/dat_prefix/dat_prefix_fname_cost_suffix.out
  character(len=*), dimension(ncnsv), parameter :: fname_cost_suffix=(/ &
       'N   ', &
       'Z   ', &
       'P   ', &
       'PrPr', &
       'BAC ', &
       'BPr ', &
       'sDOC' /) !Luo_MB

  ! fname_sed_suffix         suffix for sediment trap observational datafile
  ! fname_var_suffix         suffix for cost function inverse weights filename
  character(len=*), parameter :: &
       fname_sed_suffix='ST', &
       fname_var_suffix='vars.dat'

  ! ndat_tot         number of data points in observations
  integer(kind=int_kind), allocatable :: ndat_tot_glob(:)
  integer(kind=int_kind), public :: ndat_tot

  ! time_cdat        time of observational data
  ! z_cdat           depth of observational data
  ! a_loc            observations
  real(kind=dbl_kind), dimension(ndat_max,ncnsv) :: &
       time_cdat,z_cdat,a_loc
  real(kind=dbl_kind), allocatable :: a_norm_glob(:,:)
  real(kind=dbl_kind), dimension(ncnsv), public :: a_norm

  ! sediment trap flux variables
  ! btime_csed       begin time for trap deployments
  ! z_csed
  ! a_csed
  ! dtime_csed       time of deployment
  ! etime_csed       ending time for trap deployments
  real(kind=dbl_kind), allocatable :: &
       btime_csed_glob(:,:),z_csed_glob(:,:),etime_csed_glob(:,:)
  real(kind=dbl_kind), dimension(ndat_max), public :: btime_csed,z_csed, &
       etime_csed
  real(kind=dbl_kind), allocatable :: &
       a_csed_glob(:,:),dtime_csed_glob(:,:)
  real(kind=dbl_kind), dimension(ndat_max), public :: a_csed
  real(kind=dbl_kind), dimension(ndat_max) :: dtime_csed

  ! ndat_sed         number of data points for trap data
  integer(kind=int_kind), allocatable:: ndat_sed_glob(:)
  integer(kind=int_kind), public :: ndat_sed
  ! a_norm_sed
  real(kind=dbl_kind), allocatable :: a_norm_sed_glob(:)
  real(kind=dbl_kind), public :: a_norm_sed


  ! mask             mask (against observations) for cost function
  ! a_cdat           observations on model grid
  real(kind=dbl_kind), allocatable:: &
       mask_glob(:,:,:,:),a_cdat_glob(:,:,:,:)

  real(kind=dbl_kind), dimension(nz_max,ncnsv,nt_max), public :: mask,a_cdat
  

contains

  !-------------------------------------------------------------------------
  ! This subroutine reads in data to be included in cost function
  ! computation
  !-------------------------------------------------------------------------
  subroutine readdata
    use common_mod, only : datdir,nsite,dat_prefix_glob
    use const, only : c0,c1,SecPerDay
    use forcing, only : cdays_glob
    use grid, only : zifc_glob,nz_glob,delt_glob
    implicit none


    ! iterators,counters
    integer(kind=int_kind) id,istep,ixc,iz,idat,isite
    integer(kind=int_kind), dimension(nz_max,nt_max) :: ncount

    integer(kind=int_kind), dimension(ncnsv) :: ndat

    ! unit numbers
    integer(kind=int_kind), parameter :: iun=36

    allocate (w_cost_glob(ncnsv,nsite))
    allocate (w_cost_sed_glob(nsite))
    allocate (ndat_tot_glob(nsite))
    allocate (a_norm_glob(ncnsv,nsite))
    allocate ( &
       btime_csed_glob(ndat_max,nsite), &
       z_csed_glob(ndat_max,nsite), &
       etime_csed_glob(ndat_max,nsite), &
       a_csed_glob(ndat_max,nsite), &
       dtime_csed_glob(ndat_max,nsite) )
    allocate (ndat_sed_glob(nsite))
    allocate (a_norm_sed_glob(nsite))
    allocate ( &
       mask_glob(nz_max,ncnsv,nt_max,nsite), &
       a_cdat_glob(nz_max,ncnsv,nt_max,nsite) )

    ! load in data files for scalars in cost function

    siteloop: do isite=1,nsite

       do ixc=1,ncnsv
          open(unit=iun, &
               file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
               trim(dat_prefix_glob(isite))//'_'// &
               trim(fname_cost_suffix(ixc))//'.dat', &
               status='old')
          read(iun,*)ndat(ixc)
          ! verify ndat large enough
          call check_ndat(ndat(ixc))
          do id=1,ndat(ixc)
             read(iun,*)time_cdat(id,ixc),z_cdat(id,ixc),a_loc(id,ixc)
          enddo
          close(iun)
       enddo

       ! load in data file for sediment flux
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'// &
            trim(fname_sed_suffix)//'.dat', &
            status='old')
       read(iun,*)ndat_sed_glob(isite)
       call check_ndat(ndat_sed_glob(isite))
       do id=1,ndat_sed_glob(isite)
          read(iun,*)btime_csed_glob(id,isite),z_csed_glob(id,isite), &
               dtime_csed_glob(id,isite),a_csed_glob(id,isite)
          etime_csed_glob(id,isite) = btime_csed_glob(id,isite)+ &
               dtime_csed_glob(id,isite)
       end do
       close(iun)

       ! read in cost function weights (inverse)
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_var_suffix), &
            status='old')
       ! read in the inverse weights
       do ixc=1,ncnsv
          read(iun,*)w_cost_glob(ixc,isite)
       end do
       read(iun,*)w_cost_sed_glob(isite)
       close(iun)
       ! invert to weight
!       w_cost_glob(:,isite) = c1/w_cost_glob(:,isite)
!       w_cost_sed_glob(isite) = c1/w_cost_sed_glob(isite)


       ! create mask and map data to model grid
       mask_glob(:,:,:,isite) = c0
       a_cdat_glob(:,:,:,isite) = c0


       ! bin data into model grid, averaging as necessary
       varloop: do ixc=1,ncnsv
          ncount=0
          dataloop: do idat=1,ndat(ixc)

!             iz=int(z_cdat(idat,ixc)/dz)+1            
             iz = 1
             do while (z_cdat(idat,ixc) > zifc_glob(iz+1,isite) .and. &
                  iz < nz_glob(isite))
                iz = iz+1
             end do
             ! check to see if cost data is below grid, disregard if so.
             if (z_cdat(idat,ixc) > zifc_glob(nz_glob(isite)+1,isite)) then
                cycle dataloop
             end if

             istep=int((time_cdat(idat,ixc)-cdays_glob(1,isite))/ &
                  delt_glob(isite)*SecPerDay)+1
             mask_glob(iz,ixc,istep,isite)=c1
             a_cdat_glob(iz,ixc,istep,isite)= &
                  a_cdat_glob(iz,ixc,istep,isite)+a_loc(idat,ixc)
             ncount(iz,istep)=ncount(iz,istep)+1

          end do dataloop

          where (ncount .gt. 0)
             a_cdat_glob(:,ixc,:,isite)=a_cdat_glob(:,ixc,:,isite)/ncount
          endwhere
       end do varloop

       ! compute normalizations
       a_norm_glob(:,isite) = c1/((ncnsv+ncnsv_sed)*ndat)
       a_norm_sed_glob(isite) = c1/((ncnsv+ncnsv_sed)*ndat_sed_glob(isite))
       ndat_tot_glob(isite) = 1

    end do siteloop

  end subroutine readdata

  
  !-------------------------------------------------------------------------
  ! check_ndat checks ndat_test against ndat_max to ensure enough storage
  ! allocated.
  !-------------------------------------------------------------------------
  subroutine check_ndat(ndat_test)

    implicit none

    ! ndat_test          number of data points being read in
    integer(kind=int_kind), intent(in) :: ndat_test

    if (ndat_test .gt. ndat_max) then
       write(6,*) 'ERROR: ndat_max not large enough, need to recompile'
       write(6,*) 'ndat_max = ',ndat_max,' but need ndat = ',ndat_test
       stop
    endif

  end subroutine check_ndat

  
  !-------------------------------------------------------------------------
  ! print out standard output files of variables at data points
  ! corresponding to locations of cost function data
  !-------------------------------------------------------------------------
  subroutine output(isite,cpred,csed_pred)
    use common_mod, only : homedir,dat_prefix_glob
    use const, only : c0
    use forcing, only : cdays
    use grid, only : zmid,nz_glob,nt_glob
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------
    integer(kind=int_kind) :: isite
    real(kind=dbl_kind), dimension(:,:,:) :: cpred
    real(kind=dbl_kind), dimension(:) :: csed_pred

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    ! counters
    integer(kind=int_kind) ixc,it,iz,ised

    ! unit number
    integer(kind=int_kind), parameter :: iun=91

    ! output cost variable model equivalents
    do ixc=1,ncnsv
       open(unit=iun, &
            file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'// &
            trim(fname_cost_suffix(ixc))//'.out')
       write(iun,'(i10)') nint(sum(mask(:,ixc,:)))
       do it=1,nt_glob(isite)
          do iz=1,nz_glob(isite)
             if (mask(iz,ixc,it) .gt.c0) then
                write(iun,'(f8.3,f8.1,f15.9)')cdays(it), &
                     zmid(iz),cpred(iz,ixc,it)
             endif
          enddo
       enddo
       close(iun)
    enddo
    
    ! write out model equivalents of sediment trap flux
    open(unit=iun, &
         file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
         trim(dat_prefix_glob(isite))//'_'//trim(fname_sed_suffix)//'.out')
    do ised=1,ndat_sed
       write(iun,'(f8.3,f8.3,f10.3,f15.9)') &
            btime_csed(ised),z_csed(ised), &
            dtime_csed(ised),csed_pred(ised)
    enddo
    close(iun)

  end subroutine output


subroutine calccost(fcost,cpred,csed_pred)
  use const, only : c0
  use common_mod, only : ioflag,dat_prefix_loc
  implicit none

  !----------------------------------------------------------------------  
  ! Arguments
  !----------------------------------------------------------------------

  ! fcost         cost function value
  ! cpred         predicted values of cost variables
  ! csed_pred     predicted values of sediment flux
  real(kind=dbl_kind), intent(out) :: fcost

  real(kind=dbl_kind), dimension(:,:,:) :: cpred
  real(kind=dbl_kind), dimension(:) :: csed_pred

  !----------------------------------------------------------------------  
  ! Local variables
  !----------------------------------------------------------------------

  ! f_cost        cost function (each quantity)
  ! f_cost_sed    cost function for sediment trap fluxes
  real(kind=dbl_kind), dimension(ncnsv) :: f_cost
  real(kind=dbl_kind) :: f_cost_sed

  ! counters/indices
  integer(kind=int_kind) :: ixc,ised,iz,it

  !----------------------------------------------------------------------
  ! Calculate cost function for scalars
  !----------------------------------------------------------------------
  f_cost=c0
  do ixc=1,ncnsv
!     f_cost(ixc) = f_cost(ixc)+a_norm(ixc)*sum(mask(:,ixc,:)* &
!          (w_cost(ixc)*(a_cdat(:,ixc,:)-cpred(:,ixc,:)))**2)
     do  it=1,nt_max
            do iz = 1,nz_max
                  if (mask(iz,ixc,it).GT.0) then
                       f_cost(ixc) =  f_cost(ixc) + mask(iz,ixc,it)* &
                                       (w_cost(ixc)*(1-cpred(iz,ixc,it)/a_cdat(iz,ixc,it)))**2
                  end if
             end do
      end do
     f_cost(ixc) =   f_cost(ixc) * a_norm(ixc)
!      f_cost(ixc) = a_norm(ixc)*sum(mask(:,ixc,:)* &
!          (w_cost(ixc)*(1-cpred(:,ixc,:)/a_cdat(:,ixc,:)))**2)
!     f_cost(ixc) = a_norm(ixc)*sum(mask(:,ixc,:)* &
!          (w_cost(ixc)*(a_cdat(:,ixc,:)-cpred(:,ixc,:)))**2)
  enddo
  
  !----------------------------------------------------------------------
  ! Compute the contribution from the sediment trap results
  !----------------------------------------------------------------------

!  f_cost_sed = a_norm_sed*sum((w_cost_sed* &
!       (a_csed(1:ndat_sed)-csed_pred(1:ndat_sed)))**2)

  f_cost_sed = c0

  if(ncnsv_sed.gt.0) then
      do ised = 1,ndat_sed
    	 f_cost_sed = f_cost_sed + a_norm_sed*(w_cost_sed* &
          (a_csed(ised)-csed_pred(ised))/a_csed(ised))**2
 !   	 f_cost_sed = f_cost_sed + a_norm_sed*(w_cost_sed* &
 !         (a_csed(ised)-csed_pred(ised)))**2

! use logarithm
!     f_cost_sed = f_cost_sed + a_norm_sed*(w_cost_sed* &
!          (log(a_csed(ised))-log(csed_pred(ised))))**2
     end do
  endif


  !----------------------------------------------------------------------
  ! Normalize the cost values based on number of data points
  !----------------------------------------------------------------------

  fcost=(sum(f_cost) + f_cost_sed)/ndat_tot
  ! debugging only 
  !fcost=f_cost_sed/ndat_tot


    write(*,*)trim(dat_prefix_loc)//' '//'Normalized cost = ', &
       sngl(fcost)
    write(*,*)trim(dat_prefix_loc)//' '//'N cost = ', &
         sngl(f_cost(xcnit)/ndat_tot)
    write(*,*)trim(dat_prefix_loc)//' '//'Z cost = ', &
         sngl(f_cost(xczoo)/ndat_tot)
    write(*,*)trim(dat_prefix_loc)//' '//'P cost = ', &
         sngl(f_cost(xcchl)/ndat_tot)
    write(*,*)trim(dat_prefix_loc)//' '//'PP cost =', &
         sngl(f_cost(xcppr)/ndat_tot)
    write(*,*)trim(dat_prefix_loc)//' '//'BAC cost =', &
         sngl(f_cost(xcbac)/ndat_tot) !Luo_MB
    write(*,*)trim(dat_prefix_loc)//' '//'BPr cost =', &
         sngl(f_cost(xcbpr)/ndat_tot) !Luo_MB
    write(*,*)trim(dat_prefix_loc)//' '//'sDOC cost =', &
         sngl(f_cost(xcdoc)/ndat_tot) !Luo_MB
    write(*,*)trim(dat_prefix_loc)//' '//'ST cost =', &
         sngl(f_cost_sed/ndat_tot)
    if (ioflag) then
       write(cost_un,'(6(1x,1PG13.6))') fcost,f_cost/ndat_tot, &
            f_cost_sed/ndat_tot
    end if

  end subroutine calccost



  subroutine setup_cost(isite)
    implicit none

    ! isite           location index
    integer(kind=int_kind) :: isite
    
    ndat_sed = ndat_sed_glob(isite)
    ndat_tot = ndat_tot_glob(isite)
    a_norm = a_norm_glob(:,isite)
    a_norm_sed = a_norm_sed_glob(isite)
    a_csed = a_csed_glob(:,isite)
    w_cost = w_cost_glob(:,isite)
    w_cost_sed = w_cost_sed_glob(isite)
    mask = mask_glob(:,:,:,isite)
    a_cdat = a_cdat_glob(:,:,:,isite)
    btime_csed = btime_csed_glob(:,isite)
    dtime_csed = dtime_csed_glob(:,isite)
    etime_csed = etime_csed_glob(:,isite)
    z_csed = z_csed_glob(:,isite)

  end subroutine setup_cost


end module cost
