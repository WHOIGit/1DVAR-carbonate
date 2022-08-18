!----------------------------------------------------------------------------
!     CVS:$Id: cost.F90,v 1.19 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module cost
  use common_mod, only : nsite,nz_max,nt_max
  use const, only : c1
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
  integer, parameter, public :: ndat_max=3000

  ! indices for cost function variables
  ! xcnit            nitrate (mmolN/m3)
  ! xczoo            zooplankton (mmolC/m3)
  ! xcphy            phytoplankton Nitrogen Biomass (mmol N/m3)
  !xcchl              phytoplankton chlorophyll a (mg Chl a /m3)
  ! xcppr            primary productivity (mgC/m3/day)
  ! xcbac           Bacterial Biomass (mmolC/m3) !Luo_MB
  ! xcbpr            Bacterial Production (mmolC/m3/day) !Luo_MB
  ! xcdoc           semi-labile DOC Concentration (mmolC/m3) !Luo_MB
    integer, parameter, public :: xcnit=1,xcpo4=2,xczoo=3,xcphy=4,xcchl=5, &
                     xcppr=6, xcbac=7,xcbpr=8,xcdoc=9,xcdon=10, &
                      xcdop=11,xcpoc=12, xcpon=13,xcpop=14  !Luo_MB

  ! ncnsv            number of cost function variables
  ! ncnsv_sed        number of export-type cost function variables
  !                  NOTE:  ncnsv_sed has limited usage, changing
  !                         from value of 1 will require some recoding
  integer, parameter, public :: ncnsv=14 !Luo_MB
  
  integer, parameter, public :: xcstc=1,xcstn=2,xcstp=3
  integer, parameter, public :: ncnsv_sed=3

  ! cost_un          unit number for cost output values
  integer, parameter, public :: cost_un=40

  ! weights for cost function - must be in order of xc index params above
  ! units are inverse of units listed above.
  double precision, allocatable:: w_cost_glob(:,:)
  double precision, allocatable :: w_cost_sed_glob(:,:)

  double precision, dimension(ncnsv), public :: w_cost
  double precision, dimension(ncnsv_sed), public :: w_cost_sed

  ! input files are formed as (note double use of prefix):
  ! datdir/dat_prefix/dat_prefix_fname_cost_suffix.dat
  ! for corresponding output files:
  ! homedir/dat_prefix/dat_prefix_fname_cost_suffix.out
  character(len=*), dimension(ncnsv), parameter :: fname_cost_suffix=(/ &
       'NO3 ', &
       'PO4 ', &
       'Z   ', &
       'PHY ', &
       'CHL ', &
       'PrPr', &
       'BAC ', &
       'BPr ', &
       'sDOC', &
       'sDON', &
       'sDOP', &
       'POC ', &
       'PON ', &
       'POP '/) !Luo_MB
       
  ! fname_sed_suffix         suffix for sediment trap observational datafile
  ! fname_var_suffix         suffix for cost function inverse weights filename
  character(len=*), dimension(ncnsv_sed),parameter :: &
       fname_sed_suffix=(/ &
            'STc', &
            'STn', &
            'STp' /)  
  character(len=*), parameter :: fname_var_suffix='vars.dat'

  ! ndat_tot         number of data points in observations
  integer, allocatable :: ndat_tot_glob(:)
  integer, public :: ndat_tot

  ! time_cdat        time of observational data
  ! z_cdat           depth of observational data
  ! a_loc            observations
  double precision, dimension(ndat_max,ncnsv) :: &
       time_loc,z_loc,a_loc
  double precision, dimension(:,:,:),allocatable ::  a_cdat_glob
  integer, dimension(:,:,:), allocatable :: z_cdat_glob, time_cdat_glob
  double precision, dimension(ncnsv, ndat_max), public ::   a_cdat
  integer, dimension(ncnsv, ndat_max), public :: z_cdat, time_cdat
  double precision, allocatable :: a_norm_glob(:,:)
  double precision, dimension(ncnsv), public :: a_norm

  ! sediment trap flux variables
  ! btime_csed       begin time for trap deployments
  ! z_csed
  ! a_csed
  ! dtime_csed       time of deployment
  ! etime_csed       ending time for trap deployments
  double precision, allocatable :: &
       btime_csed_glob(:,:,:),z_csed_glob(:,:,:),etime_csed_glob(:,:,:)
  double precision, dimension(ncnsv_sed,ndat_max), public :: btime_csed,z_csed, &
       etime_csed
  double precision, allocatable :: &
       a_csed_glob(:,:,:),dtime_csed_glob(:,:,:)
  double precision, dimension(ncnsv_sed,ndat_max), public :: a_csed
  double precision, dimension(ncnsv_sed,ndat_max) :: dtime_csed

  ! ndat_sed         number of data points for trap data
  integer, allocatable:: ndat_glob(:,:) 
  integer, allocatable:: ndat_sed_glob(:,:)
  integer, dimension(ncnsv_sed), public :: ndat_sed
  integer, dimension(ncnsv), public :: ndat
  ! a_norm_sed
  double precision, allocatable :: a_norm_sed_glob(:,:)
  double precision, dimension(ncnsv_sed), public :: a_norm_sed


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
    integer id,istep,ixc,iz,idat,isite
    integer, dimension(nz_max,nt_max) :: ncount

    ! unit numbers
    integer, parameter :: iun=36

    allocate (w_cost_glob(ncnsv,nsite))
    allocate (w_cost_sed_glob(ncnsv_sed,nsite))
    allocate (ndat_tot_glob(nsite))
    allocate (a_norm_glob(ncnsv,nsite))
    allocate ( &
       btime_csed_glob(ncnsv_sed,ndat_max,nsite), &
       z_csed_glob(ncnsv_sed,ndat_max,nsite), &
       etime_csed_glob(ncnsv_sed,ndat_max,nsite), &
       a_csed_glob(ncnsv_sed,ndat_max,nsite), &
       dtime_csed_glob(ncnsv_sed,ndat_max,nsite) )
    allocate (ndat_glob(ncnsv,nsite))
    allocate (ndat_sed_glob(ncnsv_sed,nsite))
    allocate (a_norm_sed_glob(ncnsv_sed,nsite))
    allocate(time_cdat_glob(ncnsv, ndat_max, nsite), &
                    z_cdat_glob(ncnsv, ndat_max, nsite), &
                    a_cdat_glob(ncnsv, ndat_max, nsite) )

    ! load in data files for scalars in cost function

    siteloop: do isite=1,nsite

       do ixc=1,ncnsv
          open(unit=iun, &
               file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
               trim(dat_prefix_glob(isite))//'_'// &
               trim(fname_cost_suffix(ixc))//'.dat', &
               status='old')
          read(iun,*) ndat_glob(ixc,isite)
          ! verify ndat large enough
          call check_ndat(ndat_glob(ixc,isite))
          do id=1,ndat_glob(ixc,isite)
             read(iun,*)time_loc(id,ixc),z_loc(id,ixc),a_loc(id,ixc)
          enddo
          close(iun)
       enddo

       ! load in data file for sediment flux
       do ixc=1,ncnsv_sed
          open(unit=iun, &
               file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
               trim(dat_prefix_glob(isite))//'_'// &
               trim(fname_sed_suffix(ixc))//'.dat', &
               status='old')
          read(iun,*)ndat_sed_glob(ixc,isite)
          call check_ndat(ndat_sed_glob(ixc,isite))
          do id=1,ndat_sed_glob(ixc,isite)
             read(iun,*)btime_csed_glob(ixc,id,isite),z_csed_glob(ixc,id,isite), &
                  dtime_csed_glob(ixc,id,isite),a_csed_glob(ixc,id,isite)
                  etime_csed_glob(ixc,id,isite) = btime_csed_glob(ixc,id,isite)+ &
                  dtime_csed_glob(ixc,id,isite)
           end do
          close(iun)
       end do

       ! read in cost function weights (inverse)
       open(unit=iun, &
            file=trim(datdir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'//trim(fname_var_suffix), &
            status='old')
       ! read in the inverse weights
       do ixc=1,ncnsv
          read(iun,*)w_cost_glob(ixc,isite)
       end do
       do ixc=1,ncnsv_sed
          read(iun,*)w_cost_sed_glob(ixc,isite)
       end do
       close(iun)

       ! bin data into model grid, averaging as necessary
       varloop: do ixc=1,ncnsv
          ncount=0
          dataloop: do idat=1,ndat_glob(ixc,isite)
          
             iz = 1
             do while (z_loc(idat,ixc) > zifc_glob(iz+1,isite) .and. &
                  iz < nz_glob(isite))
                iz = iz+1
             end do
             ! check to see if cost data is below grid, disregard if so.
             if (z_cdat(idat,ixc) > zifc_glob(nz_glob(isite)+1,isite)) then
                cycle dataloop
             end if

             istep=int((time_loc(idat,ixc)-cdays_glob(1,isite))/ &
                  delt_glob(isite)*SecPerDay)+1
             time_cdat_glob(ixc, idat, isite) = istep
             z_cdat_glob(ixc, idat, isite) = iz  
             a_cdat_glob(ixc, idat, isite) =  &
                   a_cdat_glob(ixc,idat,isite)+a_loc(idat,ixc) 
             ncount(ixc,idat)=ncount(ixc,idat)+1

          end do dataloop

          where (ncount .gt. 0)
             a_cdat_glob(:,:,isite)=a_cdat_glob(:,:,isite)/ncount
          endwhere
       end do varloop

       ! compute normalizations
       a_norm_glob(:,isite) = c1/((ncnsv+ncnsv_sed)*ndat_glob(:,isite))
       a_norm_sed_glob(:,isite) = c1/((ncnsv+ncnsv_sed)*ndat_sed_glob(:,isite))
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
    integer, intent(in) :: ndat_test

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
    integer :: isite
    double precision, dimension(:,:) :: cpred
    double precision, dimension(:,:) :: csed_pred

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    ! counters
    integer ixc,it,iz,ised,idat

    ! unit number
    integer, parameter :: iun=91

    ! output cost variable model equivalents
    do ixc=1,ncnsv
       open(unit=iun, &
            file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'// &
            trim(fname_cost_suffix(ixc))//'.out')
       write(iun,'(i10)') ndat(ixc)
       do idat=1,ndat(ixc)
                write(iun,'(f8.3,f8.1,f15.9)') cdays(time_cdat(ixc,idat)), &
                     zmid(z_cdat(ixc,idat)),cpred(ixc,idat)
       enddo
       close(iun)
    enddo
    
    ! write out model equivalents of sediment trap flux
    do ixc=1,ncnsv_sed
        open(unit=iun, &
             file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
             trim(dat_prefix_glob(isite))//'_'//trim(fname_sed_suffix(ixc))//'.out')
        write(iun,'(i10)') ndat_sed(ixc)
        do ised=1,ndat_sed(ixc)
           write(iun,'(f8.3,f8.3,f10.3,f15.9)') &
                btime_csed(ixc,ised),z_csed(ixc,ised), &
                dtime_csed(ixc,ised),csed_pred(ixc,ised)
        enddo
        close(iun)
    end do

  end subroutine output


subroutine calccost(fcost,cpred,csed_pred)
  use const, only : c0
  use common_mod, only : ioflag,dat_prefix_loc
  use grid, only: nz,nt
  implicit none

  !----------------------------------------------------------------------  
  ! Arguments
  !----------------------------------------------------------------------

  ! fcost         cost function value
  ! cpred         predicted values of cost variables
  ! csed_pred     predicted values of sediment flux
  double precision, intent(out) :: fcost

  double precision, dimension(:,:) :: cpred
  double precision, dimension(:,:) :: csed_pred

  !----------------------------------------------------------------------  
  ! Local variables
  !----------------------------------------------------------------------
    
  ! f_cost        cost function (each quantity)
  ! f_cost_sed    cost function for sediment trap fluxes
  double precision, dimension(ncnsv) :: f_cost
  double precision, dimension(ncnsv_sed) :: f_cost_sed

  ! counters/indices
  integer :: ixc,ised,idat
  integer :: id_isnan
  !----------------------------------------------------------------------
  ! Calculate cost function for scalars
  !----------------------------------------------------------------------
  f_cost=c0
  do ixc=1,ncnsv
     if (w_cost(ixc).GT.c0) then
              do idat = 1,ndat(ixc)
              f_cost(ixc) = f_cost(ixc) + &
              a_norm(ixc) * (  w_cost(ixc)* (a_cdat(ixc,idat)-cpred(ixc,idat)) )**2
             enddo 
    endif
  enddo
  !----------------------------------------------------------------------
  ! Compute the contribution from the sediment trap results
  !----------------------------------------------------------------------

!  f_cost_sed = a_norm_sed*sum((w_cost_sed* &
!       (a_csed(1:ndat_sed)-csed_pred(1:ndat_sed)))**2)

  f_cost_sed = c0

  if(ncnsv_sed.gt.0) then
     do ixc = 1, ncnsv_sed 
          do ised = 1,ndat_sed(ixc)
             f_cost_sed(ixc) = f_cost_sed(ixc) + &
                   a_norm_sed(ixc)*( w_cost_sed(ixc)* &
                  (a_csed(ixc,ised)-csed_pred(ixc,ised)) )**2
         end do
     end do
  endif


  !----------------------------------------------------------------------
  ! Normalize the cost values based on number of data points
  !----------------------------------------------------------------------

  fcost=(sum(f_cost) + sum(f_cost_sed))/ndat_tot
  ! debugging only 
  !fcost=f_cost_sed/ndat_tot


    write(6,*)trim(dat_prefix_loc)//' '//'Normalized cost = ', &
       sngl(fcost)
    write(6,*)trim(dat_prefix_loc)//' '//'NO3 cost = ', &
         sngl(f_cost(xcnit)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'PO4 cost =', &
         sngl(f_cost(xcpo4)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'Z cost = ', &
         sngl(f_cost(xczoo)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'PHY cost = ', &
         sngl(f_cost(xcphy)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'CHL cost =', &
         sngl(f_cost(xcchl)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'PP cost =', &
         sngl(f_cost(xcppr)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'BAC cost =', &
         sngl(f_cost(xcbac)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'BPr cost =', &
         sngl(f_cost(xcbpr)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'sDOC cost =', &
         sngl(f_cost(xcdoc)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'sDON cost =', &
         sngl(f_cost(xcdon)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'sDOP cost =', &
         sngl(f_cost(xcdop)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'POC cost =', &
         sngl(f_cost(xcpoc)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'PON cost =', &
         sngl(f_cost(xcpon)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'POP cost =', &
         sngl(f_cost(xcpop)/ndat_tot) !Luo_MB
    write(6,*)trim(dat_prefix_loc)//' '//'STc cost =', &
         sngl(f_cost_sed(xcstc)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'STn cost =', &
         sngl(f_cost_sed(xcstn)/ndat_tot)
    write(6,*)trim(dat_prefix_loc)//' '//'STp cost =', &
         sngl(f_cost_sed(xcstp)/ndat_tot)
    if (ioflag) then
       write(cost_un,'(12(1x,1PG13.6))') fcost,f_cost/ndat_tot, &
            f_cost_sed/ndat_tot
    end if

  end subroutine calccost



  subroutine setup_cost(isite)
    implicit none

    ! isite           location index
    integer :: isite
    
    ndat = ndat_glob(:, isite)
    ndat_sed = ndat_sed_glob(:,isite)
    ndat_tot = ndat_tot_glob(isite)
    a_norm = a_norm_glob(:,isite)
    a_norm_sed = a_norm_sed_glob(:,isite)
    a_csed = a_csed_glob(:,:,isite)
    w_cost = w_cost_glob(:,isite)
    w_cost_sed = w_cost_sed_glob(:,isite)
    btime_csed = btime_csed_glob(:,:,isite)
    dtime_csed = dtime_csed_glob(:,:,isite)
    etime_csed = etime_csed_glob(:,:,isite)
    z_csed = z_csed_glob(:,:,isite)
    time_cdat = time_cdat_glob(:,:,isite)
    z_cdat = z_cdat_glob(:,:,isite)
    a_cdat = a_cdat_glob(:,:,isite)

  end subroutine setup_cost


end module cost
