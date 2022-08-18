!----------------------------------------------------------------------------
!     CVS:$Id: io.F90,v 1.13 2005/02/16 15:27:40 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module io
  use types, only : ncf_float,ncf_fillval_float
  implicit none

  ! save everything
  save

  ! file unit numbers
  integer, parameter :: param_un=41
  integer, parameter :: results_un=33

  ! output file names
  character(len=*), parameter :: &
       fname_cost_out='costs.out', &
       fname_param_out='params.out', &
       fname_results_out='accepted_iter.out', &
       fname_eco_nc='eco_out.nc'

  ! eco_ncid          netcdf file id
  ! eco_tvarid        netcdf variable it for time variable
  integer :: eco_ncid,eco_tvarid

contains

  !-------------------------------------------------------------------------
  ! Initial I/O for driver - read in directory paths, filenames, etc.
  !-------------------------------------------------------------------------
  subroutine driver_init

    use common_mod, only : homedir,datdir,ecopar_fname
    use common_mod, only : adjpar_fname,optpar_fname
    use common_mod, only : nsite,nsite_max,dat_prefix_glob
    implicit none

    integer, parameter :: stdin=5

    namelist /driver_io_nml/ &
         homedir, &
         datdir, &
         ecopar_fname, &
         adjpar_fname, &
         optpar_fname, &
         nsite, &
         dat_prefix_glob
    
    !-----------------------------------------------------------------------
    !   default namelist settings
    !-----------------------------------------------------------------------
    homedir='../output/'
    datdir='../forcing/'
    ecopar_fname='eco_params.in'
    adjpar_fname='adj_params.in'
    optpar_fname='opt_params.in'
    nsite=nsite_max
    dat_prefix_glob = (/'HOT1989', 'HOT1990', 'HOT1991', 'HOT1992', &
                                     'HOT1993', 'HOT1994', 'HOT1995', 'HOT1996', &
                                     'HOT1997', 'HOT1998', 'HOT1999', 'HOT2000', &
                                     'HOT2001', 'HOT2002', 'HOT2003', 'HOT2004', &
                                     'HOT2005', 'AS     ', 'EqPac  '   /)

    read(unit=stdin, nml=driver_io_nml)

    ! verify nsite <= nsite_max
    if (nsite .gt. nsite_max) then
       write(6,*) 'ERROR: nsite_max not large enough.'
       write(6,*) 'change value of nsite_max and recompile.'
       write(6,*) 'nsite_max = ',nsite_max,' but need nsite = ',nsite
       stop
    endif

  end subroutine driver_init


  !-------------------------------------------------------------------------
  ! I/O initialization - open files, etc.
  !-------------------------------------------------------------------------
  subroutine io_init(lcostout)
    use common_mod, only : homedir,ioflag,iotempflag
    use cost, only : cost_un
    use eco_params, only : param_names,get_opt_param_names,nparams_opt
    implicit none

    logical :: lcostout
    
    character(len=20), dimension(nparams_opt) :: opt_param_names
    ioflag = lcostout

    if (ioflag) then
       ! open file for output of cost function values
       open(unit=cost_un,file=trim(homedir)//'/'//trim(fname_cost_out), &
            status='replace')
       
       ! open output file for parameter values
       open(unit=param_un,file=trim(homedir)//'/'//trim(fname_param_out), &
            status='replace')
       write(param_un,'(100(1x,a17,:))') param_names


       ! open output file for results
       open(unit=results_un,file=trim(homedir)//'/'//trim(fname_results_out), &
            status='replace')
       call get_opt_param_names(opt_param_names)
       write(results_un,'(102(1x,a20,:))') 'iter','tot_cost',opt_param_names

    end if
    
  end subroutine io_init

  !-------------------------------------------------------------------------
  ! I/O finalization - close files, etc.
  !-------------------------------------------------------------------------
  subroutine io_final
    use cost, only : cost_un
    use common_mod, only : ioflag,iotempflag
    implicit none

    if (ioflag) then
       close(param_un)
       close(cost_un)
       close(results_un)
    end if

  end subroutine io_final


  !-------------------------------------------------------------------------
  ! write out ecosystem state variables and diagnostics
  !-------------------------------------------------------------------------

subroutine write_eco_out(isite,bio,diag)
    use grid, only : nt,nz,ntsout,zmid,delt
    use eco_common, only : bio_ncvaratts,diag_ncvaratts
    use eco_params, only : NumStateVar,NumDiagVar,fname_bio_suffix
    use forcing, only :cdays
    use common_mod, only : homedir,dat_prefix_glob
    use netcdf
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------
    integer :: isite               ! location index
    double precision, dimension(:,:,:) :: bio,diag

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    integer isv,inz,iun


    double precision, allocatable :: ncvar_temp(:,:)
    
    double precision, allocatable :: cdays_temp(:)
    
    integer ntsout_start, nt_temp, intsout

    call create_ncfiles(isite)
    
    ntsout_start = max(1,int(ntsout/2))
    nt_temp = int(nt/ntsout)
    allocate (ncvar_temp(nz,nt_temp))
    allocate (cdays_temp(nt_temp))
    
    cdays_temp = cdays(ntsout_start:(ntsout_start+(nt_temp-1)*ntsout):ntsout)
    call check_nc(nf90_put_var(ncid=eco_ncid, &
         varid=eco_tvarid,values=cdays_temp,start=(/1/),count=(/nt_temp/)))
    do isv=1,NumStateVar
       ncvar_temp = bio(1:nz,isv,1:nt:ntsout)
       do intsout = 2,ntsout
           ncvar_temp = ncvar_temp + bio(1:nz,isv,intsout:nt:ntsout)
       end do
       ncvar_temp =  ncvar_temp/ntsout
       call check_nc(nf90_put_var(ncid=eco_ncid, &
            varid=bio_ncvaratts(isv)%id, &
            values=ncvar_temp))
    end do

    do isv=1,NumDiagVar
       ncvar_temp = diag(:,isv,1:nt:ntsout)
       do intsout = 2,ntsout
           ncvar_temp = ncvar_temp + diag(1:nz,isv,intsout:nt:ntsout)
       end do
       ncvar_temp =  ncvar_temp/ntsout 
       call check_nc(nf90_put_var(ncid=eco_ncid, &
            varid=diag_ncvaratts(isv)%id, &
            values=ncvar_temp))
    end do

    call check_nc(nf90_close(eco_ncid))
    
    do  isv=1,NumStateVar
            open(unit=iun, &
            file=trim(homedir)//'/'//trim(dat_prefix_glob(isite))//'/'// &
            trim(dat_prefix_glob(isite))//'_'// &
            trim(fname_bio_suffix(isv))//'.init')
           do  inz=1,nz
             !! [HK] 
             !! Print *, inz, nt, bio(inz,isv,floor(nt-(86400/delt)*90))   
             write(iun, '(f7.1,e12.4)') zmid(inz), bio(inz,isv,floor(nt-(86400/delt)*90))
           end do
           close(iun)
    end do

  end subroutine write_eco_out 


  subroutine create_ncfiles(isite)
    use common_mod, only : homedir,dat_prefix_glob
    use eco_common, only : bio_ncvaratts,diag_ncvaratts
    use eco_params, only : NumStateVar,NumDiagVar
    use grid, only : nz,zmid
    use netcdf
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------
    integer :: isite               ! location index
    
    
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    ! [zt]dimid         dimension ids for depth and time
    ! [zt]varid         variable ids for depth and time
    integer :: zdimid,tdimid,zvarid

    ! counters
    integer :: isv


    ! create scalar output file
    call check_nc(nf90_create(path=trim(homedir)//'/'// &
         trim(dat_prefix_glob(isite))//'/'// &
         trim(dat_prefix_glob(isite))//'_'// &
         trim(fname_eco_nc), &
         cmode=nf90_clobber,ncid=eco_ncid))
    ! define dimensions for scalar file
    call create_dimensions(eco_ncid,zdimid,tdimid,zvarid,eco_tvarid)

    ! create variables for scalar output
    do isv=1,NumStateVar
       call check_nc(nf90_def_var(ncid=eco_ncid, &
            name=bio_ncvaratts(isv)%short_name, &
            xtype=bio_ncvaratts(isv)%type,dimids=(/zdimid,tdimid/), &
            varid=bio_ncvaratts(isv)%id))
       call check_nc(nf90_put_att(ncid=eco_ncid,varid=bio_ncvaratts(isv)%id, &
            name="long_name",values=trim(bio_ncvaratts(isv)%long_name)));
       call check_nc(nf90_put_att(ncid=eco_ncid,varid=bio_ncvaratts(isv)%id, &
            name="units",values=trim(bio_ncvaratts(isv)%units)));
       call check_nc(nf90_put_att(ncid=eco_ncid,varid=bio_ncvaratts(isv)%id, &
            name="_FillValue",values=NF90_FILL_FLOAT));
       call check_nc(nf90_put_att(ncid=eco_ncid,varid=bio_ncvaratts(isv)%id, &
            name="missing_value",values=NF90_FILL_FLOAT));
    end do

    ! create variables for diagnostics
    do isv=1,NumDiagVar
       call check_nc(nf90_def_var(ncid=eco_ncid, &
            name=diag_ncvaratts(isv)%short_name, &
            xtype=diag_ncvaratts(isv)%type,dimids=(/zdimid,tdimid/), &
            varid=diag_ncvaratts(isv)%id))
       call check_nc(nf90_put_att(ncid=eco_ncid, &
            varid=diag_ncvaratts(isv)%id, &
            name="long_name",values=trim(diag_ncvaratts(isv)%long_name)));
       call check_nc(nf90_put_att(ncid=eco_ncid, &
            varid=diag_ncvaratts(isv)%id, &
            name="units",values=trim(diag_ncvaratts(isv)%units)));
       call check_nc(nf90_put_att(ncid=eco_ncid, &
            varid=diag_ncvaratts(isv)%id, &
            name="_FillValue",values=NF90_FILL_FLOAT));
       call check_nc(nf90_put_att(ncid=eco_ncid, &
            varid=diag_ncvaratts(isv)%id, &
            name="missing_value",values=NF90_FILL_FLOAT));
    end do


    ! put output file into data mode
    call check_nc(nf90_enddef(eco_ncid))
    
    ! TODO: remove need for nz?
    ! write dimension data to depth
    call check_nc(nf90_put_var(ncid=eco_ncid, &
         varid=zvarid, values=zmid(1:nz),start=(/1/), &
         count=(/nz/)))

  end subroutine create_ncfiles


  subroutine check_nc(status)

    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
       write(6,*) trim(nf90_strerror(status))
    end if

  end subroutine check_nc
    
  subroutine create_dimensions(ncid,zdimid,tdimid,zvarid,tvarid)

    use grid, only : nz,nt
    use netcdf
    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(out) :: zvarid,tvarid
    integer, intent(out) :: zdimid,tdimid

    call check_nc(nf90_def_dim(ncid=ncid, name="depth", &
         len=nz, dimid=zdimid))
    call check_nc(nf90_def_dim(ncid=ncid, name="time", &
         len=NF90_UNLIMITED, dimid=tdimid))

    ! create variables corresponding to above dimensions
    call check_nc(nf90_def_var(ncid=ncid, name="depth", &
         xtype=ncf_float,dimids=(/zdimid/), &
         varid=zvarid))
    call check_nc(nf90_put_att(ncid=ncid,varid=zvarid, &
         name="units",values="meters"));
    call check_nc(nf90_def_var(ncid=ncid, name="time", &
         xtype=ncf_float,dimids=(/tdimid/), &
         varid=tvarid))
    ! TODO: get units correct on time - should be "days since ..."
    call check_nc(nf90_put_att(ncid=ncid,varid=tvarid, &
         name="units",values="days"));

  end subroutine create_dimensions


  subroutine write_matrix(fname,matrix,header)

    implicit none
    
    double precision, dimension(:,:) :: matrix
    character(len=*) :: fname
    character(len=*), dimension(:) :: header

    integer, parameter :: un = 43

    ! msize              size of matrix (should be square)
    ! i                  counter
    integer :: msize, i

    msize = size(matrix(:,1))
    open(unit=un,file=fname,status='replace')
    write(un,'(100(1x,a20,:))') header(:)
    do i=1,msize
       write(un,'(100(1x,1PG20.12,:))') matrix(i,:)
    end do
    close(un)
    
  end subroutine write_matrix




end module io
