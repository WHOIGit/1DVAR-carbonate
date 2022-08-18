!----------------------------------------------------------------------------
!     CVS:$Id: model_mod.F90,v 1.18 2005/01/31 21:59:43 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module model_mod
  implicit none

contains



  !-------------------------------------------------------------------------
  ! map from biological parameter list to optimizable parameter list
  !-------------------------------------------------------------------------
  subroutine bio_to_opt(bioparams,optparams)

    use common_mod, only : bioparams_default
    use const, only : c0
    use eco_params, only : bio_opt_map,nparams_opt,nparams_bio

    implicit none

    double precision, dimension(:) :: optparams, bioparams

    integer :: ipar

    do ipar=1,nparams_opt
       if (bioparams_default(bio_opt_map(ipar)) /= c0) then
          optparams(ipar)=log(bioparams(bio_opt_map(ipar))/ &
               bioparams_default(bio_opt_map(ipar)))
       else
          optparams(ipar)=log(bioparams(bio_opt_map(ipar)))
       endif
    end do
       
  end subroutine bio_to_opt

  !-------------------------------------------------------------------------
  ! map from optimizable parameter list to biological parameter list
  !-------------------------------------------------------------------------
  subroutine opt_to_bio(optparams,bioparams)

    use common_mod, only : bioparams_default,bioparams0
    use const, only : c0
    use eco_params, only : bio_opt_map,nparams_opt,nparams_bio
    implicit none

    double precision, dimension(:) :: optparams, bioparams

    integer :: ipar

    bioparams=bioparams0
    do ipar=1,nparams_opt
       if (bioparams_default(bio_opt_map(ipar)) /= c0) then
          bioparams(bio_opt_map(ipar)) = &
               exp(optparams(ipar))*bioparams_default(bio_opt_map(ipar))
       else
          bioparams(bio_opt_map(ipar)) = exp(optparams(ipar))
       endif
    end do

  end subroutine opt_to_bio


! A few comments on TAMC - some intrinsics, like SUM, LOG, and SQRT can fail
! to generate results with TAMC.  Most of the time, this seems only to be 
! a problem when the calls are performed from within subroutines contained
! within a module, or when a parent subroutine is contained within a 
! module.  The TAMC preprocessor directives were set up to handle this
! situation - with some success.

  subroutine model(optparams,fcost)

    use common_mod, only : bio0,nsite,nz_max,nt_max, iotemp_un, iotempflag
    use common_mod, only : setup_common
    use const, only : c0
    use cost, only : ncnsv,ncnsv_sed,ndat_max
    use cost, only : calccost,setup_cost
    use eco_derivs, only : get_costpred
    use eco_derivs, only : setup_bio
    use eco_params, only : NumStateVar,NumDiagVar,nparams_bio
    use forcing, only : setup_forcing
    use grid, only : nz,delt,nt,ntsout
    use grid, only : setup_grid
    use cost, only : output
    use io, only : write_eco_out


    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! optparams           optimizable biological parameters (various units)
    ! bio0                initial conditions on state variables
    double precision, dimension(:) :: optparams

    ! fcost               value of cost function
    double precision :: fcost
    

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    ! TODO: make allocatable?
    double precision, dimension(nz_max+1,NumStateVar,nt_max) :: bio
    double precision, dimension(nz_max,NumDiagVar,nt_max) :: diag

    ! counters
    integer :: istep,isite,iout

    double precision, dimension(ncnsv, ndat_max) :: cpred
    double precision, dimension(ncnsv_sed,ndat_max) :: csed_pred

    double precision, dimension(nz_max+1,NumStateVar) :: bio_prev,bio_next,bio_temp
    double precision, dimension(nz_max,NumDiagVar) :: diag_prev,diag_next,diag_temp

    double precision, dimension(nparams_bio) :: bioparams

    double precision, dimension(nsite) :: fcost_loc
    
!    integer :: ispin
    
    fcost_loc = c0

    ! set bioparams from optimizable set
    call opt_to_bio(optparams,bioparams)

    siteloop: do isite=1,nsite

       call setup_common(isite)
       call setup_grid(isite)
       call setup_forcing(isite)
       call setup_bio(isite)
       call setup_cost(isite)
                     
       !allocate(bio(nz+1,NumStateVar,nt))
       !allocate(diag(nz,NumDiagVar,nt))
       ! initialize biological array from initial conditions
       bio(1:nz+1,:,1)=bio0
       ! initialize diagnostics to zero.  
       diag(1:nz,:,1)=c0

!       do ispin=1,3
!ADJ LOOP = SEQUENTIAL
       do istep=2,nt
       
          bio_prev(1:nz+1,:) = bio(1:nz+1,:,istep-1)
          diag_prev(1:nz,:) = diag(1:nz,:,istep-1)


! TODO: find better way than 1:nz stuff ??
          !tamc call ts_init(bio_sm,diag_sm,istep)
          call ts_integrate(istep,bio_prev(1:nz+1,:),bio_next(1:nz+1,:), &
               diag_prev(1:nz,:),diag_next(1:nz,:), &
               bioparams)
          call ts_final(istep,bio_next(1:nz+1,:),diag_next(1:nz,:))
          bio(1:nz+1,:,istep) = bio_next(1:nz+1,:)
          diag(1:nz,:,istep) = diag_next(1:nz,:)/delt
       enddo

       call get_costpred(cpred,csed_pred,bio(1:nz+1,:,:), &
            diag(1:nz,:,:),bioparams)
       call calccost(fcost_loc(isite),cpred,csed_pred)

       call output(isite,cpred,csed_pred)
       call write_eco_out(isite,bio(1:nz+1,:,:),diag(1:nz,:,:))
                  
       if (iotempflag) then
           close(iotemp_un)
       end if

    end do siteloop

    fcost = sum(fcost_loc)
    write(6,*) 'Total Cost = ',fcost


  end subroutine model

  subroutine CALFUN(nparams_opt,optparams,fcost)

     use common_mod, only : bio0,nsite,nz_max,nt_max, iotemp_un, iotempflag,ioflag
    use common_mod, only : setup_common
    use const, only : c0
    use cost, only : ncnsv,ncnsv_sed,ndat_max
    use cost, only : calccost,setup_cost
    use eco_derivs, only : get_costpred
    use eco_derivs, only : setup_bio
    use eco_params, only : NumStateVar,NumDiagVar,nparams_bio
    use forcing, only : setup_forcing
    use grid, only : nz,delt,nt,ntsout
    use grid, only : setup_grid
    use cost, only : output
    use io, only : write_eco_out, param_un


    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! optparams           optimizable biological parameters (various units)
    ! bio0                initial conditions on state variables
    double precision, dimension(:) :: optparams

    ! fcost               value of cost function
    double precision :: fcost
    
     integer :: nparams_opt

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    ! TODO: make allocatable?
    double precision, dimension(:,:,:), allocatable :: bio, diag

    ! counters
    integer :: istep,isite,iout

    double precision, dimension(:,:), allocatable :: cpred
    double precision, dimension(ncnsv_sed,ndat_max) :: csed_pred

    double precision, dimension(nz_max+1,NumStateVar) :: bio_prev,bio_next,bio_temp
    double precision, dimension(nz_max,NumDiagVar) :: diag_prev,diag_next,diag_temp

    double precision, dimension(nparams_bio) :: bioparams

    double precision, dimension(nsite) :: fcost_loc
    
!    integer :: ispin
    
    fcost_loc = c0

    ! set bioparams from optimizable set
    call opt_to_bio(optparams,bioparams)
    
    if (ioflag) then   !Luo
          write(param_un,'(100(1x,1PG16.8E3,:))') bioparams
    end if

    siteloop: do isite=1,nsite

       call setup_common(isite)
       call setup_grid(isite)
       call setup_forcing(isite)
       call setup_bio(isite)
       call setup_cost(isite)
      
       allocate(bio(nz+1,NumStateVar,nt/ntsout))
       allocate(diag(nz,NumDiagVar,nt/ntsout))
       ! initialize biological array from initial conditions
       bio_prev(1:nz+1,:)=bio0
       ! initialize diagnostics to zero.  
       diag_prev(1:nz,:)=c0
       iout = 0


!       do ispin=1,3
!ADJ LOOP = SEQUENTIAL
       do istep=2,nt

! TODO: find better way than 1:nz stuff ??
          !tamc call ts_init(bio_sm,diag_sm,istep)
          call ts_integrate(istep,bio_prev(1:nz+1,:),bio_next(1:nz+1,:), &
               diag_prev(1:nz,:),diag_next(1:nz,:), &
               bioparams)
          call ts_final(istep,bio_next(1:nz+1,:),diag_next(1:nz,:))
          bio_prev = bio_next
          diag_prev = diag_next/delt
          bio_temp = bio_temp + bio_prev
          diag_temp = diag_temp + diag_prev
          if (mod(istep, ntsout).eq.0) then
                  iout = iout + 1
                  bio(1:nz+1, :, iout) = bio_temp(1:nz+1,:)/ntsout
                  diag(1:nz, :, iout) = diag_temp(1:nz,:)/ntsout
                  bio_temp = c0
                  diag_temp = c0
          end if
       enddo

       call get_costpred(cpred(1:nz,:),csed_pred,bio(1:nz+1,:,:), &
            diag(1:nz,:,:),bioparams)
       call calccost(fcost_loc(isite),cpred(1:nz,:),csed_pred)
        
       if (iotempflag) then
           close(iotemp_un)
       end if

    end do siteloop

    fcost = sum(fcost_loc)
    write(6,*) 'Total Cost = ',fcost

  end subroutine CALFUN
  
  !-------------------------------------------------------------------------
  ! calculate derivatives for each timestep
  !-------------------------------------------------------------------------
  subroutine ts_derivs(istep,bio_prev,dydt_bio,diag_prev,dydt_diag,bioparams)
    use common_mod, only : nz_max
    use const, only : c0
    use eco_params, only : numstatevar
    use eco_derivs, only : bioderivs
    use grid, only : nz
    use physderivs_mod, only : physderivs_init,physderivs
    use eco_derivs, only : bioderivs_init
    
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    integer istep

    double precision, dimension(:,:), intent(in) :: bio_prev
    double precision, dimension(:,:), intent(out) :: dydt_bio,dydt_diag
    double precision, dimension(:), intent(in) :: bioparams
    double precision, dimension(:,:), intent(in) :: diag_prev

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    
    ! par 
    double precision, dimension(nz) :: par
    double precision, dimension(nz+1) :: par_ifc
    double precision, dimension(nz+1,numstatevar) :: dydt_phys
    

    dydt_bio=c0
    dydt_diag=c0
    dydt_phys=c0
    call bioderivs_init(istep,par,par_ifc,bio_prev,diag_prev)
    call bioderivs(istep,par,par_ifc,bio_prev,dydt_bio,dydt_diag,bioparams)
    call physderivs_init(istep,bio_prev)
    call physderivs(istep,bio_prev,dydt_phys,bioparams)
    dydt_bio = dydt_bio + dydt_phys

  end subroutine ts_derivs


  !-------------------------------------------------------------------------
  ! Step forward one time step
  !-------------------------------------------------------------------------
  subroutine ts_integrate(istep,bio_prev,bio_next,diag_prev,diag_next, &
       bioparams)
    use common_mod, only : nz_max
    use const, only : p5
    use eco_params, only : NumStateVar,NumDiagVar
    use grid, only : nz,delt
   

    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! istep             current time step
    integer istep

    double precision, dimension(:,:),intent(in) :: bio_prev
    double precision, dimension(:,:),intent(out) :: bio_next,diag_next
    double precision, dimension(:),intent(in) :: bioparams
    double precision, dimension(:,:), intent(in) :: diag_prev

    !-----------------------------------------------------------------------
    ! Local variables 
    !-----------------------------------------------------------------------
    double precision, dimension(nz+1,NumStateVar) :: dydt_bio, &
         bio_temp,dydt_bio_temp
    double precision, dimension(nz,NumDiagVar) :: dydt_diag, &
         dydt_diag_temp,diag_temp


    ! TODO: remove nz references everywhere with new dz framework
    !       they shouldn't be needed

    ! do integration
!    if(istep.eq.2) then
       ! Euler step 
       call ts_derivs(istep,bio_prev,dydt_bio,diag_prev,dydt_diag,bioparams)
       bio_next(1:nz,:)=bio_prev(1:nz,:)+ &
            delt*dydt_bio(1:nz,:) 
       diag_next=delt*dydt_diag
!    else
!       call ts_derivs(istep,bio_prev,dydt_bio,diag_prev,dydt_diag,bioparams)
!       bio_temp(1:nz,:)=bio_prev(1:nz,:)+ &
!            delt*dydt_bio(1:nz,:)
!       diag_temp = dydt_diag
!       call ts_derivs(istep,bio_temp,dydt_bio_temp,diag_temp, &
!            dydt_diag_temp,bioparams)
!       bio_next(1:nz,:)=bio_prev(1:nz,:)+ &
!            delt*p5*(dydt_bio_temp(1:nz,:)+dydt_bio(1:nz,:))
!       diag_next=delt*p5*(dydt_diag_temp+dydt_diag)       
!    endif

  end subroutine ts_integrate


  !-------------------------------------------------------------------------
  ! Perform any steps necessary at end of every integration time step.
  !-------------------------------------------------------------------------
  subroutine ts_final(istep,bio_next,diag_next)
    use const, only : c0
    use eco_params, only : NumStateVar
    use forcing, only : dml
    use grid, only : dzt,zifc,nz
    
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------
    integer :: istep
    double precision, dimension(:,:) :: bio_next,diag_next

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    integer :: ndml,isv,iz
    double precision ds

    ! redistribute mixed layer state variables
    ! calculate no. of levels in mixed layer
!    ndml=int(dml(istep)/dz)

    ndml = 0
    iz = 1
    do while (dml(istep) > zifc(iz+1) .and. iz <= nz)
       ndml = iz
       iz = iz + 1
    end do
       
    if (ndml .gt. 0) then
       do isv = 1,NumStateVar
          ! TAMC - the sum() call seems to create TAMC failures
          ! ds=sum(bio_next(1:ndml,isv),1)/ndml
          ds=c0
          do iz=1,ndml
             ds=ds+bio_next(iz,isv)*dzt(iz)
          enddo
          ds = ds/zifc(ndml+1)
          bio_next(1:ndml,isv)=ds
       enddo
    endif

  end subroutine ts_final



end module model_mod
