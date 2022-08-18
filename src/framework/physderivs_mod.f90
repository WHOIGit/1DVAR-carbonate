!----------------------------------------------------------------------------
!     CVS:$Id: physderivs_mod.F90,v 1.16 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module physderivs_mod
  use const, only : c0
  
  use common_mod, only : iotempflag, iotemp_un
  
  implicit none

  ! weight on horizontal advection
  double precision, parameter :: hadv_wt = 1.0

  ! smalln parameter for flux limited advection/sinking scheme
  double precision , parameter             :: smalln=2.23d-16


contains

  subroutine physderivs(istep,bio_prev,dydt,bioparams)
    use const, only : c0
    use eco_params, only : numstatevar
    use grid, only : nz,dzt,ntsout
    
    implicit none
    
    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev
    double precision, dimension(:,:) :: bio_prev
    double precision, dimension(:,:), intent(out) :: dydt
    double precision, dimension(:), intent(in) :: bioparams

    ! istep               current time step
    integer :: istep, ii,eu
    double precision, dimension(nz+1,numstatevar) :: dydt_tmp
    double precision, dimension(numstatevar) :: dydt_avg
    dydt_tmp = c0
    call vertmix(bio_prev,dydt_tmp,istep)
    dydt = dydt_tmp
    dydt_tmp = c0
    call horizadv(bio_prev,dydt_tmp,istep)
    dydt = dydt + dydt_tmp
    dydt_tmp = c0
    call vertadv(istep,bio_prev,dydt_tmp)
    dydt = dydt + dydt_tmp
    dydt_tmp = c0
    call set_sflux(bio_prev,dydt_tmp,istep)
    dydt = dydt + dydt_tmp
    dydt_tmp = c0
    call sink(istep,bio_prev,dydt_tmp,bioparams)
    dydt = dydt + dydt_tmp
    dydt_tmp = c0
    
    eu = 18
    if ((mod(istep,ntsout).eq.0).and.(iotempflag)) then
        do ii = 1,numstatevar
            dydt_avg(ii) = sum(dydt(1:nz, ii)*dzt(1:nz))/sum(dzt(1:nz))
        end do
	write(iotemp_un, '(50(1x, 1PG14.6))') dydt_avg
	do ii = 1,numstatevar
            dydt_avg(ii) = sum(dydt(1:eu, ii)*dzt(1:eu))/sum(dzt(1:eu))
        end do
        write(iotemp_un, '(50(1x, 1PG14.6))') dydt_avg
    end if
    !~ write(iotemp_un, *) istep, sum(dydt(1:nz,2)*dzt(1:nz)), &
                                          !~ sum(dydt(1:nz,5)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,8)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,11)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,14)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,17)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,20)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,23)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,26)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,28)*dzt(1:nz)), &
					  !~ sum(dydt(1:nz,29)*dzt(1:nz))

  end subroutine physderivs


  subroutine vertmix(bio_prev,dydt,istep)
    use eco_params, only : NumStateVar
    use grid, only : nz,delt,zmid,dzt
    use const, only : c0,c1,c2,p5
    use forcing, only : rkz

    use numeric_subs, only : tridag
    implicit none

    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev
    double precision, dimension(:,:) :: bio_prev,dydt

    ! istep             current time step
    integer :: istep

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    ! tridiagonal solver coefficients
    double precision, dimension(nz) :: tria,trib,tric,rr,bio_new
    
    ! counters
    integer :: iz,isv

    dydt = c0

    iz=1
    tria(iz) = c0
    tric(iz) = -p5*delt*rkz(iz+1,istep)/(zmid(iz+1)-zmid(iz))/dzt(iz)
    trib(iz) = c1 - tric(iz)
    do iz=2,nz
       tria(iz) = -p5*delt*rkz(iz,istep)/(zmid(iz)-zmid(iz-1))/dzt(iz)
       tric(iz) = -p5*delt*rkz(iz+1,istep)/(zmid(iz+1)-zmid(iz))/dzt(iz)
       trib(iz) = c1-tria(iz)-tric(iz)
    end do

    trib(nz) = c1-tria(nz)

    do isv=1,NumStateVar
       !     set up right hand side for tridag

       iz=1
       rr(iz) = bio_prev(iz,isv) + &
            (p5*delt/(zmid(iz+1)-zmid(iz))/dzt(iz))*rkz(iz+1,istep) * &
            (bio_prev(iz+1,isv)-bio_prev(iz,isv))
       do iz=2,nz-1
          rr(iz) = bio_prev(iz,isv) + &
               (p5*delt/(zmid(iz)-zmid(iz-1))/dzt(iz))*rkz(iz,istep)* &
               (bio_prev(iz-1,isv)-bio_prev(iz,isv)) + &
               (p5*delt/(zmid(iz+1)-zmid(iz))/dzt(iz))*rkz(iz+1,istep)* &
               (bio_prev(iz+1,isv)-bio_prev(iz,isv))
       end do

       rr(nz) = bio_prev(nz,isv) + &
               (p5*delt/(zmid(nz)-zmid(nz-1))/dzt(nz))*rkz(nz,istep)* &
               (bio_prev(nz-1,isv)-bio_prev(nz,isv)) + &
               (delt/(zmid(nz+1)-zmid(nz))/dzt(nz))*rkz(nz+1,istep)* &
               (bio_prev(nz+1,isv)-bio_prev(nz,isv))

       call tridag(tria,trib,tric,rr,bio_new,nz)

       do iz=1,nz
          dydt(iz,isv)=(bio_new(iz)-bio_prev(iz,isv))/delt
       enddo
    enddo

  end subroutine vertmix


  subroutine horizadv(bio_prev,dydt,istep)
    use common_mod, only : lhoriz_adv,horiz_adv
    use const, only : c0
    use eco_params, only : NumStateVar
    use grid, only : nz
    implicit none


    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev
    double precision, dimension(:,:) :: bio_prev,dydt

    ! istep             current time step
    integer :: istep

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    ! counters, indices
    integer :: iz

    dydt = c0

    !#ifdef TAMC
!    do iz=1,nz
    do iz=16,nz    ! restrict to below euphotic zone

       ! TODO:
    ! NOTE:  removed restriction to below euphotic zone to accomodate
    ! differences in depth resolution.  Restriction should be moved to
    ! the input file (set variable horiz_adv = 0). - JAD

       where(lhoriz_adv)
          dydt(iz,:) = -hadv_wt*horiz_adv(iz,:,istep)* &
                bio_prev(iz,:)
       end where
       
    end do

  end subroutine horizadv


  subroutine set_sflux(bio_prev,dydt,istep)

    use common_mod, only : aeoflux
    use const, only : c0,HourPerDay
    use eco_common, only : aeonsv
    use grid, only : dzt
    implicit none


    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev
    double precision, dimension(:,:) :: bio_prev,dydt
    integer :: istep

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------


    dydt = c0
    where (aeonsv .ne. 0)
       dydt(1,:)=aeoflux(:,istep)/(dzt(1)*HourPerDay)
    endwhere

  end subroutine set_sflux


  subroutine sink(istep,bio_prev,dydt,bioparams)
    use const, only : c0,c1,c2,rc6,SecPerDay
    use eco_common, only : wnsvflag
    use eco_params, only : NumStateVar
    use grid, only : nz,rdzv,rdzt,delt
    use eco_params, only : iwnsvo
    implicit none


    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev
    double precision, dimension(:), intent(in) :: bioparams
    double precision, dimension(:,:) :: bio_prev,dydt
    integer :: istep

    !-----------------------------------------------------------------------
    ! Subroutine adapted from MITgcm routine.  See vertadv subroutine
    ! in this module for details.
    !-----------------------------------------------------------------------
    
    ! loop counters and field indices
    integer :: k, km2, km1, kp1, isv

    ! abbreviations and flux limiters
    double precision    :: Rjp, Rj, Rjm, wFld, wP, wM, cfl
    double precision    :: psiPRj, psiMRj, d0, d1
    double precision, dimension(nz+1) :: wFlux
    double precision :: wFluxkp1 ! auxillary variable


    ! following added by JAD to get adjoint 
    double precision :: psi_work
    double precision, dimension(NumStateVar) :: wnsv
!---------------------------------------------------------------------
    wnsv = bioparams(iwnsvo) * wnsvflag
    
    dydt=c0

    do isv=1, numstatevar

       ! some initializations
       wFlux          = c0
       wFluxkp1       = c0
       wFld           = c0

       do k=nz+1,2,-1

          ! use the "right" sign: downward velocity wFld is negative, 
          ! if sinking (i.e., downward) velocity ws is positive 
          wFld = -wnsv(isv)/SecPerDay

          ! some abbreviations
          wP   = wFld+abs(wFld)
          wM   = wFld-abs(wFld)

          kp1 = min(k+1,nz+1)   ! use nz+1 to include bbc
          km1 = max(k-1,1)
          km2 = max(k-2,1)

          Rjp=bio_prev(km2,isv)-bio_prev(km1,isv)
          Rj =bio_prev(km1,isv)-bio_prev(k,isv)
          Rjm=bio_prev(k,isv)-bio_prev(kp1,isv)


          ! compute Courant number cfl
          cfl=abs(wFld*delt*rdzv(k))
          ! DST3 parameters
          d0=(c2-cfl)*(c1-cfl)*rc6
          d1=(c1-cfl*cfl)*rc6

          psiPRj=d0*Rj+d1*Rjm
          psiMRj = d0*Rj+d1*Rjp

          if (Rj > c0) then

             if (Rj < psiPRj) then
                psiPRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjm
             if (psi_work < psiPRj) then
                psiPRj = psi_work
             end if
             if (c0 > psiPRj) then
                psiPRj = c0
             end if

             if (Rj < psiMRj) then
                psiMRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjp
             if (psi_work < psiMRj) then
                psiMRj = psi_work
             end if
             if (c0 > psiMRj) then
                psiMRj = c0
             end if

          else

             if (Rj > psiPRj) then
                psiPRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjm
             if (psi_work > psiPRj) then
                psiPRj = psi_work
             end if
             if (c0 < psiPRj) then
                psiPRj = c0
             end if

             if (Rj > psiMRj) then
                psiMRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjp
             if (psi_work > psiMRj) then
                psiMRj = psi_work
             end if
             if (c0 < psiMRj) then
                psiMRj = c0
             end if

          end if

          wflux(k)=                                      &
               ( 0.5*wP*(bio_prev(k,isv)   + psiPRj )   &
               + 0.5*wM*(bio_prev(km1,isv) - psiMRj ) )

          ! sink due to sinking for layer/cell k
          ! minus sign because this has been moved to the right hand side
          if (k .le. nz) then
             dydt(k,isv) = (wFluxkp1 - wflux(k))*rdzt(k)
          end if
          ! store flux at level kp1 for the next cycle
          wfluxkp1 = wflux(k)
       end do
       k=1
       wflux(k)= c0
       dydt(k,isv) = (wFluxkp1-wflux(k))*rdzt(k)
    end do

  end subroutine sink



  subroutine vertadv(istep,bio_prev,dydt)
    use const, only : c0,c1,c2,rc6,SecPerDay
    use eco_params, only : NumStateVar
    use forcing, only : wvel
    use grid, only : nz,rdzv,rdzt,delt

    implicit none


    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    ! dydt                time derivatives for bio_prev

    double precision, dimension(:,:) :: bio_prev,dydt

    integer :: istep


!---------------------------------------------------------------------
!      subroutine recom_sinking_dst3fl(dt,k_max,dz,ws,c,sink)
!---------------------------------------------------------------------    
!     code adopted from the MITgcm routine: 
!     /==========================================================\
!     | SUBROUTINE GAD_DST3_ADV_R                                |
!     | o Compute Vertical advective Flux of Tracer using        |
!     |   3rd Order DST Sceheme with flux limiting               |
!     |==========================================================|
!
!     This routine assumes that all vertical points are "wet" and that there
!     are no boundaries other than the surface and the absolute bottom 
!     (k_max) of the the domain. If this is not the case then one has to
!     include a masking field for the tracer field c. For simplicity, this 
!     is omitted here.
!
!     Author: Martin Losch (Alfred Wegener Institute, AWI), July, 2003
!     Modifications: Markus Schartau 
!     (Marine Sciences Research Center, Stony Brook; January, 2005)
!     Modifications: Jeff Dusenberry, February, 2005
!     Modifications: John Dunne, February, 2005

    ! loop counters and field indices
    integer :: k, km2, km1, kp1, isv

    ! abbreviations and flux limiters
    double precision    :: Rjp, Rj, Rjm, wFld, wP, wM, cfl
    double precision    :: psiPRj, psiMRj, d0, d1
    double precision , dimension(nz+1) :: wFlux
    double precision :: wFluxkp1 ! auxillary variable


    ! following added by JAD to get adjoint 
    double precision :: psi_work

!---------------------------------------------------------------------
    dydt=c0

    do isv=1, numstatevar

       ! some initializations
       wFlux          = c0
       wFluxkp1       = c0
       wFld           = c0

       do k=nz+1,1,-1

          ! use the "right" sign: downward velocity wFld is negative, 
          ! if sinking (i.e., downward) velocity ws is positive 
          wFld = -wvel(k,istep)

          ! some abbreviations
          wP   = wFld+abs(wFld)
          wM   = wFld-abs(wFld)

          kp1 = min(k+1,nz+1)   ! use nz+1 to include bbc
          km1 = max(k-1,1)
          km2 = max(k-2,1)

          Rjp=bio_prev(km2,isv)-bio_prev(km1,isv)
          Rj =bio_prev(km1,isv)-bio_prev(k,isv)
          Rjm=bio_prev(k,isv)-bio_prev(kp1,isv)


          ! compute Courant number cfl
          cfl=abs(wFld*delt*rdzv(k))
          ! DST3 parameters
          d0=(c2-cfl)*(c1-cfl)*rc6
          d1=(c1-cfl*cfl)*rc6

          ! compute flux limiters psiP, and psiM
          ! compute product psiP*Rj ("psiPRj") and psiM*Rj ("psiMRj")
          ! instead of thetaP and psiP to avoid division by zero issue
          !jad          thetaP=Rjm/(1.D-20+Rj)
          !jad          psiP=d0+d1*thetaP
          psiPRj=d0*Rj+d1*Rjm
          psiMRj = d0*Rj+d1*Rjp

          ! don't use min/max functions as in original code - TAMC seems 
          ! to not quite get them correct.
          if (Rj > c0) then

             if (Rj < psiPRj) then
                psiPRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjm
             if (psi_work < psiPRj) then
                psiPRj = psi_work
             end if
             if (c0 > psiPRj) then
                psiPRj = c0
             end if

             if (Rj < psiMRj) then
                psiMRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjp
             if (psi_work < psiMRj) then
                psiMRj = psi_work
             end if
             if (c0 > psiMRj) then
                psiMRj = c0
             end if

          else

             if (Rj > psiPRj) then
                psiPRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjm
             if (psi_work > psiPRj) then
                psiPRj = psi_work
             end if
             if (c0 < psiPRj) then
                psiPRj = c0
             end if

             if (Rj > psiMRj) then
                psiMRj = Rj
             end if
             psi_work = (c1-cfl)/(smalln+cfl)*Rjp
             if (psi_work > psiMRj) then
                psiMRj = psi_work
             end if
             if (c0 < psiMRj) then
                psiMRj = c0
             end if

          end if

          wflux(k)=                                     &
               ( 0.5*wP*(bio_prev(k,isv)   + psiPRj )   &
               + 0.5*wM*(bio_prev(km1,isv) - psiMRj ) )

          ! sink due to sinking for layer/cell k
          ! minus sign because this has been moved to the right hand side
          if (k .le. nz) then
             dydt(k,isv) = (wFluxkp1 - wflux(k) + bio_prev(k,isv) &
                             * (wvel(kp1,istep) - wvel(k,istep)))*rdzt(k)
          end if
           
          ! store flux at level kp1 for the next cycle
          wfluxkp1 = wflux(k)
       end do
    end do
  end subroutine vertadv



  subroutine physderivs_init(istep,bio_prev)
    use common_mod, only : bbcnsv, flag_bbcnsv
    use const, only : c0,c1,c2,mc1
    use grid, only : nz

    implicit none


    !-----------------------------------------------------------------------
    ! Arguments
    !-----------------------------------------------------------------------

    ! bio_prev(nz,nsv)    state variable profiles (various units)
    double precision, dimension(:,:) :: bio_prev
    integer :: istep

    where(bbcnsv(:,istep) .lt. mc1)
       bio_prev(nz+1,:)=bio_prev(nz,:)- &
            (bio_prev(nz-1,:)-bio_prev(nz,:))
       ! constrain to be non-negative
       where (bio_prev(nz+1,:) < c0)
          bio_prev(nz+1,:) = c0
       end where
       !bio_prev(NumDepth+1,:) = max(bio_prev(NumDepth+1,:),c0)
    elsewhere
        where(flag_bbcnsv(:) .eq. c1)
               bio_prev(nz+1,:)=bbcnsv(:,istep-1)
        elsewhere(flag_bbcnsv(:) .ge.c2)
               bio_prev(nz+1,:)=bio_prev(nz,:)+bbcnsv(:,istep-1)
        endwhere
       where (bio_prev(nz+1,:) < c0)
          bio_prev(nz+1,:) = c0
       end where   
    endwhere
 
!    bio_prev(nz+1,:)=bbcnsv(:,istep-1)

  end subroutine physderivs_init





end module physderivs_mod
