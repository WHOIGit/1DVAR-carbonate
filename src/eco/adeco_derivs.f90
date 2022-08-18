!----------------------------------------------------------------------------
!     CVS:$Id: adeco_derivs.F90,v 1.9 2005/01/31 22:09:32 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module adeco_derivs
!-------------------------------------------------------------------------
! Code adapted from routines generated by 
! Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2
! Modified by Jeff Dusenberry
!-------------------------------------------------------------------------
  implicit none


contains

! ***Required***
  subroutine adbioderivs(istep,par,par_ifc,bio_prev,bioparams,adpar, &
       adpar_ifc,adbio_prev,addydt_bio,addydt_diag,adbioparams)
    use eco_common, only : rlt
    use grid, only : nz
    use adeco_common, only : rltb
    use adderivs_mod, only : adderivs
    implicit none

!==============================================
! define arguments
!==============================================
    double precision, dimension(:,:) :: adbio_prev,addydt_bio, &
         addydt_diag,bio_prev
    double precision, dimension(:) :: adbioparams,adpar,bioparams,par, &
         adpar_ifc,par_ifc
    integer :: istep
!==============================================
! define local variables
!==============================================
    integer :: iz

!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
    do iz = nz, 1, -1
       rlt = par(iz)
       call adderivs(bio_prev(iz,:),bioparams,adbio_prev(iz,:), &
            addydt_bio(iz,:),addydt_diag(iz,:),istep,iz,adbioparams)
       adpar(iz) = adpar(iz)+rltb
       rltb = 0.
    end do


  end subroutine adbioderivs


! ***Required***
  subroutine adget_costpred(bio,diag,bioparams,adcpred,adcsed_pred,adbio, &
       addiag,adbioparams)
    use cost, only : ncnsv,xcppr,xcphy,xcchl,xczoo,xcnit,xcpo4,xcbac,xcbpr, &
                            xcdoc,xcdon,xcdop,xcpoc,xcpon,xcpop,xcstc,xcstn,xcstp !Luo_MB
    use cost, only : btime_csed,etime_csed,ndat_max,ndat_sed,z_csed
    use cost, only: ndat, time_cdat, z_cdat, a_cdat
    use const, only : c0,c1,c2,secperhour,secperday
    use eco_common, only : wnsvflag
    use eco_params, only : iDETc,iNH4,iNO3,ipp,iremin,iremin_prf_n,&
                  iremin_prf_p,iMZc,mwC,mwN,mwP,numstatevar
    use eco_params, only: iPO4,iBAc,iprBAc,iSDOMc,iSDOMn,iSDOMp,&
                                     iSPc,iTRc,iUNc,iPRTc,iDETc, &
				     iSPn,iSPp,iTRn,iTRp,iUNn,iUNp,iPRTn,iPRTp, &
				     iDETn,iDETp,iSPchl,iTRchl,iUNchl !Luo_MB
    use forcing, only : cdays
    use grid, only : nz,nt,zifc,delt,ntsout
    use eco_params, only: iwnsvo
    implicit none

!==============================================
! define arguments
!==============================================
    double precision, dimension(:,:) ::  adcsed_pred,adcpred
    double precision, dimension(:,:,:) :: bio,diag,adbio,addiag
    double precision, dimension(:) :: bioparams,adbioparams

!==============================================
! define local variables
!==============================================
    integer ised,ised_loc,idat
    integer it
    integer it2, it_temp
!double precision nsed_pred(ndat_max)
    integer, dimension(ndat_max) :: nsed_pred
    double precision t2trap
    double precision z2trap
    integer :: &
         ntsperday      , & ! number of timesteps per day
         ntsperhalfday      ! number of timesteps per half day

   double precision, dimension(numstatevar) :: wnsv,adwnsv



!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
wnsv = bioparams(iwnsvo) * wnsvflag
!nsed_pred = c0
    nsed_pred = 0
    do it = 1, nt
       do ised = 1, ndat_sed(xcstc)
          ised_loc = nz+1
          if  (z_csed(xcstc,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstc,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstc,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETc)
          if (t2trap .gt. btime_csed(xcstc,ised) .and. &
               t2trap .le. etime_csed(xcstc,ised)) then
!nsed_pred(ised) = nsed_pred(ised)+c1
             nsed_pred(ised) = nsed_pred(ised)+1
          endif
       end do
    end do
!where (nsed_pred(1:ndat_sed) .gt. c0)
    where (nsed_pred(1:ndat_sed(xcstc)) .gt. 0)
       adcsed_pred(xcstc,1:ndat_sed(xcstc)) = &
                adcsed_pred(xcstc,1:ndat_sed(xcstc))/nsed_pred(1:ndat_sed(xcstc))
    endwhere
    do it = nt, 1, -1
       it_temp = it
       do ised = 1, ndat_sed(xcstc)
          ised_loc = nz+1
          if  (z_csed(xcstc,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstc,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstc,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETc)
          if (t2trap .gt. btime_csed(xcstc,ised) .and. &
               t2trap .le. etime_csed(xcstc,ised)) then                   
             adbio(ised_loc-1,iDETc,it_temp) = adbio(ised_loc-1,iDETc,it_temp)+ &
                  adcsed_pred(xcstc,ised)*wnsv(iDETc)*mwC* &
                  exp(-(z2trap*bioparams(iremin)/wnsv(iDETc)))
             adbioparams(iremin) = adbioparams(iremin)- &
                  adcsed_pred(xcstc,ised)*wnsv(iDETc)*mwC* &
                  bio(ised_loc-1,iDETc,it_temp)*(z2trap/wnsv(iDETc))* &
                  exp(-(z2trap*bioparams(iremin)/wnsv(iDETc)))
             adwnsv(iDETc) = adwnsv(iDETc) + adcsed_pred(xcstc,ised)  * &
                  mwC*bio(ised_loc-1,iDETc,it_temp)* &
                  ( exp(-z2trap*bioparams(iremin)/wnsv(iDETc)) + &
                  z2trap*bioparams(iremin)/wnsv(iDETc) *&
                  exp(-z2trap*bioparams(iremin)/wnsv(iDETc))  )
          endif
       end do
    end do
    
    nsed_pred = 0
    do it = 1, nt
       do ised = 1, ndat_sed(xcstn)
          ised_loc = nz+1
          if  (z_csed(xcstn,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstn,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstn,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETn)
          if (t2trap .gt. btime_csed(xcstn,ised) .and. &
               t2trap .le. etime_csed(xcstn,ised)) then
!nsed_pred(ised) = nsed_pred(ised)+c1
             nsed_pred(ised) = nsed_pred(ised)+1
          endif
       end do
    end do
!where (nsed_pred(1:ndat_sed) .gt. c0)
    where (nsed_pred(1:ndat_sed(xcstn)) .gt. 0)
       adcsed_pred(xcstn,1:ndat_sed(xcstn)) = &
                adcsed_pred(xcstn,1:ndat_sed(xcstn))/nsed_pred(1:ndat_sed(xcstn))
    endwhere
    do it = nt, 1, -1
       it_temp = it
       do ised = 1, ndat_sed(xcstn)
          ised_loc = nz+1
          if  (z_csed(xcstn,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstn,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstn,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETn)
          if (t2trap .gt. btime_csed(xcstn,ised) .and. &
               t2trap .le. etime_csed(xcstn,ised)) then               
             adbio(ised_loc-1,iDETn,it_temp) = adbio(ised_loc-1,iDETn,it_temp)+ &
                  adcsed_pred(xcstn,ised)*wnsv(iDETn)*mwN* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn)))
             adbioparams(iremin) = adbioparams(iremin)- &
                  adcsed_pred(xcstn,ised)*wnsv(iDETn)*mwN* &
                  bio(ised_loc-1,iDETn,it_temp)*(z2trap*bioparams(iremin_prf_n)/wnsv(iDETn))* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn)))
             adbioparams(iremin_prf_n) = adbioparams(iremin_prf_n)- &
                  adcsed_pred(xcstn,ised)*wnsv(iDETn)*mwN* &
                  bio(ised_loc-1,iDETn,it_temp)*(z2trap*bioparams(iremin)/wnsv(iDETn))* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn)))
             adwnsv(iDETn) = adwnsv(iDETn) + adcsed_pred(xcstn,ised) * &
                  mwN*bio(ised_loc-1,iDETn,it_temp) * &
                  (  exp(-z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn))  + &
                  z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn) * &
                  exp(-z2trap*bioparams(iremin)*bioparams(iremin_prf_n)/wnsv(iDETn)) )
          endif
       end do
    end do
    
    nsed_pred = 0
    do it = 1, nt
       do ised = 1, ndat_sed(xcstp)
          ised_loc = nz+1
          if  (z_csed(xcstp,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstp,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstp,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETp)
          if (t2trap .gt. btime_csed(xcstp,ised) .and. &
               t2trap .le. etime_csed(xcstp,ised)) then
!nsed_pred(ised) = nsed_pred(ised)+c1
             nsed_pred(ised) = nsed_pred(ised)+1
          endif
       end do
    end do
!where (nsed_pred(1:ndat_sed) .gt. c0)
    where (nsed_pred(1:ndat_sed(xcstp)) .gt. 0)
       adcsed_pred(xcstp,1:ndat_sed(xcstp)) = &
                adcsed_pred(xcstp,1:ndat_sed(xcstp))/nsed_pred(1:ndat_sed(xcstp))
    endwhere
    do it = nt, 1, -1
       it_temp = it
       do ised = 1, ndat_sed(xcstp)
          ised_loc = nz+1
          if  (z_csed(xcstp,ised).lt.zifc(nz+1)) then
             do while (z_csed(xcstp,ised).lt.zifc(ised_loc))
                      ised_loc = ised_loc - 1
             end do
          endif  
          z2trap = z_csed(xcstp,ised)-zifc(ised_loc)
          t2trap = cdays(it)+z2trap/wnsv(iDETp)
          if (t2trap .gt. btime_csed(xcstp,ised) .and. &
               t2trap .le. etime_csed(xcstp,ised)) then               
             adbio(ised_loc-1,iDETp,it_temp) = adbio(ised_loc-1,iDETp,it_temp)+ &
                  adcsed_pred(xcstp,ised)*wnsv(iDETp)*mwP* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp)))
             adbioparams(iremin) = adbioparams(iremin)- &
                  adcsed_pred(xcstp,ised)*wnsv(iDETp)*mwP* &
                  bio(ised_loc-1,iDETp,it_temp)*(z2trap*bioparams(iremin_prf_p)/wnsv(iDETp))* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp)))
             adbioparams(iremin_prf_p) = adbioparams(iremin_prf_p)- &
                  adcsed_pred(xcstp,ised)*wnsv(iDETp)*mwP* &
                  bio(ised_loc-1,iDETp,it_temp)*(z2trap*bioparams(iremin)/wnsv(iDETp))* &
                  exp(-(z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp)))
             adwnsv(iDETp) = adwnsv(iDETp) + adcsed_pred(xcstp,ised) * &
                   mwP*bio(ised_loc-1,iDETp,it_temp) * &
                  (  exp(-z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp))  + &
                  z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp) * &
                  exp(-z2trap*bioparams(iremin)*bioparams(iremin_prf_p)/wnsv(iDETp)) )
          endif
       end do
    end do
    
    adcsed_pred = c0
    

    ntsperday = nint(secperday/delt)
    ntsperhalfday = nint(secperday/delt/c2)
    
    adcpred(xcppr,:) = adcpred(xcppr,:) * secperday / ntsperday
    do idat = 1,ndat(xcppr)
        do it = time_cdat(xcppr,idat)-ntsperhalfday+ntsperday-1, &
	           time_cdat(xcppr,idat)-ntsperhalfday, -1
                   it_temp = it
		   addiag(z_cdat(xcppr, idat),iPP,it_temp) = &
		       addiag(z_cdat(xcppr, idat),iPP,it_temp) + adcpred(xcppr, idat)
        enddo 
    enddo
    
    adcpred(xcbpr,:) = adcpred(xcbpr,:) * secperday / ntsperday !Luo_MB
    do idat = 1,ndat(xcbpr)
        do it = time_cdat(xcbpr,idat)-ntsperhalfday+ntsperday-1, &
	           time_cdat(xcbpr,idat)-ntsperhalfday, -1
                   it_temp = it
		   addiag(z_cdat(xcbpr, idat),iprBAc,it_temp) = &
		       addiag(z_cdat(xcbpr, idat),iprBAc,it_temp) + adcpred(xcbpr, idat)
        enddo 
    enddo
    
    do idat = 1, ndat(xczoo)
            adbio(z_cdat(xczoo,idat), iMZc, time_cdat(xczoo,idat)) = &
               adbio(z_cdat(xczoo,idat), iMZc, time_cdat(xczoo,idat)) + adcpred(xczoo, idat)
            adcpred(xczoo, idat)  = c0   
    enddo 
    
    do idat = 1, ndat(xcnit)
           adbio(z_cdat(xcnit,idat), iNO3, time_cdat(xcnit,idat)) = &
               adbio(z_cdat(xcnit,idat), iNO3, time_cdat(xcnit,idat)) + adcpred(xcnit, idat)  
           adcpred(xcnit, idat) = c0 
    enddo 
        
    do idat = 1, ndat(xcphy)
         adbio(z_cdat(xcphy,idat), iSPn, time_cdat(xcphy,idat)) = &
               adbio(z_cdat(xcphy,idat), iSPn, time_cdat(xcphy,idat)) + adcpred(xcphy, idat)
		 adbio(z_cdat(xcphy,idat), iUNn, time_cdat(xcphy,idat)) = &
               adbio(z_cdat(xcphy,idat), iUNn, time_cdat(xcphy,idat)) + adcpred(xcphy, idat)
         adcpred(xcphy, idat) = c0
    enddo
    
    do idat = 1, ndat(xcchl)
         adbio(z_cdat(xcchl,idat), iSPchl, time_cdat(xcchl,idat)) = &
               adbio(z_cdat(xcchl,idat), iSPchl, time_cdat(xcchl,idat)) + adcpred(xcchl, idat)
         adbio(z_cdat(xcchl,idat), iTRchl, time_cdat(xcchl,idat)) = &
               adbio(z_cdat(xcchl,idat), iTRchl, time_cdat(xcchl,idat)) + adcpred(xcchl, idat) 
         adbio(z_cdat(xcchl,idat), iUNchl, time_cdat(xcchl,idat)) = &
               adbio(z_cdat(xcchl,idat), iUNchl, time_cdat(xcchl,idat)) + adcpred(xcchl, idat)
         adcpred(xcchl, idat) = c0
    enddo
    
    do idat = 1, ndat(xcpo4)
          adbio(z_cdat(xcpo4,idat), iPO4, time_cdat(xcpo4,idat)) = &
               adbio(z_cdat(xcpo4,idat), iPO4, time_cdat(xcpo4,idat)) + adcpred(xcpo4, idat)
	     adcpred(xcpo4, idat) = c0
    enddo

    do idat = 1, ndat(xcbac)
         adbio(z_cdat(xcbac,idat), iBAc, time_cdat(xcbac,idat)) = &
               adbio(z_cdat(xcbac,idat), iBAc, time_cdat(xcbac,idat)) + adcpred(xcbac, idat)
	     adcpred(xcbac, idat) = c0
    enddo 
    
    do idat = 1, ndat(xcdoc)
         adbio(z_cdat(xcdoc,idat), iSDOMc, time_cdat(xcdoc,idat)) = &
               adbio(z_cdat(xcdoc,idat), iSDOMc, time_cdat(xcdoc,idat)) + adcpred(xcdoc, idat)
	     adcpred(xcdoc, idat) = c0
    enddo 
    
    do idat = 1, ndat(xcdon)
         adbio(z_cdat(xcdon,idat), iSDOMn, time_cdat(xcdon,idat)) = &
               adbio(z_cdat(xcdon,idat), iSDOMn, time_cdat(xcdon,idat)) + adcpred(xcdon, idat)
	     adcpred(xcdon, idat) = c0
    enddo
    
    do idat = 1, ndat(xcdop)
         adbio(z_cdat(xcdop,idat), iSDOMp, time_cdat(xcdop,idat)) = &
               adbio(z_cdat(xcdop,idat), iSDOMp, time_cdat(xcdop,idat)) + adcpred(xcdop, idat)
	     adcpred(xcdop, idat) = c0
    enddo
    
    do idat = 1, ndat(xcpoc)
         adbio(z_cdat(xcpoc,idat), iSPc, time_cdat(xcpoc,idat)) = &
               adbio(z_cdat(xcpoc,idat), iSPc, time_cdat(xcpoc,idat)) + adcpred(xcpoc, idat)
         adbio(z_cdat(xcpoc,idat), iUNc, time_cdat(xcpoc,idat)) = &
               adbio(z_cdat(xcpoc,idat), iUNc, time_cdat(xcpoc,idat)) + adcpred(xcpoc, idat)
         adbio(z_cdat(xcpoc,idat), iTRc, time_cdat(xcpoc,idat)) = &
               adbio(z_cdat(xcpoc,idat), iTRc, time_cdat(xcpoc,idat)) + adcpred(xcpoc, idat)
         adbio(z_cdat(xcpoc,idat), iPRTc, time_cdat(xcpoc,idat)) = &
               adbio(z_cdat(xcpoc,idat), iPRTc, time_cdat(xcpoc,idat)) + adcpred(xcpoc, idat)
         adbio(z_cdat(xcpoc,idat), iDETc, time_cdat(xcpoc,idat)) = &
               adbio(z_cdat(xcpoc,idat), iDETc, time_cdat(xcpoc,idat)) + adcpred(xcpoc, idat)    
	     adcpred(xcpoc, idat) = c0
    enddo
    
    do idat = 1, ndat(xcpon)
         adbio(z_cdat(xcpon,idat), iSPn, time_cdat(xcpon,idat)) = &
               adbio(z_cdat(xcpon,idat), iSPn, time_cdat(xcpon,idat)) + adcpred(xcpon, idat)
         adbio(z_cdat(xcpon,idat), iUNn, time_cdat(xcpon,idat)) = &
               adbio(z_cdat(xcpon,idat), iUNn, time_cdat(xcpon,idat)) + adcpred(xcpon, idat)
         adbio(z_cdat(xcpon,idat), iTRn, time_cdat(xcpon,idat)) = &
               adbio(z_cdat(xcpon,idat), iTRn, time_cdat(xcpon,idat)) + adcpred(xcpon, idat)
         adbio(z_cdat(xcpon,idat), iPRTn, time_cdat(xcpon,idat)) = &
               adbio(z_cdat(xcpon,idat), iPRTn, time_cdat(xcpon,idat)) + adcpred(xcpon, idat)
         adbio(z_cdat(xcpon,idat), iDETn, time_cdat(xcpon,idat)) = &
               adbio(z_cdat(xcpon,idat), iDETn, time_cdat(xcpon,idat)) + adcpred(xcpon, idat)   
	     adcpred(xcpon, idat) = c0
    enddo
    
    do idat = 1, ndat(xcpop)
         adbio(z_cdat(xcpop,idat), iSPp, time_cdat(xcpop,idat)) = &
               adbio(z_cdat(xcpop,idat), iSPp, time_cdat(xcpop,idat)) + adcpred(xcpop, idat)
         adbio(z_cdat(xcpop,idat), iUNp, time_cdat(xcpop,idat)) = &
               adbio(z_cdat(xcpop,idat), iUNp, time_cdat(xcpop,idat)) + adcpred(xcpop, idat)
         adbio(z_cdat(xcpop,idat), iTRp, time_cdat(xcpop,idat)) = &
               adbio(z_cdat(xcpop,idat), iTRp, time_cdat(xcpop,idat)) + adcpred(xcpop, idat)
         adbio(z_cdat(xcpop,idat), iPRTp, time_cdat(xcpop,idat)) = &
               adbio(z_cdat(xcpop,idat), iPRTp, time_cdat(xcpop,idat)) + adcpred(xcpop, idat)
         adbio(z_cdat(xcpop,idat), iDETp, time_cdat(xcpop,idat)) = &
               adbio(z_cdat(xcpop,idat), iDETp, time_cdat(xcpop,idat)) + adcpred(xcpop, idat)   
	     adcpred(xcpop, idat) = c0
    enddo
    
    adbioparams(iwnsvo) = adbioparams(iwnsvo)+sum(adwnsv*wnsvflag)
    adwnsv = 0.d0
    
    adcpred = c0
    adcsed_pred = c0

  end subroutine adget_costpred



!  subroutine adbioderivs_init(bio_prev,istep,adbio_prev,adpar )
  subroutine adbioderivs_init(istep,bio_prev,diag_prev,adpar,adpar_ifc, &
       adbio_prev,addiag_prev)
    use adlight, only : adcalc_light
    use const, only : c0,p5
    use eco_params, only : iSPchl, iTRchl, iUNchl
    use forcing, only : qi
    use grid, only : nz
    implicit none

!==============================================
! define arguments
!==============================================
    double precision, dimension(:,:) :: bio_prev,adbio_prev
    double precision, dimension(:,:) :: diag_prev,addiag_prev
    double precision, dimension(:) :: adpar,adpar_ifc
    integer :: istep

!==============================================
! define local variables
!==============================================
    double precision, dimension(nz) :: adchl,chl, adchl_temp, chl_temp

!----------------------------------------------
! RESET LOCAL ADJOINT VARIABLES
!----------------------------------------------
    adchl = c0
    
!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
     chl = bio_prev(:,iSPchl) + bio_prev(:,iTRchl) + bio_prev(:,iUNchl)
 !   chl = max(chl_temp,c0)
    call adcalc_light(qi(istep),chl,adchl,adpar,adpar_ifc)
 !   adchl_temp = adchl_temp + adchl * (p5+sign(p5, chl_temp-c0))
 !   adchl = c0
    adbio_prev(:,iSPchl) = adbio_prev(:,iSPchl) + adchl
    adbio_prev(:,iTRchl) = adbio_prev(:,iTRchl) + adchl
    adbio_prev(:,iUNchl) = adbio_prev(:,iUNchl) + adchl

    adchl = c0

  end subroutine adbioderivs_init

end module adeco_derivs