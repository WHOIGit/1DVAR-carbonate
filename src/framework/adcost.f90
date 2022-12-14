!----------------------------------------------------------------------------
!     CVS:$Id: adcost.F90,v 1.10 2005/04/28 17:34:20 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module adcost
  !-------------------------------------------------------------------------
  ! Code adapted from routines generated by 
  ! Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2
  ! Modified by Jeff Dusenberry
  !-------------------------------------------------------------------------
  implicit none

contains

  
  subroutine adcalccost( cpred, csed_pred, adfcost, adcpred, adcsed_pred )
    use const, only : c0,c1
    use cost, only : cost_un
    use cost, only : ncnsv,ncnsv_sed,ndat_max,ndat
    use cost, only : a_cdat,a_norm,ndat_tot,w_cost
    use cost, only : a_csed,a_norm_sed,ndat_sed,w_cost_sed
    use cost, only : xcphy,xcchl,xcnit,xcppr,xczoo
    use grid, only: nz,nt
!    use cost, only: f_cost, f_cost_sed
    implicit none

    !==============================================
    ! define arguments
    !==============================================
    double precision, dimension(:,:) :: adcpred, cpred
    double precision, dimension(:,:) ::  adcsed_pred, csed_pred
    double precision :: adfcost

    !==============================================
    ! define local variables
    !==============================================
    double precision, dimension(ncnsv) :: adf_cost
    double precision, dimension(ncnsv_sed) :: adf_cost_sed
    integer ixc,ised,iz,it,idat
    integer :: id_isnan

    !----------------------------------------------
    ! RESET LOCAL ADJOINT VARIABLES
    !----------------------------------------------
    adf_cost = c0
    adf_cost_sed = c0

    !----------------------------------------------
    ! ROUTINE BODY
    !----------------------------------------------
    adf_cost_sed = adf_cost_sed+adfcost/ndat_tot
    adf_cost = adf_cost+adfcost/ndat_tot
    adfcost = c0
    do ixc = 1, ncnsv_sed
        do ised = 1, ndat_sed(ixc)
           adcsed_pred(ixc,ised) = adcsed_pred(ixc,ised)-2*adf_cost_sed(ixc)*a_norm_sed(ixc)* &
                w_cost_sed(ixc)*w_cost_sed(ixc)*(a_csed(ixc,ised)-csed_pred(ixc,ised))
        end do
    end do
    
    do ixc=1,ncnsv
     if (w_cost(ixc).GT.c0) then
              do idat = 1,ndat(ixc)
	          adcpred(ixc,idat) = adcpred(ixc,idat) - 2*adf_cost(ixc) * a_norm(ixc) * &
	                w_cost(ixc) * w_cost(ixc) * (a_cdat(ixc,idat)-cpred(ixc,idat))
             enddo 
    endif
    
  enddo

  end subroutine adcalccost


end module adcost
