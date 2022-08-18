!----------------------------------------------------------------------------
!     CVS:$Id: adlight.F90,v 1.6 2005/01/26 19:10:22 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module adlight
  !-------------------------------------------------------------------------
  ! Code adapted from routines generated by 
  ! Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2
  ! Modified by Jeff Dusenberry
  !-------------------------------------------------------------------------
  implicit none

contains

  subroutine adcalc_light(q_surf,chl,adchl,adpar,adpar_ifc)
    use const, only : c0,c2
    use grid, only : dzt,nz,zmid
    use light, only : kchl,kw
    implicit none

    !==============================================
    ! define arguments
    !==============================================
    double precision, dimension(:) :: chl,adchl,adpar,adpar_ifc
    double precision :: q_surf

    !==============================================
    ! define local variables
    !==============================================
    double precision adkb
    integer iz
    double precision kb
    ! TODO: verify nz not causing memory problems here
    double precision, dimension(nz+1) :: par_ifc

    !----------------------------------------------
    ! RESET LOCAL ADJOINT VARIABLES
    !----------------------------------------------
    adkb = c0
    
    !q_surf2=max(0.d0,-12.3214d0+1.2347d0*q_surf+0.0016d0*q_surf*q_surf)
    par_ifc(1) = q_surf
    do iz = 1, nz-1
       kb = kw+kchl*chl(iz)
       par_ifc(iz+1) = par_ifc(iz)*exp(-(kb*dzt(iz)))
    end do
    do iz = nz, 1, -1
       kb = kw+kchl*chl(iz)
       adkb = adkb-adpar_ifc(iz+1)*par_ifc(iz)*dzt(iz)*exp(-(kb*dzt(iz)))
       adpar_ifc(iz) = adpar_ifc(iz)+adpar_ifc(iz+1)*exp(-(kb*dzt(iz)))
       adpar_ifc(iz+1) = c0
       adkb = adkb-adpar(iz)*par_ifc(iz)*(dzt(iz)/c2)*exp(-(kb*(dzt(iz)/c2)))
       adpar_ifc(iz) = adpar_ifc(iz)+adpar(iz)*exp(-(kb*(dzt(iz)/c2)))
       adpar(iz) = c0
       adchl(iz) = adchl(iz)+adkb*kchl
       adkb = c0
    end do

  end subroutine adcalc_light




end module adlight