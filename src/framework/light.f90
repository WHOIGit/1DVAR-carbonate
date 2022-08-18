!----------------------------------------------------------------------------
!     CVS:$Id: light.F90,v 1.10 2005/01/19 19:57:34 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------

module light
  
  implicit none
  save

  ! kw                   attenuation coefficient for water (1/m)
  ! kchl                 attenuation coefficient for chl (1/(mg Chl/m3)/m)
  double precision, parameter :: kw=0.038d0,kchl=0.05d0


contains




subroutine calc_light(q_surf,chl,par,par_ifc)
  use const, only : c2
  use grid, only : zmid,nz,dzt
  implicit none

  ! q_surf           surface PAR
  double precision q_surf

  ! chl              chlorophyll profile, model grid
  ! par              photosynthetically active radiation at midpoints
  !                    (units same as q_surf)
  ! par_ifc          PAR at level interfaces (units same as q_surf)

  double precision, dimension(:) :: chl,par,par_ifc

  double precision :: kb
  
  integer :: iz

  !q_surf2=max(0.d0,-12.3214d0+1.2347d0*q_surf+0.0016d0*q_surf*q_surf)
  par_ifc(1)=q_surf

!ADJ LOOP = SEQUENTIAL
  do iz=1,nz
     kb=kw+kchl*chl(iz)
     par(iz)=par_ifc(iz)*exp(-kb*(dzt(iz)/c2))
     par_ifc(iz+1)=par_ifc(iz)*exp(-kb*dzt(iz))
  enddo


end subroutine calc_light



end module light
