!----------------------------------------------------------------------------
!     CVS:$Id: const.F90,v 1.7 2005/01/31 21:59:43 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module const

  implicit none

  double precision, parameter :: c0=0.0,c1=1.0, &
       c1p5=1.5,c2=2.0,c2p5=2.5, &
       c3=3.0, &
       c4=4.0,c6=6.0,c10=10.0,p5=0.5
  double precision, parameter :: mc1=-c1,rc6=c1/c6

  double precision, parameter :: &
       SecPerDay=86400.0, &
       HourPerDay=24.0, &
       SecPerHour=3600.0, &
       dpy=365.0, &
       dps=c1/SecPerDay

end module const
