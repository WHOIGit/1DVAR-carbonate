!----------------------------------------------------------------------------
!     CVS:$Id: kinds_mod.F90,v 1.2 2004/07/30 19:17:41 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module kinds_mod

  implicit none

  integer, parameter :: &     
       int_kind  = kind(1), &
       log_kind  = kind(.true.), &
       real_kind = selected_real_kind(6), &
       dbl_kind  = selected_real_kind(13)
       
end module kinds_mod
