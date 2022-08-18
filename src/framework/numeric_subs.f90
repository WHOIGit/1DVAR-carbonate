!----------------------------------------------------------------------------
!     CVS:$Id: numeric_subs.F90,v 1.3 2004/07/30 19:17:41 duse Exp $
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module numeric_subs
  implicit none

contains

!---------------------------------------------------------------------+
!     SUBROUTINE TRIDAG: Solves for a vector U of length N the      |
!     tridiagonal linear set given by the Crank-Nicholson scheme.    |
!     A, B, C and R are input vectors and are not modified.  This   |
!     subroutine is from Press et al., Numerical Recipes, p. 40.
!
! Adapted for use with TAMC - Jeff Dusenberry 04 May 2004
!---------------------------------------------------------------------+
  subroutine tridag(a,b,c,r,u,n)

    implicit none

    integer n

    double precision, dimension(n) :: gam,xu,bL
    double precision, dimension(:) :: a,b,c,r,u

    ! TAMC - change bet from scalar to array to avoid recomputation in
    ! adjoint
    double precision, dimension(n) :: bet

    integer :: j

    ! TAMC - turn off error checking logic for computation of adjoint
    if(b(1).eq.0.d0) then
       write(*,*) 'first diagonal entry is zero'
       write(*,*) ' stop in tridag.F'
       stop
    endif

    bet(1)=b(1)
    xu(1)=r(1)/bet(1)

    do j=2,n
       gam(j) = c(j-1)/bet(j-1)
       bet(j)=b(j)-a(j)*gam(j)
       if(bet(j).eq.0.d0) then
          write(*,*) j,'th diagonal entry is zero'
          write(*,*) ' stop in tridag.F'
          stop         
       endif
       xu(j)=(r(j)-a(j)*xu(j-1))/bet(j)
    enddo

    bL(n)=xu(n)
    do j=n-1,1,-1
       bL(j)=xu(j)-gam(j+1)*xu(j+1)
       xu(j)=bL(j)
    enddo

    do j=1,n
       u(j)=bL(j)
    enddo

  end subroutine tridag


end module numeric_subs
