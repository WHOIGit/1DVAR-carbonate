module hessian_mod

  use eco_params, only : nparams_opt
  implicit none

  double precision, dimension(nparams_opt,nparams_opt), save :: &
       hessian_inv,eye

contains

  subroutine hessian_init

    use const, only : c0,c1
    implicit none

    
    ! Local variables    
    ! counter
    integer i

    ! initialize to identity matrix
    eye = c0
    do i=1,nparams_opt
       eye(i,i) = c1
    end do

    hessian_inv = eye

  end subroutine hessian_init


  subroutine update_hessian(s,y)
    use const, only : c1,c2
    use eco_params, only : nparams_opt
    implicit none

    ! Arguments
    ! s               step
    ! y               change in gradient
    double precision, dimension(:) :: s, y
    
    
    ! Local variables
    !double precision :: rho
    !double precision, dimension(nparams_opt) :: sstar
    double precision, dimension(nparams_opt,nparams_opt) :: psy,pys,pss
    double precision :: sumsy,sumyy,sumss,eps

    integer :: i,j

    ! it looks like sbar, ybar in n1qn3 are already normalized to rho ???
    !rho = c1 / dot_product(s,y)
    !sstar = rho * s
    
    ! TODO: move to const
    eps = epsilon(c1)

    sumsy = sum(s*y)
    sumyy = sum(y*y)
    sumss = sum(s*s)

    ! don't update hessian if sumsy too small
    if (sumsy**2 > eps*sumyy*sumss) then

       do i=1,nparams_opt
          do j=1,nparams_opt
             !psy(i,j) = sstar(i)*y(j)
             !pss(i,j) = sstar(i)*s(j)
             psy(i,j) = s(i)*y(j)
             pss(i,j) = s(i)*s(j)
          end do
       end do

       pys = transpose(psy)

       ! kind of a brute force approach here...
       ! update hessian, according to BFGS
       hessian_inv = matmul(matmul((eye-psy),hessian_inv),(eye-pys)) + pss
       ! alternative - change order
       !hessian_inv = hessian_inv - matmul(psy,hessian_inv) - &
       !     matmul(hessian_inv,pys) + matmul(matmul(psy,hessian_inv),pys)

!       write(6,*) 'update_hessian: hessian_inv = ',hessian_inv
   
!    else

!       write(6,*) 'update_hessian:  hessian_inv not updated'

    end if
    


  end subroutine update_hessian


end module hessian_mod
