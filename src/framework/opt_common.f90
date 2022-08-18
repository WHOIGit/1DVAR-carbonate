!----------------------------------------------------------------------------
!     CVS:$Id: opt_common.F90,v 1.0 2008/05/09 10:58 LUO
!     CVS:$Name:  $
!----------------------------------------------------------------------------
module opt_common
  implicit none

  ! save everything
  save

  double precision :: RHOBEG,RHOEND
  
  integer :: NPT, IPRINT,MAXFUN
  
  double precision, allocatable :: W(:)
  
  contains

  subroutine opt_init
    
    use common_mod, only : optpar_fname
    
    use eco_params, only : nparams_opt
    
    implicit none

    integer, parameter :: iun=47

!    NPT is the number of interpolation conditions. Its value must be in the
!    interval [N+2,(N+1)(N+2)/2].
    NPT = 2 * nparams_opt + 1
!   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!   amount of printing. Specifically, there is no output if IPRINT=0 and
!   there is output only at the return if IPRINT=1. Otherwise, each new
!   value of RHO is printed, with the best vector of variables so far and
!   the corresponding value of the objective function. Further, each new
!   value of F with its variables are output if IPRINT=3.  
    IPRINT = 3
!    The array W will be used for working space. Its length must be at least
!    (NPT+11)*(NPT+N)+N*(3*N+11)/2.    
    allocate (W((NPT+11)*(NPT+nparams_opt)+ &
                     nparams_opt*(3*nparams_opt+11)/2 + 100))
!   MAXFUN must be set to an upper bound on the number of calls of CALFUN. 
!   RHOBEG and RHOEND must be set to the initial and final values of a trust
!   region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!   RHOBEG should be about one tenth of the greatest expected change to a
!   variable, and RHOEND should indicate the accuracy that is required in
!    the final values of the variables. 
    namelist /opt_parms_nml/ &
         MAXFUN, &
         RHOBEG, &
         RHOEND

    !-----------------------------------------------------------------------
    !   default namelist settings
    !-----------------------------------------------------------------------
         MAXFUN = 3000
  
         RHOBEG = 0.5
  
         RHOEND = 1d-2

    !-----------------------------------------------------------------------
    !   read in parameter values
    !-----------------------------------------------------------------------
    open(unit=iun, file=trim(optpar_fname), status='old')
    read(unit=iun, nml=opt_parms_nml)
    close(unit=iun)

  end subroutine opt_init

end module opt_common
