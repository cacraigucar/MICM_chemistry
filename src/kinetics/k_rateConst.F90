module k_rateConst
!--------------------------------------------------------
! k rate constants for 3component chemistry
!--------------------------------------------------------

use machine,         only: r8 => kind_phys

implicit none
private

public :: k_rateConst_init
public :: k_rateConst_run
public :: k_rateConst_finalize

! k_rateConst are computed at the beginning of the 
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! k_rateConst will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

contains

!> \section arg_table_k_rateConst_init Argument Table
!! [ errmsg ]
!!    standard_name = ccpp_error_message
!!    long_name = CCPP error message
!!    units = 1
!!    dimensions = ()
!!    type = character
!!    kind = len=512
!!    intent = out
!!    optional = F
!! [ errflg ]
!!    standard_name = ccpp_error_flag
!!    long_name = CCPP error flag
!!    units = flag
!!    dimensions = ()
!!    type = integer
!!    intent = out
!!    optional = F
!!
  subroutine k_rateConst_init(errmsg, errflg)
      
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg=''
    errflg=0

    ! Nothing for the 3component chemistry to do at init time currently

  end  subroutine k_rateConst_init

  !---------------------------
  ! Compute k_rateConst, given M, P, T
  ! Execute once for the chemistry-time-step advance
  !---------------------------
!> \section arg_table_k_rateConst_run Argument Table
!! [ k_rateConst ]
!!    standard_name = gasphase_rate_constants
!!    long_name = k rate constants
!!    units = s-1
!!    dimensions = (number_of_kinetics_reactions)
!!    type = real
!!    kind = kind_phys
!!    intent = out
!!    optional = F
!! [ c_m ]
!!    standard_name = total_number_density
!!    long_name = total number density
!!    units = molecules cm-3
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ rh ]
!!    standard_name = water_vapor_relative_humidity
!!    long_name = relative humidity
!!    units = percent
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ c_h2o ]
!!    standard_name = water_vapor_mole_fraction
!!    units = mole/mole
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ TEMP ]
!!    standard_name = temperature
!!    long_name = mid-point layer temperature
!!    units = K
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ errmsg ]
!!    standard_name = ccpp_error_message
!!    long_name = CCPP error message
!!    units = none
!!    dimensions = ()
!!    type = character
!!    kind = len=512
!!    intent = out
!!    optional = F
!! [ errflg ]
!!    standard_name = ccpp_error_flag
!!    long_name = CCPP error flag
!!    units = flag
!!    dimensions = ()
!!    type = integer
!!    intent = out
!!    optional = F
!!
  subroutine k_rateConst_run(k_rateConst, c_m, rh, c_h2o, temp, errflg, errmsg)
  
    real(r8),           intent(inout) :: k_rateConst(:)
    real(r8),           intent(in)  :: c_m
    real(r8),           intent(in)  :: rh 
    real(r8),           intent(in)  :: c_h2o     
    real(r8),           intent(in)  :: TEMP
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg


    ! retrieve the temperature used by tuv (the photolysis level)

    errmsg=''
    errflg=0
  
! These are probably set by the Chemistry Cafe
#include "k_rateConst.inc"


  end subroutine k_rateConst_run
  
  subroutine k_rateConst_finalize
  end subroutine k_rateConst_finalize

#include "rate_functions.inc"
  
end module k_rateConst
