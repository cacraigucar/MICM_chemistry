!--------------------------------------------------------------------------------
! utility to compute background air mass quantities
!--------------------------------------------------------------------------------
module mass_quantities_util
  use machine, only : rk => kind_phys
  use const_props_mod, only : const_props_type
  implicit none
  
  integer :: o2_ndx, n2_ndx
  real(rk), allocatable :: molar_mass(:)
  
contains

!> \section arg_table_mass_quantities_util_init Argument Table
!! [ cnst_info ]
!!    standard_name = chemistry_constituent_info
!!    long_name = chemistry_constituent_info
!!    units = DDT
!!    dimensions = (number_of_constituents)
!!    type = const_props_type
!!    intent = in
!!    optional = F
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
  subroutine mass_quantities_util_init(cnst_info, errmsg, errflg)
    type(const_props_type), intent(in)  :: cnst_info(:)
    character(len=512),     intent(out) :: errmsg
    integer,                intent(out) :: errflg

    integer :: ncnst, i

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    o2_ndx =-1
    n2_ndx =-1

    ncnst = size(cnst_info)
    allocate(molar_mass(ncnst))

    do i = 1,ncnst
       molar_mass(i) = cnst_info(i)%get_wght()
       if ( cnst_info(i)%get_name() == 'O2' ) o2_ndx = i
       if ( cnst_info(i)%get_name() == 'N2' ) n2_ndx = i
    end do

  end subroutine mass_quantities_util_init

!> \section arg_table_mass_quantities_util_run Argument Table
!! [ press ]
!!    standard_name = air_pressure
!!    long_name = ambient pressure
!!    units = Pa
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ temp ]
!!    standard_name = temperature
!!    long_name = ambient temperature
!!    units = K
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ vmr ]
!!    standard_name = concentration
!!    long_name = species concentration
!!    units = mole mole-1
!!    dimensions = (number_of_constituents)
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ density ]
!!    standard_name = total_number_density
!!    long_name = total number density
!!    units = molecules cm-3
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = out
!!    optional = F
!! [ mbar ]
!!    standard_name = mean_molec_mass
!!    long_name = mean molecular mass
!!    units = g mole-1
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = out
!!    optional = F
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
  subroutine mass_quantities_util_run( press, temp, vmr, density, mbar, errmsg, errflg )

    real(rk), intent(in)            :: press
    real(rk), intent(in)            :: temp
    real(rk), intent(in)            :: vmr(:)
    real(rk), intent(out)           :: density
    real(rk), intent(out)           :: mbar
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    real(rk) :: n2_vmr
    real(rk), parameter :: molar_mass_n2 = 28.0134_rk ! g/mole
    real(rk), parameter :: kboltz= 1.38064852e-16_rk ! boltzmann constant (erg/K)
   
    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    density = 10._rk*press/(kboltz*temp)

    if (o2_ndx>0) then
       mbar = sum( vmr(:)*molar_mass(:) )
       if (n2_ndx<-1) then
          n2_vmr = 1._rk - sum(vmr(:))
          mbar = mbar + n2_vmr*molar_mass_n2
       endif
    else
       mbar = 28.966_rk ! set to constant if the major species VMRs are unknown
    endif
    
  end subroutine mass_quantities_util_run
  
end module mass_quantities_util
