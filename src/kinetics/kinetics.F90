
module kinetics

  use kinetics_module, only: kinetics_type
  use const_props_mod, only: const_props_type
  use machine,         only: rk => kind_phys
  
  implicit none

  private
  public :: kinetics_init 
  public :: kinetics_run
  public :: kinetics_finalize


  
contains

!> \section arg_table_kinetics_init Argument Table
!! [ nkRxt ]
!!    standard_name = number_of_kinetics_reactions
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
!! [ njRxt ]
!!    standard_name = number_of_photolysis_reactions
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
!! [ nTotRxt ]
!!    standard_name = total_number_of_chemical_reactions
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
!! [ ncnst ]
!!    standard_name = number_of_constituents
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
!! [ nSpecies ]
!!    standard_name = number_of_prognostic_chemistry_constituents
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
!! [ cnst_info ]
!!    standard_name = chemistry_constituent_info
!!    units = DDT
!!    dimensions = (number_of_constituents)
!!    type = const_props_type
!!    intent = out
!! [ model_name ]
!!    standard_name = chemistry_model_name
!!    units = 1
!!    dimensions = ()
!!    type = character
!!    kind = len=80
!!    intent = out
!! [ theKinetics ]
!!    standard_name = kinetics_data
!!    long_name = chemistry kinetics
!!    units = DDT
!!    dimensions = ()
!!    type = kinetics_type
!!    intent = out
!! [ errmsg ]
!!    standard_name = ccpp_error_message
!!    long_name = CCPP error message
!!    units = 1
!!    dimensions = ()
!!    type = character
!!    kind = len=512
!!    intent = out
!! [ errflg ]
!!    standard_name = ccpp_error_flag
!!    long_name = CCPP error flag
!!    units = flag
!!    dimensions = ()
!!    type = integer
!!    intent = out
!!
  subroutine kinetics_init( nkRxt, njRxt, nTotRxt, ncnst, nSpecies, cnst_info, model_name, theKinetics, errmsg, errflg )

    !--- arguments
    integer,                        intent(out) :: nkRxt
    integer,                        intent(out) :: njRxt
    integer,                        intent(out) :: nTotRxt
    integer,                        intent(out) :: ncnst
    integer,                        intent(out) :: nSpecies
    type(const_props_type), intent(out) :: cnst_info(:)
    character(len=*),               intent(out) :: model_name
    type(kinetics_type),            intent(out) :: theKinetics
    character(len=512),             intent(out) :: errmsg
    integer,                        intent(out) :: errflg

    character(len=120) :: jsonfile 

#include "model_name.inc"

    ! Read in the constituent information from the appropriate json file. 
    jsonfile = '../../MICM_chemistry/generated/'//trim(model_name)//'/molec_info.json'
    call json_loader_read( jsonfile, nSpecies, cnst_info )

    call theKinetics%rateConst_init( nTotRxt )

    errmsg = ''
    errflg = 0

  end subroutine kinetics_init

!> \section arg_table_kinetics_run Argument Table
!! [ theKinetics ]
!!    standard_name = kinetics_data
!!    long_name = chemistry kinetics
!!    units = DDT
!!    dimensions = ()
!!    type = kinetics_type
!!    intent = inout
!!    optional = F
!! [ k_rateConst ]
!!    standard_name = gasphase_rate_constants
!!    long_name = gas phase rates constants
!!    units = s-1
!!    dimensions = (:)
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ j_rateConst ]
!!    standard_name = photo_rate_constants
!!    long_name = photochemical rates constants
!!    units = s-1
!!    dimensions = (:)
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ c_m ]
!!    standard_name = total_number_density
!!    long_name = total number density
!!    units = molecules/cm3
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
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
  subroutine kinetics_run( theKinetics, k_rateConst, j_rateConst, c_m, errmsg, errflg )

    !--- arguments
    type(kinetics_type), pointer, intent(inout)      :: theKinetics
    real(rk),           intent(in)    :: k_rateConst(:)
    real(rk),           intent(in)    :: j_rateConst(:)
    real(rk),           intent(in)    :: c_m ! total number density
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    !--- local variables
    integer :: Ierr

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !--- set the gas phase rate constants
    call theKinetics%rateConst_update( k_rateConst, j_rateConst, c_m)

  end subroutine kinetics_run

!> \section arg_table_kinetics_finalize Argument Table
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
  subroutine kinetics_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    errmsg = ''
    errflg = 0

  end subroutine kinetics_finalize

end module kinetics
