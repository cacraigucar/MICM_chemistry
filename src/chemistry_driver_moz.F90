!-----------------------------------------------------------------------------------------------
! chemistry driver which uses MOZART solver method
!-----------------------------------------------------------------------------------------------
module chemistry_driver_moz

use kinetics_module,  only       : kinetics_type
use kinetics,         only       : kinetics_init, kinetics_run
use machine,          only       : r8 => kind_phys
use const_props_mod,  only       : const_props_type
use Mozart_Solver,    only: MozartSolver

implicit none

type(MozartSolver) :: theSolver
type(kinetics_type), allocatable :: theKinetics 

contains

!> \section arg_table_chemistry_driver_moz_init Argument Table
!! [ Time ]
!!    standard_name = model_time
!!    units = s
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ dt ]
!!    standard_name = time_step_for_physics
!!    long_name = time_step_for_physics
!!    units = s
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ ncnst ]
!!    standard_name = number_of_constituents
!!    units = count
!!    dimensions = ()
!!    type = integer
!!    intent = out
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
subroutine chemistry_driver_moz_init(Time, dt, ncnst, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  real(r8), intent(in)            :: Time
  real(r8), intent(in)            :: dt
  integer,  intent(out)           :: ncnst      ! number prognostic constituents
  character(len=512),intent(out)  :: errmsg
  integer, intent(out)            :: errflg          ! error index from CPF

  type(const_props_type),allocatable :: cnst_info(:)
  real(r8)                        :: TimeStart, TimeEnd
  integer                         :: nTotRxt    ! total number of chemical reactions

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(r8) :: rcntrl(20)     ! real control array for ODE solver
  real(r8), allocatable :: absTol(:), relTol(:)
  character(len=40) :: model_name
  
  write(0,*) ' Entered chemistry_driver_moz_init'
  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  TimeStart = Time
  TimeEnd = Time + dt

!-----------------------------------------------------------
!  initialize the kinetics
!-----------------------------------------------------------
  allocate( theKinetics )
  call kinetics_init(nTotRxt, ncnst, theKinetics, cnst_info, errmsg, errflg)

!-----------------------------------------------------------
!  initialize ode solver "control" variable defaults
!-----------------------------------------------------------
  allocate(absTol(ncnst))
  allocate(relTol(ncnst))

  absTol(:) = 1.e-8_r8
  relTol(:) = 1.e-3_r8
  icntrl(:) = 0
  rcntrl(:) = 0._r8

!-----------------------------------------------------------
!  set ode solver "control" variables for MOZART solver
!-----------------------------------------------------------
  icntrl(1) = 1                                 ! autonomous, F depends only on Y
  rcntrl(2) = dt                                ! Hmax
  rcntrl(3) = .01_r8*dt                         ! Hstart

  write(*,*) ' '
  write(*,*) 'icntrl settings'
  write(*,'(10i6)') icntrl(1:10)
  write(*,*) 'rcntrl settings'
  write(*,'(1p,10(1x,g0))') rcntrl(1:10)
  write(*,*) ' '

  call theSolver%Initialize( Tstart=TimeStart, Tend=TimeEnd, AbsTol=AbsTol, RelTol=RelTol, &
                                       ICNTRL=icntrl, RCNTRL=rcntrl, Ierr=errflg )


end subroutine chemistry_driver_moz_init

!> \section arg_table_chemistry_driver_moz_run Argument Table
!! [ vmr ]
!!    standard_name = concentration
!!    long_name = species concentration
!!    units = mole mole-1
!!    dimensions = (number_of_constituents)
!!    type = real
!!    kind = kind_phys
!!    intent = inout
!!    optional = F
!! [ Time]
!!    standard_name = model_time
!!    units = s
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ dt ]
!!    standard_name = time_step_for_physics
!!    long_name = time_step_for_physics
!!    units = s
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ j_rateConst ]
!!    standard_name = photo_rate_constants
!!    long_name = photochemical rates constants
!!    units = s-1
!!    dimensions = (horizontal_dimension)
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ k_rateConst ]
!!    standard_name = gasphase_rate_constants
!!    long_name = k rate constants
!!    units = s-1
!!    dimensions = (horizontal_dimension)
!!    type = real
!!    kind = kind_phys
!!    intent = in
!!    optional = F
!! [ c_m ]
!!    standard_name = total_number_density
!!    long_name = total number density
!!    units = molecules cm-3
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
subroutine chemistry_driver_moz_run(vmr, Time, dt, j_rateConst,  k_rateConst, c_m, errmsg, errflg)


  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------

  real(kind=r8), intent(inout)    :: vmr(:)                ! "working" concentration passed thru CPF
  real(r8), intent(in)            :: Time
  real(r8), intent(in)            :: dt
  real(kind=r8), intent(in)       :: j_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=r8), intent(in)       :: k_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=r8), intent(in)       :: c_m                   ! total number density
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF

  integer            :: i, k, n
  real(r8)                        :: TimeStart, TimeEnd

  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  TimeStart = Time
  TimeEnd = Time + dt

!-----------------------------------------------------------
!  update the kinetics
!-----------------------------------------------------------
  call kinetics_run(theKinetics, k_rateConst, j_rateConst, c_m, errmsg, errflg)
  if (errflg /= 0) return

  call theKinetics%rateConst_print()

!-----------------------------------------------------------
!  solve the current timestep's chemistry
!-----------------------------------------------------------
  call theSolver%Run( Tstart=TimeStart, Tend=TimeEnd, y=vmr, theKinetics=theKinetics, Ierr=errflg  )

  if (errflg /= 0) then
     errmsg = 'ERROR: theSolver%Run'
     return
  end if
  
end subroutine chemistry_driver_moz_run

!> \section arg_table_chemistry_driver_moz_finalize Argument Table
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
subroutine chemistry_driver_moz_finalize(errmsg,errflg)
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF
  errmsg = ''
  errflg = 0
  
  deallocate( theKinetics )
  
end subroutine chemistry_driver_moz_finalize

end module chemistry_driver_moz

