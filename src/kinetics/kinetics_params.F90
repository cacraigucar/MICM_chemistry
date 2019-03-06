module kinetics_params

 use machine, only: kind_phys

 integer :: ncnst
 integer :: nkRxt
 integer :: njRxt
 integer :: nTotRxt

 real(kind_phys) :: total_dens

end module kinetics_params
