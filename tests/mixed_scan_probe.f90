program mixed_scan_probe
  use vars
  use functions, only: read_in_params,read_scan_input,read_guess_input
  use disprels, only: refine_guess,om_scan,om_double_scan
  use scan_support, only: prepare_scan_start
  use adaptive_kinetic, only: finish_bessel_diagnostics
  use bessel_cache, only: write_bessel_stats
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  integer :: i,j
  integer(int64),allocatable :: before(:),after(:)
  real,allocatable :: state(:)
  call read_in_params
  allocate(state(7*nspec+2))
  call physical_state()
  before=transfer(state,[0_int64],size(state))
  call read_scan_input
  call prepare_scan_start()
  call read_guess_input
  call refine_guess
  if (option==1) then
     do i=1,nscan
        call om_scan(i)
        do j=1,nroot_max
           write(*,'(a,2i4,2es26.17e3)') 'FINAL_ROOT ',i,j,wroots(:,j)
        enddo
        call write_bessel_stats(6)
     enddo
  else if (option==2) then
     call om_double_scan
     call write_bessel_stats(6)
  else
     error stop 'Unsupported probe option'
  endif
  call physical_state()
  after=transfer(state,[0_int64],size(state))
  if (any(before/=after)) error stop 'Scan changed physical plasma state'
  if (has_cold_reference() .and. allocated(solver_spec)) error stop 'Cold scan allocated thermal reference state'
  write(*,'(a)') 'PHYSICAL_STATE_UNCHANGED'
  call finish_bessel_diagnostics
contains
  subroutine physical_state()
    integer :: j
    state(1:2)=[betap,vtp]
    do j=1,nspec
       state(3+7*(j-1):2+7*j)=[spec(j)%tau_s,spec(j)%mu_s,spec(j)%alph_s,spec(j)%Q_s, &
            spec(j)%D_s,spec(j)%vv_s,spec(j)%nu_ns]
    enddo
  end subroutine
end program
