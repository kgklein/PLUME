program kinetic_core_probe
  !! Isolated direct susceptibility experiment. Never calls input validation for
  !! unsupported physical states, never calls roots, never reports a mixed det.
  use vars, only: z_backend
  use vars, only: bessel_mode
  use vars, only: specie,nspec,spec,susc,lam,low_n,pi,betap,kpar,kperp,vtp,normalization
  use vars, only: solver_betap,solver_vtp,solver_kpar,solver_kperp,solver_spec,sync_solver_state
  use auxiliary_normalization
  use auxiliary_kinetic_kernel, only: probe_disp0
  use disprels, only: disp,bessel,zet_in
  use kinetic_core, only: kinetic_normalization,calc_chi_kinetic
  use cold_susceptibility, only: invalid_complex,finite_complex
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(auxiliary_scale) :: scale
  type(kinetic_normalization) :: core_state
  type(specie) :: converted
  complex :: w, unused
  real :: wr,wi,b_star,h,nan
  integer :: i,j,k,ios,mode,status,species_status
  integer(int64), allocatable :: physical_bits(:)
  integer(int64) :: global_bits(4)

  z_backend='legacy'
  bessel_mode='fixed' ! Historical Stage 3B comparison, explicit harmonic convention.
  read(*,*) nspec,kpar,kperp,vtp
  allocate(spec(nspec),solver_spec(nspec),susc(nspec,3,3))
  do i=1,nspec
     read(*,*) spec(i)%tau_s,spec(i)%mu_s,spec(i)%alph_s,spec(i)%Q_s,spec(i)%D_s,spec(i)%vv_s
     spec(i)%nu_ns=0.
  enddo
  betap=spec(1)%tau_s
  normalization='d_ref'
  low_n=.false.
  pi=4.*atan(1.)
  ! Bit comparisons also permit deliberately poisoned cold anisotropy fields.
  physical_bits=transfer(spec,[0_int64],size=7*nspec)
  global_bits=transfer([betap,kpar,kperp,vtp],global_bits)
  nan=real(invalid_complex())
  do
     read(*,*,iostat=ios) wr,wi,b_star,mode
     if (ios<0) exit
     if (ios/=0) stop 1
     w=cmplx(wr,wi)
     status=0
     if (mode==1) then
        ! Independent production baseline, valid warm-reference inputs only.
        unused=disp(w)
     else if (mode==2) then
        ! Original beta and reference anisotropy; extracted kernel versus actual
        ! production body in exactly the same internal representation.
        call sync_solver_state
        unused=probe_disp0(w,spec(1)%alph_s)
     else
        solver_betap=nan; solver_vtp=nan; solver_kpar=nan; solver_kperp=nan
        call make_auxiliary_scale(kpar,kperp,vtp,b_star,scale,status)
        if (status==0) then
           call scale%export(solver_betap,solver_vtp,solver_kpar,solver_kperp,h)
           do i=1,nspec
              call auxiliary_species(spec(i),scale,converted,species_status)
              if (spec(i)%tau_s==0.) then
                 ! A real cold population is not an auxiliary population. Its
                 ! susceptibility is outside this experiment; no cold tensor or
                 ! total dielectric/determinant is emitted by this probe.
                 if (species_status/=auxiliary_not_warm) error stop 'Cold species converted'
                 solver_spec(i)=spec(i)
              else
                 status=species_status
                 if (status/=0) exit
                 solver_spec(i)=converted
              endif
           enddo
        endif
        if (status==0) then
           ! Cold rows must be initialized for disp0's unused accumulation only.
           susc=cmplx(0.,0.)
           if (mode==3) then
              core_state=kinetic_normalization(solver_betap,solver_vtp,solver_kpar,solver_kperp,h)
              do i=1,nspec
                 if (spec(i)%tau_s==0.) cycle
                 call calc_chi_kinetic(solver_spec(i),w,core_state,susc(i,:,:),status,bessel,zet_in)
                 if (status/=0) exit
              enddo
           else
              unused=probe_disp0(w,h)
           endif
        endif
     endif
     if (status==0) then
        do i=1,nspec
           if (spec(i)%tau_s<=0.) cycle
           do j=1,3
              do k=1,3
                 if (.not. finite_complex(susc(i,j,k))) status=4
              enddo
           enddo
        enddo
     endif
     if (any(transfer(spec,[0_int64],size=7*nspec)/=physical_bits)) &
          error stop 'Physical species state changed'
     if (any(transfer([betap,kpar,kperp,vtp],global_bits)/=global_bits)) &
          error stop 'Physical globals changed'
     write(*,'(a,i4)') 'STATUS ',status
     if (status/=0) cycle
     do i=1,nspec
        if (spec(i)%tau_s<=0.) cycle
        do j=1,3
           do k=1,3
              write(*,'(a,3i4,2es26.17e3)') 'CHI ',i,j,k,real(susc(i,j,k)),aimag(susc(i,j,k))
           enddo
        enddo
     enddo
  enddo
end program kinetic_core_probe
