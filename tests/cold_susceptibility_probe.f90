program cold_susceptibility_probe
  !! Direct-call probe: no input reader, thermal adapter or dispersion solver.
  use nrtype, only: SP, SPC
  use vars, only: specie, betap, vtp, normalization, spec, solver_spec
  use cold_susceptibility, only: calc_chi_cold
  implicit none
  type(specie) :: species
  real(SP) :: wr, wi, kpar, kperp, va
  complex(SPC) :: chi(3,3)
  integer :: ierr, ios, row, column

  normalization = 'd_ref'
  betap = 0.0_SP
  vtp = 0.0_SP
  ! Neither global species array is needed, even for a fully cold reference.
  if (allocated(spec) .or. allocated(solver_spec)) stop 1
  write(*,'(a,2i4)') 'PRECISION ',precision(0.0_SP),storage_size(0.0_SP)
  do
     read(*,*,iostat=ios) wr, wi, kpar, kperp, va, species%mu_s, &
          species%Q_s, species%D_s, species%vv_s, &
          species%tau_s, species%alph_s, species%nu_ns
     if (ios < 0) exit
     if (ios /= 0) stop 2
     call calc_chi_cold(species, cmplx(wr,wi,kind=SPC), kpar, kperp, va, chi, ierr)
     write(*,'(a,i4)') 'STATUS ',ierr
     do row=1,3
        do column=1,3
           write(*,'(2es27.17e3)') real(chi(row,column),SP),aimag(chi(row,column))
        enddo
     enddo
  enddo
end program cold_susceptibility_probe
