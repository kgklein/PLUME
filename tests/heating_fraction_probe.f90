program heating_fraction_probe
  use heating_energy
  implicit none
  type(heating_result) :: r
  real :: electric,magnetic,one,two,uncertainty
  integer :: ios
  do
     read(*,*,iostat=ios) electric,magnetic,one,two,uncertainty
     if (ios/=0) exit
     call reset_heating(r,2,0)
     r%electric=electric;r%magnetic=magnetic;r%species=[one,two];r%uncertainty=uncertainty
     call finish_energy(r)
     write(*,'(i3,9es26.17e3)') r%status,r%total,r%f_e,r%f_b,r%f_species,r%f_e_abs,r%fractions,r%uncertainty,r%species_total
  enddo
end program
