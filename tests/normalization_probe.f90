program normalization_probe
  use vars, only: nspec, susc, lam, low_n, spec, betap, kperp, kpar, vtp
  use functions, only: read_in_params
  use disprels, only: disp, calc_eigen, rtsec
  implicit none
  complex :: omega, determinant, root, ef(3), bf(3)
  complex, allocatable :: density(:), velocity(:,:)
  real, allocatable :: heating(:), split(:,:)
  integer :: i,j,k,iflag

  call read_in_params
  allocate(density(nspec),velocity(3,nspec),heating(nspec),split(7,nspec))
  omega=cmplx(0.012,-0.001)
  determinant=disp(omega)
  write(*,'(a,2es25.16)') 'DET ',real(determinant),aimag(determinant)
  do i=1,nspec
     do j=1,3
        do k=1,3
           write(*,'(a,3i4,2es25.16)') 'CHI ',i,j,k,real(susc(i,j,k)),aimag(susc(i,j,k))
        enddo
     enddo
  enddo
  call calc_eigen(omega,ef,bf,velocity,density,heating,split,.true.,.true.)
  write(*,'(a,2es25.16)') 'POL ', &
       2.*aimag(ef(1)*conjg(ef(2)))/(abs(ef(1))**2+abs(ef(2))**2), &
       abs(ef(3))/sqrt(sum(abs(ef)**2))
  do j=1,3
     write(*,'(a,i4,4es25.16)') 'FIELD ',j,real(ef(j)),aimag(ef(j)),real(bf(j)),aimag(bf(j))
  enddo
  do i=1,nspec
     write(*,'(a,i4,2es25.16)') 'HEAT ',i,heating(i),split(5,i)
     do j=1,7
        write(*,'(a,2i4,es25.16)') 'SPLIT ',i,j,split(j,i)
     enddo
     write(*,'(a,i4,2es25.16)') 'DENS ',i,real(density(i)),aimag(density(i))
     do j=1,3
        write(*,'(a,2i4,2es25.16)') 'FLOW ',i,j,real(velocity(j,i)),aimag(velocity(j,i))
     enddo
  enddo
  root=rtsec(disp,cmplx(0.01,-0.001),cmplx(0.011,-0.001),1.e-13,iflag)
  write(*,'(a,2es25.16,i5)') 'ROOT ',real(root),aimag(root),iflag
  determinant=disp(root)
  write(*,'(a,es25.16)') 'RESID ',abs(determinant)/maxval(abs(lam))**3
end program normalization_probe
