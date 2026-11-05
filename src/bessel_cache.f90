module bessel_cache
  !! Per-species, in-memory coefficient storage. Frequency never enters a key.
  use harmonic_selection
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: bessel_entry, fetch_bessel, clear_bessel_cache, cache, write_bessel_stats, cap_events
  type :: bessel_entry
     real :: argument=-1., tolerance=-1., tail0=1., tail2=1.
     integer :: active=0, capacity=-1, filled=-1, limit=-1, hits=0, misses=0, evaluations=0
     logical :: valid=.false., converged=.false.
     character(16) :: mode=''
     real, allocatable :: jn(:,:),jpn(:,:)
     complex, allocatable :: tsi(:),zz(:)
  end type
  type(bessel_entry), allocatable, target :: cache(:)
  type :: cap_record
     integer :: species,limit
     real :: argument,tolerance,e0,e2
  end type
  type(cap_record), allocatable :: cap_events(:)
contains
  subroutine clear_bessel_cache()
    if (allocated(cache)) deallocate(cache)
    if (allocated(cap_events)) deallocate(cap_events)
  end subroutine

  subroutine fetch_bessel(species,count,argument,mode,tolerance,limit,fixed,entry,ierr)
    integer, intent(in) :: species,count,limit,fixed
    real, intent(in) :: argument,tolerance
    character(*), intent(in) :: mode
    type(bessel_entry), pointer, intent(out) :: entry
    integer, intent(out) :: ierr
    type(bessel_entry), allocatable :: resized(:)
    type(bessel_entry) :: grown
    real :: e0,e2
    real, allocatable :: coefficients(:)
    integer :: n,capacity,j,stat,requested
    logical :: converged
    nullify(entry)
    ierr=1
    if (species<1 .or. species>count .or. count<1) return
    if (.not.ieee_is_finite(argument) .or. argument<0.) return
    if (.not.ieee_is_finite(tolerance) .or. tolerance<=0. .or. tolerance>=1.) return
    if (limit<1 .or. fixed<1) return
    if (mode/='adaptive' .and. mode/='fixed') return
    if (.not.allocated(cache)) allocate(cache(count))
    if (size(cache)/=count) then
       allocate(resized(count))
       resized(1:min(count,size(cache)))=cache(1:min(count,size(cache)))
       call move_alloc(resized,cache)
    endif
    entry=>cache(species)
    requested=limit
    if (mode=='fixed') requested=fixed
    if (entry%valid .and. entry%argument==argument .and. entry%mode==mode .and. &
         entry%tolerance==tolerance .and. entry%limit==requested) then
       entry%hits=entry%hits+1
       ierr=0
       return
    endif
    entry%misses=entry%misses+1
    call select_harmonics(argument,tolerance,requested,n,e0,e2,converged,ierr)
    if (ierr/=0) return
    if (mode=='fixed') then
       n=fixed
       call harmonic_tail_bounds(argument,n,e0,e2)
       converged=max(e0,e2)<tolerance
    endif
    capacity=entry%capacity
    if (n>capacity) capacity=max(n,capacity+max(8,capacity/2))
    ! A changed tolerance can reuse all coefficients at the same physical
    ! argument, even if it changes the active cutoff. Storage never shrinks.
    if (.not.entry%valid .or. entry%argument/=argument .or. n>entry%filled) then
       allocate(coefficients(0:n+1),stat=stat)
       ierr=5
       if (stat/=0) return
       call scaled_bessel_array(argument,n+1,coefficients,ierr)
       if (ierr/=0) return
       if (capacity>entry%capacity) then
          ! Allocate first, then move: a failed growth leaves the old entry
          ! usable for its original key and cannot leave partially allocated storage.
          allocate(grown%jn(-capacity-1:capacity+1,1),grown%jpn(-capacity:capacity,1), &
               grown%tsi(-capacity:capacity),grown%zz(-capacity:capacity),stat=stat)
          if (stat/=0) then
             ierr=5
             return
          endif
          call move_alloc(grown%jn,entry%jn)
          call move_alloc(grown%jpn,entry%jpn)
          call move_alloc(grown%tsi,entry%tsi)
          call move_alloc(grown%zz,entry%zz)
          entry%capacity=capacity
       endif
       do j=-n-1,n+1
          entry%jn(j,1)=coefficients(abs(j))
       enddo
       do j=-n,n
          entry%jpn(j,1)=.5*(entry%jn(j-1,1)+entry%jn(j+1,1))
       enddo
       entry%evaluations=entry%evaluations+1
       entry%filled=n
    endif
    entry%argument=argument; entry%tolerance=tolerance; entry%limit=requested
    entry%mode=mode; entry%active=n; entry%tail0=e0; entry%tail2=e2
    entry%converged=converged; entry%valid=.true.; ierr=0
    if (mode=='adaptive' .and. .not.converged) then
       if (.not.allocated(cap_events)) allocate(cap_events(0))
       cap_events=[cap_events,cap_record(species,limit,argument,tolerance,e0,e2)]
       write(*,'(a,i0,a,es12.4,a,es12.4,a,i0,a,2es12.4)') &
            'WARNING: Bessel cap: species ',species,', lambda=',argument,', tolerance=',tolerance, &
            ', maximum=',limit,', tail bounds=',e0,e2
    endif
  end subroutine

  subroutine scaled_bessel_array(x,last,values,ierr)
    !! Miller's downward ratio recurrence, normalized by the positive sum.
    !! Unlike repeated legacy bessim(n,x), its starting order resolves the full
    !! distribution, including when x is large and low-order values are small.
    real, intent(in) :: x
    integer, intent(in) :: last
    real, intent(out) :: values(0:last)
    integer, intent(out) :: ierr
    real, allocatable :: ratios(:)
    real :: e0,e2,product,total,compensation,y,t
    integer :: m,n,stat,start
    logical :: converged
    ierr=0; values=0.
    if (x==0.) then
       values(0)=1.
       return
    endif
    ! Extra recurrence depth below the double-precision tail, rather than a
    ! fixed-order seed that fails for broad distributions. Range errors are
    ! numerical failures, never a silently substituted physical argument.
    start=max(last,32)
    do
       call harmonic_tail_bounds(x,start,e0,e2)
       if (max(e0,e2)<epsilon(x)*.01) exit
       if (start>(huge(start)-32)/2) then
          ierr=5
          return
       endif
       start=2*start
    enddo
    call select_harmonics(x,epsilon(x)*.01,start,m,e0,e2,converged,ierr)
    if (ierr/=0) return
    m=max(last,m)+32
    allocate(ratios(1:m+1),stat=stat)
    if (stat/=0) then
       ierr=5
       return
    endif
    ratios(m+1)=0.
    do n=m,1,-1
       ratios(n)=x/(2.*real(n)+x*ratios(n+1))
    enddo
    product=1.; total=1.; compensation=0.; values(0)=1.
    do n=1,m
       product=product*ratios(n)
       if (n<=last) values(n)=product
       y=2.*product-compensation
       t=total+y
       compensation=(t-total)-y
       total=t
    enddo
    values=values/total
    if (.not.all(ieee_is_finite(values))) ierr=5
  end subroutine

  subroutine write_bessel_stats(unit)
    integer, intent(in) :: unit
    integer :: j
    write(unit,'(a)') '# species lambda N capacity converged E0_bound E2_bound hits misses array_evaluations'
    if (.not.allocated(cache)) return
    do j=1,size(cache)
       if (.not.cache(j)%valid) cycle
       write(unit,'(i6,es24.15,2i9,l3,2es24.15,3i12)') j,cache(j)%argument, &
            cache(j)%active,cache(j)%capacity,cache(j)%converged,cache(j)%tail0,cache(j)%tail2, &
            cache(j)%hits,cache(j)%misses,cache(j)%evaluations
    enddo
    if (allocated(cap_events)) then
       write(unit,'(a)') '# unresolved history: species maximum lambda tolerance E0_bound E2_bound'
       do j=1,size(cap_events)
          write(unit,'(2i9,4es24.15)') cap_events(j)%species,cap_events(j)%limit,cap_events(j)%argument, &
               cap_events(j)%tolerance,cap_events(j)%e0,cap_events(j)%e2
       enddo
    endif
  end subroutine
end module bessel_cache
