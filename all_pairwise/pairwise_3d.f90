program main
  implicit none

  integer, parameter :: dp = 8
  integer, parameter :: dimEspace = 3
  integer, parameter :: nbPoints = 1000
  integer, parameter :: nbTimes = 20
  integer, parameter :: nb = 8
  
  real(dp), allocatable, dimension(:,:) :: points
  real(dp), allocatable, dimension(:,:) :: distance
  integer :: i, j, k, ii, jj
  ! integer :: beginning, end, rate
  real(dp) :: beginning, end
  real(dp) :: time, sum_squares, sum_time
  
  allocate(points(1:dimEspace, 1:nbPoints))
  allocate(distance(1:nbPoints, 1:nbPoints))
  call random_number(points)

  do j = 1, nbPoints
     do i = 1, nbPoints
        distance(i, j) =     &
             (points(1, i) - points(1, j))**2 + &
             (points(2, i) - points(2, j))**2 + &
             (points(3, i) - points(3, j))**2  
     end do
  end do
  
  sum_squares=0.0
  sum_time = 0.0
  ! call system_clock(COUNT_RATE=rate)
  do k = 1, nbTimes
     ! call system_clock(beginning)
     call cpu_time(beginning)
     do j = 1, nbPoints
        do i = 1, nbPoints
           distance(i, j) =      &
                (points(1, i) - points(1, j))**2 + &
                (points(2, i) - points(2, j))**2 + &
                (points(3, i) - points(3, j))**2  
        end do
     end do
     ! call system_clock(end)
     call cpu_time(end)
     ! time = 1e-3*real(end - beginning) / real(rate)
     time = 1e3*(end - beginning)
     write(*,*) 'fortran time = ', time, ' ms'! , ' end = ',end, ' beginning = ',beginning ! , ' rate = ',rate
     sum_squares = sum_squares + time*time
     sum_time = sum_time + time
  end do

  ! do k = 1, nbTimes
  !    call system_clock(beginning, rate)
  !    do ii=1,nbPoints,nb ! < nb is blocking factor
  !       do jj=1,nbPoints,nb
  !          do i=ii,min(nbPoints,ii+nb-1)
  !             do j=jj,min(nbPoints,jj+nb-1)
  !                distance(i,j) =   &
  !                     (points(1, i) - points(1, j))**2 + &
  !                     (points(2, i) - points(2, j))**2 + &
  !                     (points(3, i) - points(3, j))**2  
  !             end do
  !          end do
  !       end do
  !    end do
  !    call system_clock(end)
  !    time = 1e3*real(end - beginning) / real(rate)
  !    sum_squares = sum_squares + time*time
  !    sum_time = sum_time + time
  ! end do


  
  write (*,8) sum_time/nbTimes, sqrt( sum_squares/nbTimes - sum_time/nbTimes * sum_time/nbTimes)

8 format(' Time taken: ', F6.2, ' +- ', F0.2, ' ms')

  
end program main
