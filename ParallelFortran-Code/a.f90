program a
  implicit none
  real*4 rand
  integer k,j,b(2,3), c(3)
 
  do j=1,3
     do k=1,2
        b(j,k)=floor(10*rand(0))
     end do
     c(j)=floor(10*rand(0))
  end do
  write(*,*) 'b',b
  write(*,*) 'c',c
  write(*,*) 'sum(b)', sum(b)
  write(*,*) spread(c,1,2)

  write(*,*) 3.5**2, 3.5**2.0 
end program a
