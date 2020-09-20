program tulli12
use mod_global
implicit none
integer :: i,j,k,total_trajectories,number_of_cores,ntraj,TT,t,loop
real*8 :: x1,dx1,start,finish
real*8, dimension(:),allocatable :: population_matrix
real*8, dimension(:),allocatable :: Time_matrix
!call setup_initial_values1
open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)
allocate(Time_matrix(TT))
allocate(population_matrix(TT))
do j=1,TT
population_matrix(TT)=0
enddo
do i=1,ntraj
call setup_initial_values2


!dx1=0.04
!pos(1)=20.0
!call potential
!call nonadiabaticvector
!call Rungekutta
!call gs
!write(*,*) sum(c(12,1:Ne)*conjg(c(12,1:Ne)))
!write(*,*) sum(c(1:Ne,5)*conjg(c(1:Ne,5)))
!write(*,*) matmul(Energy_hamil,transpose(Energy_hamil))
!write(*,*) (-2*real(A(1,7)*sum(v*(sum(Energy_hamil(:,1)*(matmul(Gradient_hamil(:,:,1),Energy_hamil(:,7))))/(Energy(7)-Energy(1))))))/A(1,1)
!write(*,*) sum(v*acw(1,7,:))
!write(*,*)sum(v*(sum(Energy_hamil(:,1)*matmul(Gradient_hamil(:,:,1),Energy_hamil(:,7)))/(Energy(7)-Energy(1))))
!write(*,*) sum(v*acw(1,7,:)
!write(*,*)sum(Energy_hamil(:,1)*matmul(Gradient_hamil(:,:,1),Energy_hamil(:,7)))/(Energy(7)-Energy(1))
!write(*,*) acw(1,7,:)
!write(*,*) dt*b(21,1)/A(1,1)

!do j=1,Ne
!write(*,*) S1(j,j)
!enddo
!call populations
!do j=1,Ne
!write(*,*)population
!enddo
!call populations
!write(*,*) 
!stop

call CPU_TIME(start)
call classical_evolution 


do k=1,TT
population_matrix(k)=population_matrix(k)+population_mat(k)
!population_mat(k)=0
enddo

!write(53,*) population_matrix/ntraj
!write(73,*) population_matrix(:,2)
call CPU_time(finish)
write(119,*)finish-start
enddo
!pos(1)=-10.0

!j=0
!do while (j<1000)
!pos(1)=pos(1)+dx1
!call potential

!j=j+1
!write(126,'(11ES15.5)')pos(1),sum(Energy(1:10))+U0
!write(35,*) pos(1),U0,H(1,1)
!enddo
!do j=1,Ne
 !k=j+1
!write(*,*) sum(c(j,:)*conjg(c(j,:)))
!enddo

!end do
do t=1,TT
Time_matrix(t)=10*t
enddo

write(77,*)   0.0000000000,     1.000000000
do loop=1,TT
write(77,*) Time_matrix(loop),1-population_matrix(loop)/ntraj
enddo

end program
!..............................................................................


