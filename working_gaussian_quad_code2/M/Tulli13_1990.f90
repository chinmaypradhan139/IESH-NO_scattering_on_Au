module mod_global

implicit none
  integer::nsize,INFO,total_dimensions,Hi,Ne,nold!(be careful with setting dt) 
  real*8 :: mass,tau,planks_constant,kT,Vr,Energy_of_diabat,U0,dt,dtc,total_time
  complex*16 :: population
  real*8 :: time,rnd,omega,gh,dG,Band_width,RAND
  complex*16,dimension(:,:), allocatable :: c,b,A
  real*8,dimension(:,:), allocatable :: Energy_hamil,g,Identity,H,Grd,vdij
  real*8,dimension(:,:,:), allocatable :: Gradient_hamil,acw
  real*8,dimension(:), allocatable :: Energy,acc1,E_levels,w_levels,knot_x
  real*8,dimension(:), allocatable :: pos,v,old_acc2,momentum,population_mat
  integer,dimension(:), allocatable :: lambda
  complex*16, dimension(:), allocatable :: cwork
  integer(kind=4) ::iseed
  !complex*16, dimension(20,20) ::  S1,S2,S3
  !complex*16, dimension(20,20) :: diagS1,diagS2,diagS3
  !complex*16, dimension(20) :: es1,es2,es3
  !real*8, dimension(40,40,20) :: rho_a
  

contains
!.......................................................................
subroutine setup_initial_values1
implicit none
integer(kind=4) :: O,j
integer(kind=4), dimension(:),allocatable::x
real*8 :: rhos,wt
integer:: ut,xt,ip,kn
real*8, dimension(14) :: inpot

 open(25,file='fort.23')
 do ip=1,14
 read(25,*) inpot(ip)
 enddo
 close(25)



Hi=inpot(1)
Ne=int(Hi/2)
tau=inpot(4)
Band_width=inpot(5)
mass=inpot(6)
omega=inpot(7)
gh=inpot(8)
dG=inpot(9)
KT=inpot(10)
dtc=inpot(11)
dt=inpot(12)
total_dimensions=inpot(13)
wt=inpot(14)
total_time=wt*5000
rhos=real(Hi)/Band_width
Vr=sqrt(tau/(2*3.1416))
  allocate(pos(total_dimensions))
  allocate(v(total_dimensions))
  allocate(momentum(total_dimensions))
  allocate(acc1(total_dimensions))
  allocate(old_acc2(total_dimensions))
  allocate(acw(Hi,Hi,total_dimensions))
  allocate(H(Hi,Hi))
  allocate(Energy_hamil(Hi,Hi)) 
  allocate(Energy(Hi))
  allocate(Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c(Ne,Hi))
  allocate(b(Ne,Hi))
  allocate(A(Ne,Hi))
  allocate(g(Ne,Hi))
  allocate(lambda(Ne))
  allocate(Identity(Hi,Hi))
  allocate(Grd(Hi,Hi))
  allocate(vdij(Hi,Hi))
  allocate(population_mat(int(total_time/dtc)))
  allocate(E_levels(int(Hi/2)))
  allocate(knot_x(int(Hi/2)))
  allocate(w_levels(int(Hi/2)))
call random_seed(size=O)
  allocate(x(O))
   do j=1,O
   x(j)=j**6*iseed+2777772
   enddo
  call random_seed(put=x)
 open(2,file='rndom',status='new')
 write(2,*) iseed,x,O
 


open(167,file='raw_x.txt')
do ut=1,int(Hi/2)
read(167,*)E_levels(ut)
enddo
close(167)
do kn=1,int(Hi/2)
knot_x(kn)=E_levels(int(Hi/2)-kn)
enddo


open(169,file='raw_w.txt')
do xt=1,int(Hi/2)
read(169,*)w_levels(xt)
enddo
close(169)







  

end subroutine
!.........................................................
  
subroutine gaussian_random_number(rnd0)
!USE IFPORT
   !! generates gaussian distribution with center 0, sigma 1
   !! q0+sig*rnd gives center=q0, sigma=sig
   implicit none
   integer(kind=4) :: n,j,M,O,k
   real*8,intent(out)::rnd0
   !integer(kind=4), dimension(:),allocatable::x
   !integer(kind=4), dimension(:),allocatable::y
   real*8 rnd1,rnd2,pi
   !call random_seed(size=M)
   !allocate(x(M))
   !do j=1,M
   !x(j)=j**2*iseed+1
   !enddo
    
   !call random_seed(size=O)
   !allocate(y(O))
   !do k=1,O
   !y(k)=3*iseed+21
   !enddo



   pi=dacos(-1.d0)
   !call random_seed(put=x)
   call random_number(rnd1)
! write(110,*) time,rnd1
   !call random_seed(put=y)

   call random_number(rnd2)
!write(110,*) time,rnd2
   !call random_number(rnd2)
   rnd0 = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!.............................................................
subroutine setup_initial_values2
integer :: n,m,p,y,i
real*8 :: rnd2,rnd1
!dt=1
!dtc=10  !dtc should be an integer multiple of dt
!total_time=50000
call gaussian_random_number(rnd2)
call gaussian_random_number(rnd1)
momentum(1)=sqrt(mass*KT)*rnd2
v(1)=(momentum(1)/mass)  
!v(1)=0.00
pos(1)=rnd1/(sqrt(mass*omega**2/KT))
time=0.00




c=0
do i=1,Ne
   c(i,i+1)=1
enddo 














do n=1,Ne
lambda(n)=n
enddo



!write(54,*) time,v(1)

call potential

c=matmul(c,Energy_hamil)
!do m=1,Ne
!do p=1,Hi
!write(101,"(f12.5$)") c(m,p)
!enddo
!write(101,*)
!enddo
!call populations
!write(77,*) time,1-real(population)

!do y=1,20
!write(66,*) pos(1),Energy(y)
!enddo


end subroutine 
!........................................................................
subroutine diag_wrapper(matrix,nsize,eigen_values,eigen_vectors)
real*8, intent(inout) :: matrix(nsize,nsize)
real*8, dimension(nsize,nsize) :: mat
integer LWORK,nsize
real*8, allocatable :: WORK(:)
real*8, intent(out) :: eigen_vectors(nsize,nsize),eigen_values(nsize)
mat=matrix
LWORK=3*nsize-1
allocate(WORK(LWORK))
call dsyev('V','U',nsize,mat,nsize,eigen_values,WORK,LWORK,INFO)
eigen_vectors=mat

end subroutine

!..........................................................................
subroutine logm(mat,log_mat,n)
   !! http://arxiv.org/pdf/1203.6151v4.pdf
   implicit none
   integer,intent(in):: n
   real*8,intent(in):: mat(n,n)
   real*8,intent(out):: log_mat(n,n)
   integer i
   complex*16 T(n,n),en(n),vect(n,n)
   complex*16 dd(n,n)

   call schur(mat,T,n,en,vect,nold,cwork)

   dd=0.d0
   do i=1,n
     dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
   enddo

   log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm


!...................................................................................
!subroutine schur(mat,n,eigen_value,eigen_vect,nold,cwork)

  !! Diaganalizing matrix using dsyevr. First m_values eigen values and
  !eigenvectors computed.

  !! The module's common variables should contain:



  !! Initialize nold=0 



  !! nold makes sure that everytime value of n changes, work and iwork are
  !re-allocated for optimal performance.

  !! mat is destroyed after use.



  !implicit none

  !integer,intent(in) :: n

  !integer,intent(inout) :: nold

  !complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)

  !complex*16,intent(in) :: mat(n,n)

  !complex*16 T(n,n)

  !complex*16,allocatable,intent(inout):: cwork(:)

  !real*8 rwork(n)

  !complex*16 mat_c(n,n)



  !integer lwork

  !logical:: select

  !logical bwork(n)

  !integer sdim,info,AllocateStatus



  !T=mat



  !info=0

  !sdim=0



  !if(nold.ne.n .or. .not.allocated(cwork)) then

  !if(nold.ne.n) then

    !lwork=-1

    !if(allocated(cwork))deallocate(cwork)

    !allocate(cwork(n))

    !call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)

    !lwork=int(cwork(1))

    !deallocate(cwork)

    !allocate(cwork(lwork),STAT=AllocateStatus)

    !if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"

    !nold=n

  !endif



  !lwork=size(cwork)

  !call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)

  !if(info.ne.0) then

    !write(6,*) "problem in diagonalization",info

    !stop

  !endif

!end subroutine schur
!.............................................................................

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
   !! Diaganalizing matrix using dsyevr. First m_values eigen values and
!eigenvectors computed.
   !! The module's common variables should contain:

   !! Initialize nold=0

   !! nold makes sure that everytime value of n changes, work and iwork
!are re-allocated for optimal performance.
   !! mat is destroyed after use.

   implicit none
   integer,intent(in) :: n
   integer,intent(inout) :: nold
   complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
   real*8,intent(in) :: mat(n,n)
   complex*16,intent(out) :: T(n,n)
   complex*16,allocatable,intent(inout):: cwork(:)
   real*8 rwork(n)
   complex*16 mat_c(n,n)

   integer lwork
   logical:: select
   logical bwork(n)
   integer sdim,info,AllocateStatus

   T=mat

   info=0
   sdim=0

   if(nold.ne.n .or. .not.allocated(cwork)) then
   !if(nold.ne.n) then
     lwork=-1
     if(allocated(cwork))deallocate(cwork)
     allocate(cwork(n))
     call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
     lwork=int(cwork(1))
     deallocate(cwork)
     allocate(cwork(lwork),STAT=AllocateStatus)
     if(allocatestatus.ne.0) write(6,*)"problem in schur, allocation"
     nold=n
   endif

   lwork=size(cwork)
   call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
   if(info.ne.0) then
     write(6,*) "problem in scur",info
     stop
   endif

end subroutine schur
!...................................................................................................................
subroutine potential
implicit none
real*8 :: U1
integer :: i,j,k,l



U0=0.5*mass*(omega**2)*(pos(1))**2
U1=0.5*mass*(omega**2)*(pos(1)-gh)**2+dG
H(1,1)=U1-U0
    do i=2,Hi-1
        if (i.le.int(Hi/2)) then
                H(1,i)=sqrt(Band_width*w_levels(i)/2)*Vr
                H(i,1)=sqrt(Band_width*w_levels(i)/2)*Vr
        else
                H(1,i)=sqrt(Band_width*w_levels(i-int(Hi/2))/2)*Vr
                H(i,1)=sqrt(Band_width*w_levels(i-int(Hi/2))/2)*Vr
        end if
        do j=2,Hi-1
                if (i.eq.j) then
                        if (i.le.int(Hi/2)) then
                                H(i,j)=-(Band_width/2)*(0.5+0.5*knot_x(i))
                        else
                                H(i,j)=(Band_width/2)*(0.5+0.5*knot_x(Hi-i+1))
                        end if
                else
                        H(i,j)=0.0
                end if
        end do
   enddo







nsize=Hi
call diag_wrapper(H,nsize,Energy,Energy_hamil)






Gradient_hamil(1,1,1)=-mass*(omega**2)*gh
do k=2,Hi
do l=2,Hi
Gradient_hamil(l,1,1)=0
Gradient_hamil(1,k,1)=0
if (l.eq.k) then
Gradient_hamil(l,k,1)=0
else
Gradient_hamil(l,k,1)=0
endif
enddo
enddo





























end subroutine
!........................................................................................
!subroutine potential
!implicit none
!real*8 :: U1
!integer :: i,j,k,l,r,s,x,y,u
!real*8, dimension(Hi,Hi) :: F
!real*8, dimension(Hi) :: U_zero
!real*8, dimension(Hi,Hi) :: G
!Energy_of_diabat=Band_width/real(Hi)
!U0=0.5*mass*(omega**2)*pos(1)**2
!U1=0.5*mass*(omega**2)*(pos(1)-gh)**2+dG
!H(1,1)=U1-U0
!do j=2,Hi
!do i=2,Hi
!H(i,1)=Vr
!H(1,j)=Vr
!if (i.eq.j) then
!H(i,j)=i*(Band_width/2)+U0
!H(i,j)=i*(Band_width/(Hi-2))-Band_width*(Hi+2)/(2*Hi-4)
!else
!H(i,j)=0
!endif
!enddo
!enddo
!do r=2,Hi
!H(r,1)=Vr
!H(1,r)=Vr
!enddo
!do x=1,Hi
!do y=1,Hi
!write(*,*) x,y,H(x,y)
!enddo
!enddo
!nsize=Hi
!call diag_wrapper(H,nsize,Energy,Energy_hamil)
!c=matmul(c,Energy_hamil)




!Gradient_hamil(1,1,1)=-mass*(omega**2)*gh
!do k=2,Hi
!do l=2,Hi
!Gradient_hamil(l,1,1)=0
!Gradient_hamil(1,k,1)=0
!if (l.eq.k) then
!Gradient_hamil(l,k,1)=0
!else
!Gradient_hamil(l,k,1)=0

!endif
!enddo
!enddo
!F=matmul(transpose(Energy_hamil),Gradient_hamil(:,:,1))
!Grd=matmul(F,Energy_hamil)

!end subroutine
!...........................................................................
subroutine nonadiabaticvector(t,u)
implicit none
integer :: acwi
integer, intent(in) :: t,u
do acwi=1,total_dimensions

if (t.ne.u) then
acw(t,u,acwi)=(sum(Energy_hamil(:,t)*matmul(Gradient_hamil(:,:,acwi),Energy_hamil(:,u))))/(Energy(t)-Energy(u))
else
acw(t,u,acwi)=0
endif

end do
end subroutine
!..........................................................................
subroutine force
implicit none
integer :: acc1i,z,u,q
real*8, dimension(total_dimensions,Ne) :: m_acc1
!call potential
do acc1i=1,total_dimensions

do z=1,Ne
m_acc1(acc1i,z)=-((sum((Energy_hamil(:,lambda(z)))*matmul(Gradient_hamil(:,:,acc1i),Energy_hamil(:,lambda(z)))))/mass)
end do
end do
acc1=sum(m_acc1(:,1:Ne))-omega**2*pos(1)
end subroutine
!..........................................................................
subroutine Rungekutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi 
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-sum(v*acw(p,q,:))
enddo
enddo

k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine
!.........................................................................

subroutine Rungefutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M
!real*8,intent(in), dimension(Hi,Hi) :: dij
do i=1,Hi
do j=1,Hi
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-vdij(p,q)
enddo
enddo

!k1=dt*matmul(c,M)
!k2=dt*matmul((c+k1/2),M)
!k3=dt*matmul((c+k2/2),M)
!k4=dt*matmul((c+k3),M) 


k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
!do i=1,Hi
!do j=1,Hi
!write(106,"(2ES10.2$)") M(i,j)
!enddo
!write(106,*)
!enddo
!stop
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine



!...................................................................................
subroutine gs
implicit none
integer :: i,d,r,s,l,m,o,t
!real*8, dimension(Ne,Hi) :: b,A
integer :: findloc
integer, dimension(Ne) :: p
complex*16 :: det_S1,det_S2
complex*16, dimension(Ne,Ne) ::  S1,S2
!real*8, dimension(Ne) :: es1,es2
nold=0
p=lambda
do i=1,Ne
do d=1,Hi
!!do i=1,Ne
!t=1
!do o=0,Ne
if (FINDLOC(p,d,1).eq.0) then

do r=1,Ne
S1(:,r)=c(:,p(r))
enddo
!call schur(S1,Ne,es1,diagS1,nold,cwork)
call modulus(S1,Ne,det_S1)

p(i)=d


do s=1,Ne
S2(:,s)=c(:,p(s))
enddo
!call schur(S2,Ne,es2,diagS2,nold,cwork)
call modulus(S2,Ne,det_S2)



!A(d,i)=product(es1)*product(es2)
A(d,i)=det_S1*conjg(det_S2)
b(d,i)=-2*real(A(d,i)*vdij(d,i))
!g(i,d)=dt*b(d,i)/A(i,i)
!A(i,i)=real(product(es1)*conjg(product(es1)))
A(i,i)=det_S1*conjg(det_S1)
g(i,d)=dt*real(b(d,i)/A(i,i))
!else if (i.eq.d) then
!do m=1,Ne
!S3(:,m)=c(:,p(m))
!end do
!call schur(S3,Ne,es3,diagS3,nold,cwork)

!A(i,d)=product(es3)*product(es3)
else
A(i,d)=0
b(i,d)=0
g(i,d)=0
endif

enddo
p(i)=lambda(i)
enddo
end subroutine
!...........................................................................
subroutine quantum_evolution

implicit none
integer :: k,j,r,s,i,t,u,M
real*8 :: EE,TE,KE,rnd
integer(kind=4), dimension(:), allocatable :: z
!real*8, dimension(Hi,Hi) :: old_Energy_hamil


!call random_seed(size=M)

!allocate(z(M))

!do k=1,s
!z(k)=k*iseed+2
!end do                  
!call random_seed(put=z)
call random_number(rnd)
!do while(time<50000)
!time=time+dt
!call potential
!old_Energy_hamil=Energy_hamil
!call force
!do t=1,Hi
!do u=1,Hi
!if (t.ne.u) then
!call nonadiabaticvector(t,u)
!else
!acw(t,u,:)=0
!end if  
!end do
!end do

!call rungekutta
!call gs
!call hopping_probablity
!call velocity_verlet2
!call populations




!KE=0.5*mass*v(1)**2
!EE=Energy(lambda(1))
!do i=2,Ne
!EE=Energy(lambda(i))+EE
!enddo


!TE=EE+KE










!itime=time+dt
!call potential
!call vdotd(old_Energy_hamil)
call Rungefutta
call gs
jloop : do t=1,Ne
do u=1,Hi
if (g(t,u)>rnd) then
call nonadiabaticvector(t,u) 
!write(104,*)TE
call hop(t,u)
!write(104,*)TE,"hop"
!write(103,*) time,t,u
exit jloop
endif
enddo
enddo jloop
!call populations
!do r=1,Ne-1
!is=r+1
!write(18,*)
!pos(1),sum(Energy(1:Ne))+0.5*mass*v(1)**2,
!enddo



!do j=1,Ne
 !k=j+1
!write(*,*) sum(c(j,:)*conjg(c(j,:)))
!enddo
!KE=0.5*mass*sum(v(1:total_dimensions)**2)
!EE=Energy(lambda(1))
!do i=2,Ne
!EE=Energy(lambda(i))+EE
!enddo


!TE=EE+KE

!write(80,*)pos(1),sum((c(j,:))*conjg(c(j,:)))
!write(42,'(11ES15.5)')pos(1),Energy(1),Energy(2)
!write(36,'(11ES15.5)') g(1,21:40)


!write(98,*)pos(1),TE
!write(81,*)pos(1),time





!do r=1,20
!do s=21,40
!if (g(r,s)>rnd) then
!write(21,*)g(r,s),'yeah'
!endif
!enddo
!enddo
!enddo
end subroutine
!.........................................................................
subroutine classical_evolution
integer :: n,p,r,TT,yt,i
real*8, dimension(Hi,Hi) :: old_Energy_hamil
real*8,dimension(Hi) :: signature
real*8 :: KE,TE,EE
call potential
call force



TT=int(total_time/dtc)

do while(time.le.total_time)

old_Energy_hamil=Energy_hamil
call velocity_verlet

!call vdotd(old_Energy_hamil,signature)
call signt(old_Energy_hamil,signature)
do r=1,Hi
 if (signature(r)<0) then
    Energy_hamil(:,r)=-Energy_hamil(:,r)
 endif
enddo

call vdotd(old_Energy_hamil)

n=int(dtc/dt)
do while(n>0) 
call quantum_evolution

n=n-1
enddo



call populations
yt=int(time/dtc)
population_mat(yt)=real(population)


!write(102,*) time,1-real(population)
!do p=1,Ne
!write(91,*)pos(1),real(sum(c(p,:)*conjg(c(p,:))))
!enddo

!write(34,'(11ES15.5)')pos(1),Energy(1:10)
!write(54,*) time,v(1)
!WRITE(91,*) 
!write(63,*) time,pos(1)
!write(45,"(f15.7,10i5)") time,lambda
!write(34,*)pos(1),sum(c(:,p)*conjg(c(:,p)))
!time=time+dtc
!call force

!KE=0.5*mass*(v(1)**2)
!EE=Energy(lambda(1))
!do i=2,Ne
!EE=Energy(lambda(i))+EE
!enddo
!TE=KE+EE
!write(72,*) time,TE
!write(73,*) time,EE
!write(74,*) time,KE


!write(42,'(11ES15.5)')pos(1),EE+U0

time=time+dtc


!






enddo




end subroutine
!................................................................................
subroutine hop(t,u)
implicit none
integer, intent(in) :: t,u
real*8 :: para_v
real*8, dimension(:),allocatable :: perp_v
!call force
!call nonadiabaticvector
!call Rungekutta
!call gs

allocate(perp_v(total_dimensions))
!do t=1,Hi
!do i=1,Ne
if (((sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))) & 
**2+(2*Energy(lambda(t))/mass)-(2*Energy(u)/mass))>0) then
    para_v=sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))
    perp_v=v-para_v*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))
    para_v=((para_v)/abs(para_v))*sqrt(para_v**2+(2*Energy(lambda(t))/mass)-(2*Energy(u)/mass))
    v=(para_v)*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))+perp_v
    lambda(t)=u
end if
!end do
!end do
!Note: (sum(v*acw(lambda(t),u,:))/norm(acw(lambda(t),u,:))) is the parallel component of velocity in the direction of nonadiabatic coupling vector
end subroutine
!................................................................................
subroutine velocity_verlet
real*8 :: delr(total_dimensions),delv(total_dimensions)
real*8 :: gama_dt,gamma_B,c0,c1,c2
gamma_B=2*omega
gama_dt=gamma_B*dtc
c0=dexp(-gama_dt)
c1=1.d0/gama_dt*(1.d0-c0)
c2=1.d0/gama_dt*(1.d0-c1)
call stochastic_force(delr,delv)
pos=pos+c1*dtc*v+c2*dtc*dtc*acc1+delr
old_acc2=acc1
call potential
call force
v=c0*v+(c1-c2)*dtc*old_acc2+c2*dtc*acc1+delv

end subroutine
!................................................................................
subroutine vdotd(old_Energy_hamil)
real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave
!real*8, intent(out),dimension(Hi,Hi) :: dij
real*8, dimension(Hi,Hi) :: Ut,logarithm_Ut
integer :: i

Ut=matmul(transpose(old_Energy_hamil),Energy_hamil)
call logm(Ut,logarithm_Ut,Hi)
vdij=logarithm_Ut/dtc




end subroutine
!........................................................................................
subroutine signt(old_Energy_hamil,signature)

real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave

integer :: i
real*8, intent(out),dimension(Hi) :: signature
do i=1,Hi
signature(i)=sum(old_Energy_hamil(:,i)*Energy_hamil(:,i))
enddo



end subroutine


















!....................................................................................
subroutine velocity_verlet2
pos=pos+v*dtc+0.5*acc1*dtc*dtc
old_acc2=acc1
call potential
call force
v=v+0.5*(acc1+old_acc2)*dtc
end subroutine

!...................................................................................
subroutine stochastic_force(delr,delv)
real*8, intent(out) :: delr(total_dimensions),delv(total_dimensions)
integer :: i
real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv,gdt,gamma_B
gamma_B=2*omega
gdt=gamma_B*dtc

do i=1,total_dimensions

sig_r=dtc*dsqrt(KT/mass*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
sig_v=dsqrt(KT/mass*(1-dexp(-2*gdt)))
sig_rv=(dtc*KT/mass*1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v) !correlation coefficient

call gaussian_random_number(rnd1)
call gaussian_random_number(rnd2)
delr(i)=sig_r*rnd1
delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
enddo

end subroutine stochastic_force
!....................................................................................
subroutine populations
integer :: a,i,j,k
real*8, dimension(Hi,Hi,Ne) :: rho_a
real*8, dimension(Hi,Hi) :: rho_d
real*8, dimension(Ne) :: LUMO_population
do a=1,Ne
do i=1,Hi
do j=1,Hi
!if (i.ne.j) then
rho_a(i,j,a)=c(a,i)*conjg(c(a,j))
!else if (i.eq.lambda(a)) then
!rho_a(i,i,a)=1
!else
!rho_a(i,i,a)=0
!end if
enddo
enddo
enddo

do k=1,Ne
rho_d=matmul(Energy_hamil,(matmul(rho_a(:,:,k),transpose(Energy_hamil))))
LUMO_population(k)=rho_d(1,1)
enddo


population=sum(LUMO_population(1:Ne))
!write(15,*)time,population!,lambda
end subroutine
!...................................................................................................................

subroutine modulus(matrix,n,determinant)
 IMPLICIT NONE
     complex*16, DIMENSION(n,n) :: matrix
     INTEGER, INTENT(IN) :: n
     complex*16 :: m, temp
     INTEGER :: i, j, k, l
     LOGICAL :: DetExists = .TRUE.
     complex*16,intent(out) :: determinant
     l = 1
     !Convert to upper triangular form
     
     DO k = 1, n-1
         IF (matrix(k,k) == 0) THEN
             DetExists = .FALSE.
             DO i = k+1, n
                 IF (matrix(i,k) /= 0) THEN
                     DO j = 1, n
                         temp = matrix(i,j)
                         matrix(i,j)= matrix(k,j)
                         matrix(k,j) = temp
                     END DO
                     DetExists = .TRUE.
                     l=-l
                     EXIT
                 ENDIF
             END DO
             IF (DetExists .EQV. .FALSE.) THEN
                 determinant= 0
                 return
             END IF
         ENDIF
         DO j = k+1, n
             m = matrix(j,k)/matrix(k,k)
             DO i = k+1, n
                 matrix(j,i) = matrix(j,i) - m*matrix(k,i)
             END DO
         END DO
     END DO

     !Calculate determinant by finding product of diagonal elements
     determinant= l
     DO i = 1, n
         determinant= determinant* matrix(i,i)
     END DO

END subroutine modulus
!...........................................................................................................
Subroutine Averaging









































end subroutine
!......................................................................................
end module














