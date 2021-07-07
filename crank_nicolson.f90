program crank_nicolson

use module_param
implicit none
    
integer, parameter :: jmax=1000
!integer, parameter :: nmax = 1000

real(8), parameter :: cfl = 0.01d0
real(8), parameter :: x0 = 0.0d0
real(8), parameter :: x1 = 1.0d0
real(8), parameter :: dx = 0.001

real(8), parameter :: nu=0.5

real(8), dimension(1:jmax) :: a,b,c
real(8), dimension(0:jmax) ::rhs

real(8), dimension(0:jmax) ::q,x,qt
real(8), dimension(0:jmax) ::q0
real(8) :: dt
real(8), dimension(0:jmax) :: D_eff
dt =nu*dx*dx



do j=0,jmax
    x(j)=dx * j
    D_eff(j)=-1.0d0
end do

do j=1,jmax-1
    a(j)=D_eff(j)*nu
    b(j)=1.0d0-2.0d0*D_eff(j)*nu
    c(j)=D_eff(j)*nu
end do

do j=1,jmax
    !if(x(j)<0.50d0) then
    !    q0(j)=x(j)*2.0d0
    !else
    !    q0(j)=-x(j)*2.0d0+2.0d0
    !end if
    q0(j)=sin(pi*x(j))
end do


q=q0
q(0)=0.0d0
q(jmax)=0.0d0
do n=1,1000


    do j=1,jmax-1
        !rhs(j)=0.25d0*nu*q(j-1)+q(j)-0.25d0*nu*q(j+1)
        rhs(j)=-D_eff(j)*nu*q(j-1)+(1.0d0+2.0d0*D_eff(j)*nu)*q(j)-D_eff(j)*nu*q(j+1)
    end do

    !call triv_solver(a,b,c,jmax-1,rhs)

    q=rhs
    
    q(0)=0.0d0
    q(jmax)=0.0d0
end do

do j=1,jmax
    !if(x(j)<0.50d0) then
    !    q0(j)=x(j)*2.0d0
    !else
    !    q0(j)=-x(j)*2.0d0+2.0d0
    !end if
    qt(j)=sin(pi*x(j))*exp(-pi*pi*n*dt)
end do


write(*,*) "dx,dt: ",dx,dt
write(*,*) "t: ",n*dt

open(10, file = 'output.dat', form = 'formatted')
do j = 0, jmax
    write(10,*) x(j), q0(j), q(j),qt(j), sqrt((q(j)-qt(j))**2)
end do
close(10)

!call write_vector(q,n_max)

end program

subroutine triv_solver(a,b,c,n_max,xx)
    implicit none
!!!!--------------------
!!!! Thomas method
!!!!  |b1 c1         |
!!!!  |a2 b2 c2      |
!!!!A=|  a3 b3 c3    |
!!!!  |              |
!!!!  |        an bn |
!!!! Axx=rhs
    
    integer, intent(in) :: n_max
    integer :: n_up, n_back
    real(8),intent(in) :: a(n_max),b(n_max),c(n_max)
    real(8),intent(inout) :: xx(n_max)
    real(8) :: bi

    real(8), dimension(n_max) ::cc
    
    cc(1)=c(1)/b(1)
    xx(1)=xx(1)/b(1)
    
    do n_up=2,n_max
        bi=1.0d0/(b(n_up)-a(n_up)*cc(n_up-1))
        !write(*,*) bi
        cc(n_up)=c(n_up)*bi
        xx(n_up)=(xx(n_up)-a(n_up)*xx(n_up-1))*bi
    end do
    
    do n_back=n_max-1,1,-1
        xx(n_back)=xx(n_back)-xx(n_back+1)*cc(n_back)
        !write(*,*) "n_back, xx(n_back)=", n_back, xx(n_back)
    end do
    !write(*,*) "cc="
    !call write_vector(cc,n_max)
    !write(*,*) "xx="
    !call write_vector(xx,n_max)
    
end subroutine
    
subroutine write_matrix(ma,ndmax)
    
    real(8),intent(in) :: ma(1:ndmax,1:ndmax)
    !integer :: i,j

    do i=1,ndmax
        write(*,*) ma(i,1),ma(i,2),ma(i,3)
    end do
end subroutine
    
subroutine write_vector(va,ndmax)
    
    real(8),intent(in) :: va(1:ndmax)
    !integer :: i

    do i=1,ndmax
        write(*,*) va(i)
    end do
end subroutine