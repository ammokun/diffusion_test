program adi_2d

use module_param
implicit none
    
integer, parameter :: jmax=101
integer, parameter :: kmax=101
!integer, parameter :: nmax = 1000

real(8), parameter :: cfl = 0.01d0
real(8), parameter :: x0 = 0.0d0
real(8), parameter :: x1 = 1.0d0
real(8), parameter :: dx = 0.01

real(8), parameter :: nu=0.1

real(8), dimension(1:jmax) :: a,b,c,rhs

real(8), dimension(0:jmax) ::x,y
real(8), dimension(0:jmax,0:jmax) ::q0,q
real(8) :: dt
real(8), dimension(0:jmax,0:jmax) :: D_eff
dt =nu*2*dx*dx

do j=0,jmax
    x(j)=dx*j
end do

do k=0,kmax
    y(k)=dx*k
end do

do k=0,kmax   
    do j=0,jmax
        D_eff(j,k)=1.0d0
        !q0(j,k)=exp((-(x(j)-0.5d0)**2-(y(k)-0.5d0)**2)/(2*0.10**2))
        q0(j,k)=0.0d0
    end do
end do


q=q0

do k=0,kmax
    q(0,k)=0.0d0
    q(jmax,k)=0.0d0
end do

do j=0,jmax
    q(j,0)=0.0d0
    q(j,kmax)=0.0d0
end do

do n=1,10000

    do m=1,jmax-1
        k=m
        do j=1,jmax
            a(j)=-0.50d0*D_eff(j,k)*nu
            b(j)=1.0d0+D_eff(j,k)*nu
            c(j)=-0.50d0*D_eff(j,k)*nu
        end do

        do j=1,jmax-1
            !rhs(j)=0.25d0*nu*q(j-1)+q(j)-0.25d0*nu*q(j+1)
            rhs(j)=0.50d0*nu*D_eff(j,k)*q(j,k-1)+(1.0d0-D_eff(j,k)*nu)*q(j,k)+0.50d0*nu*D_eff(j,k)*q(j,k+1)
        end do

        call triv_solver(a,b,c,jmax-1,rhs)

        do j=1,jmax-1
            q(j,k)=rhs(j)
        end do

        j=m
        do k=1,kmax
            a(k)=-0.50d0*D_eff(j,k)*nu
            b(k)=1.0d0+D_eff(j,k)*nu
            c(k)=-0.50d0*D_eff(j,k)*nu
        end do

        do k=1,kmax-1
            !rhs(j)=0.25d0*nu*q(j-1)+q(j)-0.25d0*nu*q(j+1)
            rhs(k)=0.50d0*nu*D_eff(j,k)*q(j-1,k)+(1.0d0-D_eff(j,k)*nu)*q(j,k)+0.50d0*nu*D_eff(j,k)*q(j+1,k)
        end do

        call triv_solver(a,b,c,kmax-1,rhs)

        do k=1,kmax-1
            q(j,k)=rhs(k)
        end do 
        

        do k=0,kmax
            q(0,k)=1.0d0
            q(jmax,k)=0.0d0
        end do
        
        do j=0,jmax
            q(j,0)=0.0d0
            q(j,kmax)=0.0d0
        end do
        
    end do
end do

write(*,*) "dx,dt: ",dx,dt
write(*,*) "t: ",n*dt

open(10, file = 'output.dat', form = 'formatted')
do k = 0, kmax
    do j=0,jmax
    write(10,*) x(j),y(k), q0(j,k), q(j,k)
    end do
    write(10,*) " "
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