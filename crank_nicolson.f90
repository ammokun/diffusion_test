program crank_nicolson

use module_param
implicit none
    
!integer, parameter :: nmax = 1000

real(8), parameter :: cfl = 0.01d0
real(8), parameter :: x0 = 0.0d0
real(8), parameter :: x1 = 1.0d0
!real(8), parameter :: dx = 0.01

real(8), dimension(0:imax) ::q,x,q1,source_x
real(8), dimension(0:imax) ::q0
real(8) :: dt
real(8), dimension(0:imax) :: D_eff



do j=0,imax
    x(j)=dx * j
    D_eff(j)=-1.0d0
    source_x(j)=0.0d0
end do

do j=1,imax
    q0(j)=x(j)*1.0d2
end do


q=q0
q(0)=0.0d0
q(imax)=0.0d0
do n=1,1000

    call diffudion_cn(D_eff,q1,q,source_x)
    q=q1
end do

open(10, file = 'output.dat', form = 'formatted')
do j = 0, imax
    write(10,*) x(j), q0(j), q(j)
end do
close(10)

!call write_vector(q,n_max)

end program

subroutine diffudion_cn(dif_x,nxa,nxb,source_x)
    use module_param
    implicit none

    real(8),intent(in) :: dif_x(0:imax)
    real(8),intent(in) :: nxb(0:imax)
    real(8),intent(in) :: source_x(0:imax)
    real(8),intent(out) :: nxa(0:imax)

    real(8), dimension(1:imax) :: a,b,c
    real(8), dimension(0:imax) ::rhs
    real(8) :: nu

    nu=dt_p/(2*dx*dx)

    do j=1,imax-1
        a(j)=dif_x(j)*nu
        b(j)=1.0d0-2.0d0*dif_x(j)*nu
        c(j)=dif_x(j)*nu
    end do

    do j=1,imax-1
        rhs(j)=-dif_x(j)*nu*nxb(j-1)+(1.0d0+2.0d0*dif_x(j)*nu)*nxb(j)-dif_x(j)*nu*nxb(j+1)
    end do

    call triv_solver(a,b,c,imax-1,rhs)
    
    nxa=rhs
    nxa(0)=0.0d0
    nxa(imax)=0.0d0

end subroutine

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