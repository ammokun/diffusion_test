program crank_nicolson

use module_param
implicit none
    
integer, parameter :: jmax=101
!integer, parameter :: nmax = 1000

real(8), parameter :: cfl = 0.01d0
real(8), parameter :: x0 = 0.0d0
real(8), parameter :: x1 = 1.0d0
!real(8) :: dx, dt

integer, parameter :: n_max=3
real(8), dimension(n_max) :: a,b,c
real(8), dimension(n_max) ::rhs


do i=1,n_max
    rhs(i)=1.0d0
    a(i)=2.0d0*i
    b(i)=1.0d0*i
    c(i)=2.0d0*i
end do


!call write_vector(rhs,n_max)

call triv_solver(a,b,c,n_max,rhs)
!!!!--------------------
!!!! Thomas method
!!!!  |b1 c1         |
!!!!  |a2 b2 c2      |
!!!!A=|  a3 b3 c3    |
!!!!  |              |
!!!!  |        an bn |
!!!! Axx=rhs
call write_vector(rhs,n_max)
end program

subroutine triv_solver(a,b,c,n_max,xx)
    implicit none
    
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