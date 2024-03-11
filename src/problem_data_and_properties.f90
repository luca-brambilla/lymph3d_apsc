module problem_data_and_properties

    implicit none

    contains

    ! Forcing term
    function f(lambda, mu, p)result(r)

        implicit none

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p
        real(kind=8) :: alpha, theta, c
        real(kind=8) :: lambda, mu
        real(kind=8), parameter :: pi = 4.d0*datan(1.d0) 

        call set_properties(alpha, theta, c)

        ! r(1) = 0
        ! r(2) = 0
        ! r(3) = - 9.8 * 2400

        r(1) = 3*pi**2*cos(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(lambda + 2*mu)
        r(2) = 3*pi**2*cos(pi*p(2))*sin(pi*p(1))*sin(pi*p(3))*(lambda + 2*mu)
        r(3) = 3*pi**2*cos(pi*p(3))*sin(pi*p(1))*sin(pi*p(2))*(lambda + 2*mu)

    end function f
       
    ! Dirichlet boundary data
    function gd(p,tag)result(r)
    
        integer(kind=4) :: tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p
        real(kind=8), parameter :: pi = 4.d0*datan(1.d0)

        ! r = 0.0

        if(tag==1) then
            r(1) = cos(pi*p(1))*sin(pi*p(2))*sin(pi*p(3)) 
            r(2) = sin(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            r(3) = sin(pi*p(1))*sin(pi*p(2))*cos(pi*p(3))
        endif
        if(tag==2) then
            r(1) = cos(pi*p(1))*sin(pi*p(2))*sin(pi*p(3)) 
            r(2) = sin(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            r(3) = sin(pi*p(1))*sin(pi*p(2))*cos(pi*p(3))
        endif
        if(tag==3) then
            r(1) = cos(pi*p(1))*sin(pi*p(2))*sin(pi*p(3)) 
            r(2) = sin(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            r(3) = sin(pi*p(1))*sin(pi*p(2))*cos(pi*p(3))
        endif

    end function gd

    ! Surface traction data
    function gn(lambda,mu,normal,p,tag)result(r)
    
        implicit none

        integer(kind=4) :: tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p, normal
        real(kind=8), dimension(3,3) :: stress_tensor
        real(kind=8) :: lambda, mu
        real(kind=8), parameter :: pi = 4.d0*datan(1.d0) 
        integer(kind=4) :: i, j

        ! r(1) = 0
        ! r(2) = 0
        ! r(3) = - 23520

        if(tag==1) then
            stress_tensor(1,1) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(3))*sin(pi*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*pi*cos(pi*p(2))*cos(pi*p(3))*sin(pi*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif
        if(tag==2) then
            stress_tensor(1,1) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(3))*sin(pi*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*pi*cos(pi*p(2))*cos(pi*p(3))*sin(pi*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif
        if(tag==3) then
            stress_tensor(1,1) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -pi*sin(pi*p(1))*sin(pi*p(2))*sin(pi*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*pi*cos(pi*p(1))*cos(pi*p(3))*sin(pi*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*pi*cos(pi*p(2))*cos(pi*p(3))*sin(pi*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif

        r = 0.0
        do i=1,3
            do j=1,3
                r(i) = r(i) + stress_tensor(i,j) * normal(j)
            enddo
        enddo

    end function gn
    
    ! Analytical solution
    function uex(p)result(r)

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p
        real(kind=8), parameter :: pi = 4.d0*datan(1.d0)

        r(1) = cos(pi*p(1))*sin(pi*p(2))*sin(pi*p(3)) 
        r(2) = sin(pi*p(1))*cos(pi*p(2))*sin(pi*p(3))
        r(3) = sin(pi*p(1))*sin(pi*p(2))*cos(pi*p(3))

    end function uex
    
    ! Set the properties of the numerical method
    subroutine set_properties(alpha, theta, c)

        implicit none

        real(kind=8) :: alpha, theta, c
        
        alpha = 10 ! penalty coefficient (which appears in the definition of the penalization function)
        theta = - 1 ! IP method (theta = -1 ---> SIP, theta = 0 ---> IIP, theta = 1 ---> NIP)
        c = 0 ! coefficient of the reaction term
        
    end subroutine set_properties
    
end module problem_data_and_properties