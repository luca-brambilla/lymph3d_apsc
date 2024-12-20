
!> @brief Setup for the forcing term, boundary data and exact solution. Additionally, it contains a
!> subroutine called set_properties, responsible for configuring the properties of the problem
!> (penalization coefficient, IP method, reaction coefficient).
module problem_data_and_properties

    use global_parameters

    implicit none

    contains

    !> @brief time function to be multiplied to space function to have h(x,t)=f(x)*g(t) - nodal values
    function time_function(time)result(r)

        implicit none
        real(kind=8) :: r
        real(kind=8) :: time


        r = sin(SQRT2*PI*time)

    end function time_function

    !> @brief Forcing term - nodal values
    function f(lambda, mu, p)result(r)

        implicit none

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p
        real(kind=8) :: alpha, theta, c
        real(kind=8) :: lambda, mu

        call set_properties(alpha, theta, c)

        ! r(1) = 0
        ! r(2) = 0
        ! r(3) = - 9.8 * 2400
        r(1) = 3*PI**2*cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(lambda + 2*mu)
        r(2) = 3*PI**2*cos(PI*p(2))*sin(PI*p(1))*sin(PI*p(3))*(lambda + 2*mu)
        r(3) = 3*PI**2*cos(PI*p(3))*sin(PI*p(1))*sin(PI*p(2))*(lambda + 2*mu)


    end function f

    !> @brief Forcing term depending on time takes correct density, otherwise 0 - nodal values
    function f_time(lambda, mu, p, rho)result(r)

        implicit none

        real(kind=8) :: rho

        real(kind=8), dimension(3) :: r !> result vector
        real(kind=8), dimension(3) :: p !> physical coordinates vector
        real(kind=8) :: alpha, theta, c !>
        real(kind=8) :: lambda, mu      !> elastic parameters

        call set_properties(alpha, theta, c)

        ! r(1) = 0
        ! r(2) = 0
        ! r(3) = - 9.8 * 2400

        r(1) = (-2.0*rho + 3.0*(lambda+2.0*mu))*PI**2 * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
        r(2) = (-2.0*rho + 3.0*(lambda+2.0*mu))*PI**2 * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
        r(3) = (-2.0*rho + 3.0*(lambda+2.0*mu))*PI**2 * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))


    end function f_time

    function f_null(lambda, mu, p, rho)result(r)

        implicit none

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p
        real(kind=8) :: lambda, mu, rho

        r = 0.0d0

    end function f_null

    function gd_null(p,space_fun_tag)result(r)

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r = 0.0d0

    end function gd_null

    function gn_null(lambda,mu,normal,p,space_fun_tag)result(r)
        implicit none

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p, normal
        real(kind=8) :: lambda, mu

        r = 0.0d0

    end function gn_null

    function gd_stat(p,space_fun_tag)result(r)

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r = 2.0d0

    end function gd_stat

    function gn_stat(lambda,mu,normal,p,space_fun_tag)result(r)

        implicit none

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p, normal
        real(kind=8) :: lambda, mu

        r = 0.0d0

    end function gn_stat

    !> @brief Dirichlet boundary data - nodal values
    function gd(p,space_fun_tag)result(r)

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        ! r = 0.0

        if(space_fun_tag==1) then
            r(1) = cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif
        if(space_fun_tag==2) then
            r(1) = cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif
        if(space_fun_tag==3) then
            r(1) = cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif

    end function gd

    !> @brief Dirichlet boundary data depending on time - nodal values
    function gd_time(p,space_fun_tag,time)result(r)

        real(kind=8) :: time

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        ! r = 0.0

        if(space_fun_tag==1) then
            r(1) = sin(SQRT2*PI*time) * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(SQRT2*PI*time) * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(SQRT2*PI*time) * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif
        if(space_fun_tag==2) then
            r(1) = sin(SQRT2*PI*time) * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(SQRT2*PI*time) * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(SQRT2*PI*time) * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif
        if(space_fun_tag==3) then
            r(1) = sin(SQRT2*PI*time) * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
            r(2) = sin(SQRT2*PI*time) * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            r(3) = sin(SQRT2*PI*time) * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))
        endif

    end function gd_time

    !> @brief Surface traction data - nodal values
    function gn(lambda,mu,normal,p,space_fun_tag)result(r)

        implicit none

        integer(kind=4) :: space_fun_tag
        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p, normal
        real(kind=8), dimension(3,3) :: stress_tensor
        real(kind=8) :: lambda, mu
        integer(kind=4) :: i, j

        ! r(1) = 0
        ! r(2) = 0
        ! r(3) = - 23520

        if(space_fun_tag==1) then
            stress_tensor(1,1) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(3))*sin(PI*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*PI*cos(PI*p(2))*cos(PI*p(3))*sin(PI*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif
        if(space_fun_tag==2) then
            stress_tensor(1,1) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(3))*sin(PI*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*PI*cos(PI*p(2))*cos(PI*p(3))*sin(PI*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif
        if(space_fun_tag==3) then
            stress_tensor(1,1) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(2,2) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(3,3) = -PI*sin(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))*(3*lambda + 2*mu)
            stress_tensor(1,2) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
            stress_tensor(2,1) = stress_tensor(1,2)
            stress_tensor(1,3) = 2*mu*PI*cos(PI*p(1))*cos(PI*p(3))*sin(PI*p(2))
            stress_tensor(3,1) = stress_tensor(1,3)
            stress_tensor(2,3) = 2*mu*PI*cos(PI*p(2))*cos(PI*p(3))*sin(PI*p(1))
            stress_tensor(3,2) = stress_tensor(2,3)
        endif

        r = 0.0
        do i=1,3
            do j=1,3
                r(i) = r(i) + stress_tensor(i,j) * normal(j)
            enddo
        enddo

    end function gn

    !> @brief Analytical solution - nodal values
    function uex(p)result(r)

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r(1) = cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
        r(2) = sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
        r(3) = sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))

    end function uex

        !> @brief Analytical solution - nodal values
    function uex_stat(p)result(r)

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r = 2.0d0

    end function uex_stat

    !> @brief Analytical solution depenting on time - nodal values
    function uex_time(p, time)result(r)

        real(kind=8) :: time

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r(1) = 1.0d0 + sin(SQRT2*PI*time) * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
        r(2) = 1.0d0 + sin(SQRT2*PI*time) * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
        r(3) = 1.0d0 + sin(SQRT2*PI*time) * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))

    end function uex_time

    !> @brief Initial condition for displacement - nodal values
    function ic_displacement(p)result(r)

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        !r = 0.0d0
        r = 2.0d0 + 0.0d0*p

    end function ic_displacement

    !> @brief Initial condition for velocity - nodal values
    function ic_velocity(p)result(r)

        real(kind=8), dimension(3) :: r
        real(kind=8), dimension(3) :: p

        r(1) = SQRT2*PI * cos(PI*p(1))*sin(PI*p(2))*sin(PI*p(3))
        r(2) = SQRT2*PI * sin(PI*p(1))*cos(PI*p(2))*sin(PI*p(3))
        r(3) = SQRT2*PI * sin(PI*p(1))*sin(PI*p(2))*cos(PI*p(3))

    end function ic_velocity

    !> @brief Set the properties of the numerical method
    subroutine set_properties(alpha, theta, c)

        implicit none

        real(kind=8) :: alpha, theta, c

        alpha = 10 ! penalty coefficient (which appears in the definition of the penalization function)
        theta = - 1 ! IP method (theta = -1 ---> SIP, theta = 0 ---> IIP, theta = 1 ---> NIP)
        c = 0 ! coefficient of the reaction term

    end subroutine set_properties

end module problem_data_and_properties