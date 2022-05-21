!************************************************************************************
!>
!  Test of the integration routine.
!  Integrates various functions and prints results.

    program main

    use quadrature_module, wp => quadrature_wp

    implicit none

    real(wp),parameter :: zero      = 0.0_wp
    real(wp),parameter :: one       = 1.0_wp
    real(wp),parameter :: two       = 2.0_wp
    real(wp),parameter :: three     = 3.0_wp
    real(wp),parameter :: pi        = acos(-one)
    real(wp),parameter :: tol       = 1.0e-12_wp  !! error tolerance
    !real(wp),parameter :: tol       = 100*epsilon(one)

    type,extends(integration_class_1d) :: sin_type
        integer  :: ifunc   = 0     !! which function to use
        real(wp) :: amp     = zero  !! amplitude
        real(wp) :: freq    = zero  !! frequency
        real(wp) :: phase   = zero  !! phase
        integer  :: n_evals = 0     !! number of function evaluations
    end type sin_type

    type,extends(integration_class_2d) :: my_doub
        integer :: ifunc    = 0   !! which function to use
        integer :: n_evals  = 0   !! number of function evaluations
    end type my_doub
    type,extends(integration_class_3d) :: my_trip
        integer :: ifunc    = 0   !! which function to use
        integer :: n_evals  = 0   !! number of function evaluations
    end type my_trip
    type,extends(integration_class_4d) :: my_4
        integer :: ifunc    = 0   !! which function to use
        integer :: n_evals  = 0   !! number of function evaluations
    end type my_4
    type,extends(integration_class_5d) :: my_5
        integer :: ifunc    = 0   !! which function to use
        integer :: n_evals  = 0   !! number of function evaluations
    end type my_5
    type,extends(integration_class_6d) :: my_6
        integer :: ifunc    = 0   !! which function to use
        integer :: n_evals  = 0   !! number of function evaluations
    end type my_6

    type(sin_type) :: my_int
    type(my_doub)  :: doub
    type(my_trip)  :: trip
    type(my_4)  :: int_4d
    type(my_5)  :: int_5d
    type(my_6)  :: int_6d
    real(wp) :: a, b, ans, err, answer
    real(wp) :: xl, xu, yl, yu, zl, zu, ql, qu, rl, ru, sl, su
    integer  :: ierr    !! error code
    integer  :: i       !! counter
    integer  :: meth    !! method number
    integer  :: itest   !! test number for output printing

    itest = 0

    ! csv header:
    write(*,'(A)') 'Test number, Dimension, Method, tol, ans, ierr, err, Function evaluations, Actual Error'

    do i=1,size(set_of_quadrature_methods)  !test all the methods

        meth = set_of_quadrature_methods(i)%n_points

        !============================================
        ! single integral tests
        !============================================

        !set the coefficients (a simple sin function):
        my_int%amp      = one
        my_int%freq     = one
        my_int%phase    = zero
        my_int%ifunc    = 1

        ! Case 1:
        a = zero      !lower bound
        b = pi        !upper bound
        call run_test()

        ! Case 2:
        a = zero      !lower bound
        b = pi/two    !upper bound
        call run_test()

        ! Case 3:
        a = -one      !lower bound
        b = one       !upper bound
        call run_test()

        ! case 4:
        a = -2.34567_wp        !lower bound
        b = 12.34567_wp        !upper bound
        my_int%amp   = 0.092834_wp
        my_int%freq  = 19.87_wp
        my_int%phase = 77.0001_wp
        call run_test()

        ! case 4:
        a = -2.34567_wp    !lower bound
        b = 1.34567_wp     !upper bound
        my_int%amp   = 0.092834_wp
        my_int%freq  = 1.87_wp
        my_int%phase = 7.0001_wp
        call run_test()

        !other functions:
        my_int%ifunc = 2
        a = one  !lower bound
        b = two  !upper bound
        call run_test()

        my_int%ifunc = 3
        a = -20.0_wp  !lower bound
        b = -1.0_wp   !upper bound
        call run_test()

        my_int%ifunc = 4
        a = -pi  !lower bound
        b = pi   !upper bound
        call run_test()

        my_int%ifunc = 5
        a = 0.0_wp
        b = 10.1_wp
        call run_test()

        !============================================
        ! double integral tests
        !============================================

        doub%ifunc = 1
        xl = 0.3827812_wp
        xu = 1.02928_wp
        yl = 22.91281_wp
        yu = -111.928_wp
        call run_2d_test()

        doub%ifunc = 1
        xl = 0.0_wp
        xu = 1.0_wp
        yl = 0.1_wp
        yu = 1.0_wp
        call run_2d_test()

        doub%ifunc = 2
        xl = 0.0_wp
        xu = pi
        yl = 0.0_wp
        yu = 2.0_wp*pi
        call run_2d_test()

        doub%ifunc = 2
        xl = 0.0_wp
        xu = 1.0_wp
        yl = 0.1_wp
        yu = 1.0_wp
        call run_2d_test()

        !============================================
        ! triple integral tests
        !============================================

        trip%ifunc = 1
        xl = 2.0_wp
        xu = 3.0_wp
        yl = 1.0_wp
        yu = 2.0_wp
        zl = 0.0_wp
        zu = 1.0_wp
        call run_3d_test()

        !============================================
        ! 4d integral tests
        !============================================

        int_4d%ifunc = 1
        xl = 0.0_wp
        xu = 1.0_wp
        yl = 0.0_wp
        yu = 1.0_wp
        zl = 0.0_wp
        zu = 1.0_wp
        ql = 0.0_wp
        qu = 1.0_wp
        call run_4d_test()

        !============================================
        ! 5d integral tests
        !============================================

        int_5d%ifunc = 1
        xl = 0.0_wp
        xu = 1.0_wp
        yl = 0.0_wp
        yu = 1.0_wp
        zl = 0.0_wp
        zu = 1.0_wp
        ql = 0.0_wp
        qu = 1.0_wp
        rl = 0.0_wp
        ru = 1.0_wp
        call run_5d_test()

        !============================================
        ! 6d integral tests
        !============================================



    end do

    contains
!************************************************************************************

    !*************************************************************
        subroutine run_test()

        implicit none

        itest = my_int%ifunc

        !set up the class
        call my_int%initialize(fx=test_func, xl=a, xu=b, tolx=tol, methodx=meth)

        !reset number of function evaluations:
        my_int%n_evals = 0

        !integrate the function:
        call my_int%integrate(ans, ierr, err)

        !get the true answer:
        answer = test_integral(my_int,a,b)

        !print results:
        write(*,'(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I9,A,E15.5)') &
                     itest, ',', &
                     '1D',',',&
                     trim(set_of_quadrature_methods(i)%name), ',', &
                     tol, ',', &
                     ans, ',', &
                     ierr, ',', &
                     err, ',', &
                     my_int%n_evals, ',', &
                     answer - ans

        if (abs(answer-ans)>100*tol) error stop 'TEST FAILED'

        end subroutine run_test
    !*************************************************************

    !*************************************************************
        function test_func(me,x) result(f)

        !! The function is f(x) = amp*sin(freq*x + phase)

        implicit none

        class(integration_class_1d),intent(inout)  :: me
        real(wp), intent(in)  :: x
        real(wp)              :: f

        select type (me)
        class is (sin_type)

            select case(me%ifunc)
            case(1)

                f = me%amp * sin (me%freq*x + me%phase) !a sin function

            case(2)

                f = (two*x**5 - x + three) / x**2

            case(3)

                f = three/exp(-x) - one/(three*x)

            case(4) !http://mathworld.wolfram.com/DefiniteIntegral.html

                f = log( two*cos(x/two) )

            case(5)

                f = exp(x)

            case default
                error stop 'Error in test_func: invalid value of ifunc:'
            end select

            me%n_evals = me%n_evals + 1             !number of function evaluations

        class default
            error stop 'Error in test_func: invalid class.'
        end select

        end function test_func
    !*************************************************************

    !*************************************************************
        function test_integral(me,a,b) result(f)

        !! The integral of f(x) = amp*sin(freq*x + phase) from a-->b

        implicit none

        class(integration_class_1d),intent(inout)  :: me
        real(wp), intent(in)  :: a    !! lower limit
        real(wp), intent(in)  :: b    !! upper limit
        real(wp)              :: f

        real(wp) :: c1, c2

        select type (me)
        class is (sin_type)

            select case(me%ifunc)
            case(1)

                c1 = cos( a * me%freq + me%phase )
                c2 = cos( b * me%freq + me%phase )
                f = me%amp * ( c1 - c2 ) / me%freq

            case(2)

                !for bounds: [1,2]
                f = 9.0_wp - log(2.0_wp)

            case(3)

                !for bounds: [-20,-1]
                f = three*exp(-one) - three*exp(-20.0_wp) + one/three*log(abs(20.0_wp))

            case(4) !http://mathworld.wolfram.com/DefiniteIntegral.html

                !for bounds: [-pi,pi]
                f = zero

            case(5)

                f = exp(b) - exp(a)

            case default
                error stop 'Error in test_integral: invalid value of ifunc'
            end select

        class default
            error stop 'Error in test_integral: invalid class.'
        end select

        end function test_integral
    !*************************************************************

    !*************************************************************
        subroutine run_2d_test()

        implicit none

        itest = doub%ifunc

        !set up the class
        call doub%initialize(fxy=test_2d_func,xl=xl,xu=xu,yl=yl,&
                             yu=yu,tolx=tol,toly=tol,methodx=meth,methody=meth)

        !reset number of function evaluations:
        doub%n_evals = 0

        !integrate the function:
        call doub%integrate(ans, ierr, err)

        !get the true answer:
        answer = test_2d_integral(doub,xl,xu,yl,yu)

        !print results:
        write(*,'(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I9,A,E15.5)') &
                     itest, ',', &
                     '2D',',',&
                     trim(set_of_quadrature_methods(i)%name), ',', &
                     tol, ',', &
                     ans, ',', &
                     ierr, ',', &
                     err, ',', &
                     doub%n_evals, ',', &
                     answer - ans

        if (abs(answer-ans)>100*tol) error stop 'TEST FAILED'

        end subroutine run_2d_test
    !*************************************************************

    !*************************************************************
        function test_2d_func(me,x,y) result(f)

        !! The function is f(x,y)

        implicit none

        class(integration_class_2d),intent(inout)   :: me
        real(wp), intent(in)  :: x
        real(wp), intent(in)  :: y
        real(wp)              :: f

        select type (me)
        class is (my_doub)

            select case(me%ifunc)
            case(1)

                f = x * y

            case(2)

                f  = sin(x) + cos(y)

            case default
                error stop 'Error in test_2d_func: invalid value of ifunc'
            end select

            me%n_evals = me%n_evals + 1

        class default
            error stop 'Error in test_2d_func: invalid class.'
        end select

        end function test_2d_func
    !*************************************************************

    !*************************************************************
        function test_2d_integral(me,xl,xu,yl,yu) result(f)

        !! The double integral of f(x,y) dx dy = x*y from xl->xu, yl->yu

        implicit none

        class(integration_class_2d),intent(inout)  :: me
        real(wp), intent(in)  :: xl
        real(wp), intent(in)  :: xu
        real(wp), intent(in)  :: yl
        real(wp), intent(in)  :: yu
        real(wp)              :: f

        select type (me)
        class is (my_doub)

            select case(me%ifunc)
            case(1)

                f = ( xu**2/two - xl**2/two ) * ( yu**2/two - yl**2/two )

            case(2)

                f  = (-cos(xu)+cos(xl))*(yu-yl) + (xu-xl)*(sin(yu)-sin(yl))

            case default
                error stop 'Error in test_2d_integral: invalid value of ifunc'
            end select

        class default
            error stop 'Error in test_2d_integral: invalid class.'
        end select

        end function test_2d_integral
    !*************************************************************

    !*************************************************************
        subroutine run_3d_test()

        implicit none

        itest = trip%ifunc

        !set up the class
        call trip%initialize(fxyz=test_3d_func,xl=xl,xu=xu,yl=yl,yu=yu,zl=zl,zu=zu,&
                             tolx=tol,toly=tol,tolz=tol,&
                             methodx=meth,methody=meth,methodz=meth)

        !reset number of function evaluations:
        trip%n_evals = 0

        !integrate the function:
        call trip%integrate(ans, ierr, err)

        !get the true answer:
        answer = test_3d_integral(trip,xl,xu,yl,yu,zl,zu)

        !print results:
        write(*,'(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I9,A,E15.5)') &
                     itest, ',', &
                     '3D',',',&
                     trim(set_of_quadrature_methods(i)%name), ',', &
                     tol, ',', &
                     ans, ',', &
                     ierr, ',', &
                     err, ',', &
                     doub%n_evals, ',', &
                     answer - ans

        if (abs(answer-ans)>100*tol) error stop 'TEST FAILED'

        end subroutine run_3d_test
    !*************************************************************

    !*************************************************************
        function test_3d_func(me,x,y,z) result(f)

        !! The function is f(x,y,z)

        implicit none

        class(integration_class_3d),intent(inout)   :: me
        real(wp), intent(in)  :: x
        real(wp), intent(in)  :: y
        real(wp), intent(in)  :: z
        real(wp)              :: f

        select type (me)
        class is (my_trip)

            select case(me%ifunc)
            case(1)

                ! see : http://tutorial.math.lamar.edu/Classes/CalcIII/TripleIntegrals.aspx

                f = 8.0_wp * x * y * z

            case default
                error stop 'Error in test_3d_func: invalid value of ifunc'
            end select

            me%n_evals = me%n_evals + 1

        class default
            error stop 'Error in test_3d_func: invalid class.'
        end select

        end function test_3d_func
    !*************************************************************

    !*************************************************************
        function test_3d_integral(me,xl,xu,yl,yu,zl,zu) result(f)

        !! The triple integral of f(x,y,z)

        implicit none

        class(integration_class_3d),intent(inout)  :: me
        real(wp), intent(in)  :: xl
        real(wp), intent(in)  :: xu
        real(wp), intent(in)  :: yl
        real(wp), intent(in)  :: yu
        real(wp), intent(in)  :: zl
        real(wp), intent(in)  :: zu
        real(wp)              :: f

        select type (me)
        class is (my_trip)

            select case(me%ifunc)
            case(1)

                f = 15.0_wp

            case default
                error stop 'Error in test_3d_integral: invalid value of ifunc'
            end select

        class default
            error stop 'Error in test_3d_integral: invalid class.'
        end select

        end function test_3d_integral
    !*************************************************************

    !*************************************************************
        subroutine run_4d_test()

        implicit none

        itest = int_4d%ifunc

        !set up the class
        call int_4d%initialize(fxyzq=test_4d_func,xl=xl,xu=xu,yl=yl,yu=yu,zl=zl,zu=zu,ql=ql,qu=qu,&
                               tolx=tol,toly=tol,tolz=tol,tolq=tol,&
                               methodx=meth,methody=meth,methodz=meth,methodq=meth)

        !reset number of function evaluations:
        int_4d%n_evals = 0

        !integrate the function:
        call int_4d%integrate(ans, ierr, err)

        !get the true answer:
        answer = test_4d_integral(int_4d,xl,xu,yl,yu,zl,zu,ql,qu)

        !print results:
        write(*,'(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I9,A,E15.5)') &
                     itest, ',', &
                     '4D',',',&
                     trim(set_of_quadrature_methods(i)%name), ',', &
                     tol, ',', &
                     ans, ',', &
                     ierr, ',', &
                     err, ',', &
                     int_4d%n_evals, ',', &
                     answer - ans

        if (abs(answer-ans)>100*tol) error stop 'TEST FAILED'

        end subroutine run_4d_test
    !*************************************************************

    !*************************************************************
        function test_4d_func(me,x,y,z,q) result(f)

        !! The function is f(x,y,z,q)

        implicit none

        class(integration_class_4d),intent(inout)   :: me
        real(wp), intent(in)  :: x
        real(wp), intent(in)  :: y
        real(wp), intent(in)  :: z
        real(wp), intent(in)  :: q
        real(wp)              :: f

        select type (me)
        class is (my_4)

            select case(me%ifunc)
            case(1)

                f = 1.0_wp

            case default
                error stop 'Error in test_4d_func: invalid value of ifunc'
            end select

            me%n_evals = me%n_evals + 1

        class default
            error stop 'Error in test_4d_func: invalid class.'
        end select

        end function test_4d_func
    !*************************************************************

    !*************************************************************
        function test_4d_integral(me,xl,xu,yl,yu,zl,zu,ql,qu) result(f)

        !! The 4D integral of f(x,y,z,q)

        implicit none

        class(integration_class_4d),intent(inout)  :: me
        real(wp), intent(in)  :: xl
        real(wp), intent(in)  :: xu
        real(wp), intent(in)  :: yl
        real(wp), intent(in)  :: yu
        real(wp), intent(in)  :: zl
        real(wp), intent(in)  :: zu
        real(wp), intent(in)  :: ql
        real(wp), intent(in)  :: qu
        real(wp)              :: f

        select type (me)
        class is (my_4)

            select case(me%ifunc)
            case(1)

                f = 1.0_wp

            case default
                error stop 'Error in test_4d_integral: invalid value of ifunc'
            end select

        class default
            error stop 'Error in test_4d_integral: invalid class.'
        end select

        end function test_4d_integral
    !*************************************************************


    !*************************************************************
        subroutine run_5d_test()

        implicit none

        itest = int_5d%ifunc

        !set up the class
        call int_5d%initialize(fxyzqr=test_5d_func,xl=xl,xu=xu,yl=yl,yu=yu,zl=zl,zu=zu,ql=ql,qu=qu,rl=rl,ru=ru,&
                               tolx=tol,toly=tol,tolz=tol,tolq=tol,tolr=tol,&
                               methodx=meth,methody=meth,methodz=meth,methodq=meth,methodr=meth)

        !reset number of function evaluations:
        int_5d%n_evals = 0

        !integrate the function:
        call int_5d%integrate(ans, ierr, err)

        !get the true answer:
        answer = test_5d_integral(int_5d,xl,xu,yl,yu,zl,zu,ql,qu,rl,ru)

        !print results:
        write(*,'(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I9,A,E15.5)') &
                     itest, ',', &
                     '5D',',',&
                     trim(set_of_quadrature_methods(i)%name), ',', &
                     tol, ',', &
                     ans, ',', &
                     ierr, ',', &
                     err, ',', &
                     int_5d%n_evals, ',', &
                     answer - ans

        if (abs(answer-ans)>100*tol) error stop 'TEST FAILED'

        end subroutine run_5d_test
    !*************************************************************

    !*************************************************************
        function test_5d_func(me,x,y,z,q,r) result(f)

        !! The function is f(x,y,z,q,r)

        implicit none

        class(integration_class_5d),intent(inout)   :: me
        real(wp), intent(in)  :: x
        real(wp), intent(in)  :: y
        real(wp), intent(in)  :: z
        real(wp), intent(in)  :: q
        real(wp), intent(in)  :: r
        real(wp)              :: f

        select type (me)
        class is (my_5)

            select case(me%ifunc)
            case(1)

                f = 1.0_wp

            case default
                error stop 'Error in test_5d_func: invalid value of ifunc'
            end select

            me%n_evals = me%n_evals + 1

        class default
            error stop 'Error in test_5d_func: invalid class.'
        end select

        end function test_5d_func
    !*************************************************************

    !*************************************************************
        function test_5d_integral(me,xl,xu,yl,yu,zl,zu,ql,qu,rl,ru) result(f)

        !! The 5D integral of f(x,y,z,q,r)

        implicit none

        class(integration_class_5d),intent(inout)  :: me
        real(wp), intent(in)  :: xl
        real(wp), intent(in)  :: xu
        real(wp), intent(in)  :: yl
        real(wp), intent(in)  :: yu
        real(wp), intent(in)  :: zl
        real(wp), intent(in)  :: zu
        real(wp), intent(in)  :: ql
        real(wp), intent(in)  :: qu
        real(wp), intent(in)  :: rl
        real(wp), intent(in)  :: ru
        real(wp)              :: f

        select type (me)
        class is (my_5)

            select case(me%ifunc)
            case(1)

                f = 1.0_wp

            case default
                error stop 'Error in test_5d_integral: invalid value of ifunc'
            end select

        class default
            error stop 'Error in test_5d_integral: invalid class.'
        end select

        end function test_5d_integral
    !*************************************************************


!************************************************************************************
    end program main
!************************************************************************************