!*******************************************************************************************************
!>
!  Integration of functions using adaptive Guassian quadrature.
!
!### Author
!  * Jacob Williams

    module quadrature_module

    use, intrinsic :: iso_fortran_env,  only: wp => real64  ! double precision
  !  use, intrinsic :: iso_fortran_env,  only: wp => real128 ! quad precision

    implicit none

    private

    integer,parameter,public :: quadrature_wp = wp !! export working precision

    !parameters:
    real(wp),parameter :: zero      = 0.0_wp
    real(wp),parameter :: one_half  = 0.5_wp
    real(wp),parameter :: one       = 1.0_wp
    real(wp),parameter :: two       = 2.0_wp
    real(wp),parameter :: three     = 3.0_wp
    real(wp),parameter :: four      = 4.0_wp

    type,public :: quadrature_method
        !! quadrature methods
        integer            :: n_points = 0
        character(len=100) :: name     = ''
    end type quadrature_method

    type(quadrature_method),parameter :: quad_gauss_6  = quadrature_method(6,  'Adaptive 6-point Legendre-Gauss')
    type(quadrature_method),parameter :: quad_gauss_8  = quadrature_method(8,  'Adaptive 8-point Legendre-Gauss')
    type(quadrature_method),parameter :: quad_gauss_10 = quadrature_method(10, 'Adaptive 10-point Legendre-Gauss')
    type(quadrature_method),parameter :: quad_gauss_12 = quadrature_method(12, 'Adaptive 12-point Legendre-Gauss')
    type(quadrature_method),parameter :: quad_gauss_14 = quadrature_method(14, 'Adaptive 14-point Legendre-Gauss')

    type(quadrature_method),dimension(5),parameter,public :: set_of_quadrature_methods = &
                                                                 [ quad_gauss_6 ,&
                                                                   quad_gauss_8 ,&
                                                                   quad_gauss_10,&
                                                                   quad_gauss_12,&
                                                                   quad_gauss_14 ]

    type,abstract,private :: integration_class
        private
    end type integration_class

    type,extends(integration_class),public :: integration_class_1d
        !! single integration class: for 1d integration of `f(x)`

        private

        procedure(func_1d),pointer :: fun => null()  !! function `f(x)` to be integrated

        procedure(gauss_func),pointer :: g => null()  !! the guass quadrature formula to use
        real(wp) :: a      = zero      !! lower limit of integration
        real(wp) :: b      = zero      !! upper limit of integration (may be less than a)
        real(wp) :: tol    = zero      !! the requested relative error tolerance.

        real(wp) :: val = zero   !! the value of `x`. Only used for multiple
                                 !! integration to pass to the inner integrals

        contains

        private

        procedure :: dgauss_generic  !! core integration routine. refactored from
                                     !! SLATEC with selectable quadrature method
        procedure,public :: initialize => initialize_integration_class !! to set up the class
        procedure,public :: integrate  => integrate_1d !! to integrate the function `fun`

    end type integration_class_1d

    type,extends(integration_class),public :: integration_class_2d
        !! double integration class: for 2d integration of `f(x,y)`
        private
        procedure(func_2d),pointer :: fxy => null()   !! function `f(x,y)` to be integrated
        type(integration_class_1d) :: iy     !! for the dy integration
        type(integration_class_1d) :: ix     !! for the dx integration
        contains
        private
        procedure,public :: initialize => initialize_integration_class_2d !! to set up the class
        procedure,public :: integrate  => integrate_2d !! to integrate the function `fxy`
    end type integration_class_2d

    type,extends(integration_class),public :: integration_class_3d
        !! double integration class: for 3d integration of `f(x,y,z)`
        private
        procedure(func_3d),pointer :: fxyz => null()   !! function `f(x,y,z)` to be integrated
        type(integration_class_1d) :: iz     !! for the dz integration
        type(integration_class_1d) :: iy     !! for the dy integration
        type(integration_class_1d) :: ix     !! for the dx integration
        contains
        private
        procedure,public :: initialize => initialize_integration_class_3d !! to set up the class
        procedure,public :: integrate  => integrate_3d !! to integrate the function `fxyz`
    end type integration_class_3d

    type,extends(integration_class),public :: integration_class_4d
        !! double integration class: for 4d integration of `f(x,y,z,q)`
        private
        procedure(func_4d),pointer :: fxyzq => null()   !! function `f(x,y,z,q)` to be integrated
        type(integration_class_1d) :: iq     !! for the dq integration
        type(integration_class_1d) :: iz     !! for the dz integration
        type(integration_class_1d) :: iy     !! for the dy integration
        type(integration_class_1d) :: ix     !! for the dx integration
        contains
        private
        procedure,public :: initialize => initialize_integration_class_4d !! to set up the class
        procedure,public :: integrate  => integrate_4d !! to integrate the function `fxyzq`
    end type integration_class_4d

    type,extends(integration_class),public :: integration_class_5d
        !! double integration class: for 5d integration of `f(x,y,z,q,r)`
        private
        procedure(func_5d),pointer :: fxyzqr => null()   !! function `f(x,y,z,q,r)` to be integrated
        type(integration_class_1d) :: ir     !! for the dr integration
        type(integration_class_1d) :: iq     !! for the dq integration
        type(integration_class_1d) :: iz     !! for the dz integration
        type(integration_class_1d) :: iy     !! for the dy integration
        type(integration_class_1d) :: ix     !! for the dx integration
        contains
        private
        procedure,public :: initialize => initialize_integration_class_5d !! to set up the class
        procedure,public :: integrate  => integrate_5d !! to integrate the function `fxyzqr`
    end type integration_class_5d

    type,extends(integration_class),public :: integration_class_6d
        !! double integration class: for 6d integration of `f(x,y,z,q,r,s)`
        private
        procedure(func_6d),pointer :: fxyzqrs => null()   !! function `f(x,y,z,q,r,s)` to be integrated
        type(integration_class_1d) :: is     !! for the ds integration
        type(integration_class_1d) :: ir     !! for the dr integration
        type(integration_class_1d) :: iq     !! for the dq integration
        type(integration_class_1d) :: iz     !! for the dz integration
        type(integration_class_1d) :: iy     !! for the dy integration
        type(integration_class_1d) :: ix     !! for the dx integration
        contains
        private
        procedure,public :: initialize => initialize_integration_class_6d !! to set up the class
        procedure,public :: integrate  => integrate_6d !! to integrate the function `fxyzqrs`
    end type integration_class_6d

    abstract interface

        function func_1d(me,x) result(f)
            !! 1d user function f(x)
            import :: wp,integration_class_1d
            implicit none
            class(integration_class_1d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp)                                  :: f
        end function func_1d

        function func_2d(me,x,y) result(f)
            !! 2d user function f(x,y)
            import :: wp,integration_class_2d
            implicit none
            class(integration_class_2d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp), intent(in)                      :: y
            real(wp)                                  :: f
        end function func_2d

        function func_3d(me,x,y,z) result(f)
            !! 3d user function f(x,y,z)
            import :: wp,integration_class_3d
            implicit none
            class(integration_class_3d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp), intent(in)                      :: y
            real(wp), intent(in)                      :: z
            real(wp)                                  :: f
        end function func_3d

        function func_4d(me,x,y,z,q) result(f)
            !! 4d user function f(x,y,z,q)
            import :: wp,integration_class_4d
            implicit none
            class(integration_class_4d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp), intent(in)                      :: y
            real(wp), intent(in)                      :: z
            real(wp), intent(in)                      :: q
            real(wp)                                  :: f
        end function func_4d

        function func_5d(me,x,y,z,q,r) result(f)
            !! 5d user function f(x,y,z,q,r)
            import :: wp,integration_class_5d
            implicit none
            class(integration_class_5d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp), intent(in)                      :: y
            real(wp), intent(in)                      :: z
            real(wp), intent(in)                      :: q
            real(wp), intent(in)                      :: r
            real(wp)                                  :: f
        end function func_5d

        function func_6d(me,x,y,z,q,r,s) result(f)
            !! 6d user function f(x,y,z,q,r,s)
            import :: wp,integration_class_6d
            implicit none
            class(integration_class_6d),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp), intent(in)                      :: y
            real(wp), intent(in)                      :: z
            real(wp), intent(in)                      :: q
            real(wp), intent(in)                      :: r
            real(wp), intent(in)                      :: s
            real(wp)                                  :: f
        end function func_6d

        function gauss_func(me, x, h) result(f)
            !! guass quadrature formula
            import :: wp,integration_class_1d
            implicit none
            class(integration_class_1d),intent(inout)  :: me
            real(wp), intent(in)                    :: x
            real(wp), intent(in)                    :: h
            real(wp)                                :: f
        end function gauss_func

    end interface

    contains
!*******************************************************************************************************

!********************************************************************************
!>
!  Initialize the 1D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class(me,fx,&
                                            xl,xu,tolx,methodx)

    implicit none

    class(integration_class_1d),intent(inout) :: me
    procedure(func_1d)    :: fx       !! 1d function: f(x)
    real(wp),intent(in)   :: xl       !! x integration lower bound
    real(wp),intent(in)   :: xu       !! x integration upper bound
    real(wp),intent(in)   :: tolx     !! error tolerance for dx integration
    integer,intent(in)    :: methodx  !! quadrature method to use for x

    ! select quadrature rule
    select case (methodx)
    case(6);  me%g => g6
    case(8);  me%g => g8
    case(10); me%g => g10
    case(12); me%g => g12
    case(14); me%g => g14
    case default
        error stop 'invalid quadrature method in initialize_integration_class'
    end select

    me%fun  => fx       !the function f(x) to integrate
    me%tol  = tolx      !tolerance
    me%a    = xl        !lower bound
    me%b    = xu        !upper bound

    end subroutine initialize_integration_class
!********************************************************************************

!********************************************************************************
!>
!  Initialize the 2D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class_2d(me,fxy,&
                                               xl,xu,yl,yu,&
                                               tolx,toly,&
                                               methodx,methody)

    implicit none

    class(integration_class_2d),intent(inout) :: me
    procedure(func_2d)     :: fxy      !! 2d function: f(x,y)
    real(wp),intent(in)    :: xl       !! x integration lower bound
    real(wp),intent(in)    :: xu       !! x integration upper bound
    real(wp),intent(in)    :: yl       !! y integration lower bound
    real(wp),intent(in)    :: yu       !! y integration upper bound
    real(wp),intent(in)    :: tolx     !! error tolerance for dx integration
    real(wp),intent(in)    :: toly     !! error tolerance for dy integration
    integer,intent(in)     :: methodx  !! quadrature method to use for x
    integer,intent(in)     :: methody  !! quadrature method to use for y

    procedure(func_1d),pointer :: dummy => null() !! these will be set in [[integrate_2d]]

    me%fxy => fxy  !the user-defined f(x,y) function to integrate

    ! individual integrators:
    call me%ix%initialize(dummy,xl,xu,tolx,methodx)
    call me%iy%initialize(dummy,yl,yu,toly,methody)

    end subroutine initialize_integration_class_2d
!********************************************************************************

!********************************************************************************
!>
!  Initialize the 3D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class_3d(me,fxyz,&
                                               xl,xu,yl,yu,zl,zu,&
                                               tolx,toly,tolz,&
                                               methodx,methody,methodz)

    implicit none

    class(integration_class_3d),intent(inout) :: me
    procedure(func_3d)     :: fxyz     !! 3d function: f(x,y,z)
    real(wp),intent(in)    :: xl       !! x integration lower bound
    real(wp),intent(in)    :: xu       !! x integration upper bound
    real(wp),intent(in)    :: yl       !! y integration lower bound
    real(wp),intent(in)    :: yu       !! y integration upper bound
    real(wp),intent(in)    :: zl       !! z integration lower bound
    real(wp),intent(in)    :: zu       !! z integration upper bound
    real(wp),intent(in)    :: tolx     !! error tolerance for dx integration
    real(wp),intent(in)    :: toly     !! error tolerance for dy integration
    real(wp),intent(in)    :: tolz     !! error tolerance for dz integration
    integer,intent(in)     :: methodx  !! quadrature method to use for x
    integer,intent(in)     :: methody  !! quadrature method to use for y
    integer,intent(in)     :: methodz  !! quadrature method to use for z

    procedure(func_1d),pointer :: dummy => null() !! these will be set in [[integrate_3d]]

    me%fxyz => fxyz  !the user-defined f(x,y,z) function to integrate

    ! individual integrators:
    call me%ix%initialize(dummy,xl,xu,tolx,methodx)
    call me%iy%initialize(dummy,yl,yu,toly,methody)
    call me%iz%initialize(dummy,zl,zu,tolz,methodz)

    end subroutine initialize_integration_class_3d
!********************************************************************************

!********************************************************************************
!>
!  Initialize the 4D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class_4d(me,fxyzq,&
                                               xl,xu,yl,yu,zl,zu,ql,qu,&
                                               tolx,toly,tolz,tolq,&
                                               methodx,methody,methodz,methodq)

    implicit none

    class(integration_class_4d),intent(inout) :: me
    procedure(func_4d)     :: fxyzq    !! 4d function: f(x,y,z,q)
    real(wp),intent(in)    :: xl       !! x integration lower bound
    real(wp),intent(in)    :: xu       !! x integration upper bound
    real(wp),intent(in)    :: yl       !! y integration lower bound
    real(wp),intent(in)    :: yu       !! y integration upper bound
    real(wp),intent(in)    :: zl       !! z integration lower bound
    real(wp),intent(in)    :: zu       !! z integration upper bound
    real(wp),intent(in)    :: ql       !! q integration lower bound
    real(wp),intent(in)    :: qu       !! q integration upper bound
    real(wp),intent(in)    :: tolx     !! error tolerance for dx integration
    real(wp),intent(in)    :: toly     !! error tolerance for dy integration
    real(wp),intent(in)    :: tolz     !! error tolerance for dz integration
    real(wp),intent(in)    :: tolq     !! error tolerance for dq integration
    integer,intent(in)     :: methodx  !! quadrature method to use for x
    integer,intent(in)     :: methody  !! quadrature method to use for y
    integer,intent(in)     :: methodz  !! quadrature method to use for z
    integer,intent(in)     :: methodq  !! quadrature method to use for q

    procedure(func_1d),pointer :: dummy => null() !! these will be set in [[integrate_3d]]

    me%fxyzq => fxyzq  !the user-defined f(x,y,z,q) function to integrate

    ! individual integrators:
    call me%ix%initialize(dummy,xl,xu,tolx,methodx)
    call me%iy%initialize(dummy,yl,yu,toly,methody)
    call me%iz%initialize(dummy,zl,zu,tolz,methodz)
    call me%iq%initialize(dummy,ql,qu,tolq,methodq)

    end subroutine initialize_integration_class_4d
!********************************************************************************

!********************************************************************************
!>
!  Initialize the 5D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class_5d(me,fxyzqr,&
                                               xl,xu,yl,yu,zl,zu,ql,qu,rl,ru,&
                                               tolx,toly,tolz,tolq,tolr,&
                                               methodx,methody,methodz,methodq,methodr)

    implicit none

    class(integration_class_5d),intent(inout) :: me
    procedure(func_5d)     :: fxyzqr   !! 5d function: f(x,y,z,q,r)
    real(wp),intent(in)    :: xl       !! x integration lower bound
    real(wp),intent(in)    :: xu       !! x integration upper bound
    real(wp),intent(in)    :: yl       !! y integration lower bound
    real(wp),intent(in)    :: yu       !! y integration upper bound
    real(wp),intent(in)    :: zl       !! z integration lower bound
    real(wp),intent(in)    :: zu       !! z integration upper bound
    real(wp),intent(in)    :: ql       !! q integration lower bound
    real(wp),intent(in)    :: qu       !! q integration upper bound
    real(wp),intent(in)    :: rl       !! r integration lower bound
    real(wp),intent(in)    :: ru       !! r integration upper bound
    real(wp),intent(in)    :: tolx     !! error tolerance for dx integration
    real(wp),intent(in)    :: toly     !! error tolerance for dy integration
    real(wp),intent(in)    :: tolz     !! error tolerance for dz integration
    real(wp),intent(in)    :: tolq     !! error tolerance for dq integration
    real(wp),intent(in)    :: tolr     !! error tolerance for dr integration
    integer,intent(in)     :: methodx  !! quadrature method to use for x
    integer,intent(in)     :: methody  !! quadrature method to use for y
    integer,intent(in)     :: methodz  !! quadrature method to use for z
    integer,intent(in)     :: methodq  !! quadrature method to use for q
    integer,intent(in)     :: methodr  !! quadrature method to use for r

    procedure(func_1d),pointer :: dummy => null() !! these will be set in [[integrate_3d]]

    me%fxyzqr => fxyzqr  !the user-defined f(x,y,z,q,r) function to integrate

    ! individual integrators:
    call me%ix%initialize(dummy,xl,xu,tolx,methodx)
    call me%iy%initialize(dummy,yl,yu,toly,methody)
    call me%iz%initialize(dummy,zl,zu,tolz,methodz)
    call me%iq%initialize(dummy,ql,qu,tolq,methodq)
    call me%ir%initialize(dummy,rl,ru,tolr,methodr)

    end subroutine initialize_integration_class_5d
!********************************************************************************

!********************************************************************************
!>
!  Initialize the 6D integration class.
!  Must be called before integration is performed.

    subroutine initialize_integration_class_6d(me,fxyzqrs,&
                                               xl,xu,yl,yu,zl,zu,ql,qu,rl,ru,sl,su,&
                                               tolx,toly,tolz,tolq,tolr,tols,&
                                               methodx,methody,methodz,methodq,methodr,methods)

    implicit none

    class(integration_class_6d),intent(inout) :: me
    procedure(func_6d)     :: fxyzqrs  !! 6d function: f(x,y,z,q,r,s)
    real(wp),intent(in)    :: xl       !! x integration lower bound
    real(wp),intent(in)    :: xu       !! x integration upper bound
    real(wp),intent(in)    :: yl       !! y integration lower bound
    real(wp),intent(in)    :: yu       !! y integration upper bound
    real(wp),intent(in)    :: zl       !! z integration lower bound
    real(wp),intent(in)    :: zu       !! z integration upper bound
    real(wp),intent(in)    :: ql       !! q integration lower bound
    real(wp),intent(in)    :: qu       !! q integration upper bound
    real(wp),intent(in)    :: rl       !! r integration lower bound
    real(wp),intent(in)    :: ru       !! r integration upper bound
    real(wp),intent(in)    :: sl       !! s integration lower bound
    real(wp),intent(in)    :: su       !! s integration upper bound
    real(wp),intent(in)    :: tolx     !! error tolerance for dx integration
    real(wp),intent(in)    :: toly     !! error tolerance for dy integration
    real(wp),intent(in)    :: tolz     !! error tolerance for dz integration
    real(wp),intent(in)    :: tolq     !! error tolerance for dq integration
    real(wp),intent(in)    :: tolr     !! error tolerance for dr integration
    real(wp),intent(in)    :: tols     !! error tolerance for ds integration
    integer,intent(in)     :: methodx  !! quadrature method to use for x
    integer,intent(in)     :: methody  !! quadrature method to use for y
    integer,intent(in)     :: methodz  !! quadrature method to use for z
    integer,intent(in)     :: methodq  !! quadrature method to use for q
    integer,intent(in)     :: methodr  !! quadrature method to use for r
    integer,intent(in)     :: methods  !! quadrature method to use for s

    procedure(func_1d),pointer :: dummy => null() !! these will be set in [[integrate_3d]]

    me%fxyzqrs => fxyzqrs  !the user-defined f(x,y,z,q,r,s) function to integrate

    ! individual integrators:
    call me%ix%initialize(dummy,xl,xu,tolx,methodx)
    call me%iy%initialize(dummy,yl,yu,toly,methody)
    call me%iz%initialize(dummy,zl,zu,tolz,methodz)
    call me%iq%initialize(dummy,ql,qu,tolq,methodq)
    call me%ir%initialize(dummy,rl,ru,tolr,methodr)
    call me%is%initialize(dummy,sl,su,tols,methods)

    end subroutine initialize_integration_class_6d
!********************************************************************************

!********************************************************************************
!>
!   Perform the 1D integration.

    subroutine integrate_1d (me, ans, ierr, err)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    !call the low-level routine:
    call me%dgauss_generic(me%a, me%b, me%tol, ans, ierr, err)

    end subroutine integrate_1d
!********************************************************************************

!********************************************************************************
!>
!   Perform the 2D integration.

    subroutine integrate_2d (me, ans, ierr, err)

    implicit none

    class(integration_class_2d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    ! set the two functions to the contained wrappers:
    me%iy%fun => f_of_y
    me%ix%fun => f_of_x

    ! call the low-level routine:
    call me%iy%integrate(ans, ierr, err)

    contains

        function f_of_x(ix,x) result(f)
            class(integration_class_1d),intent(inout) :: ix
            real(wp), intent(in) :: x
            real(wp)             :: f
            f = me%fxy(x,me%iy%val)
        end function f_of_x

        function f_of_y(iy,y) result(f)
            class(integration_class_1d),intent(inout)  :: iy
            real(wp), intent(in) :: y
            real(wp)             :: f
            iy%val = y  !set y value
            call me%ix%integrate(f, ierr, err)
        end function f_of_y

    end subroutine integrate_2d
!********************************************************************************

!********************************************************************************
!>
!  Perform the 3D integration.

    subroutine integrate_3d (me, ans, ierr, err)

    implicit none

    class(integration_class_3d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    ! set the functions to the contained wrappers:
    me%iz%fun => f_of_z
    me%iy%fun => f_of_y
    me%ix%fun => f_of_x

    ! call the low-level routine:
    call me%iz%integrate(ans, ierr, err)

    contains

        function f_of_x(ix,x) result(f)
            class(integration_class_1d),intent(inout) :: ix
            real(wp), intent(in) :: x
            real(wp)             :: f
            f = me%fxyz(x,me%iy%val,me%iz%val)
        end function f_of_x

        function f_of_y(iy,y) result(f)
            class(integration_class_1d),intent(inout)  :: iy
            real(wp), intent(in) :: y
            real(wp)             :: f
            iy%val = y
            call me%ix%integrate(f, ierr, err)
        end function f_of_y

        function f_of_z(iz,z) result(f)
            class(integration_class_1d),intent(inout)  :: iz
            real(wp), intent(in) :: z
            real(wp)             :: f
            iz%val = z
            call me%iy%integrate(f, ierr, err)
        end function f_of_z

    end subroutine integrate_3d
!********************************************************************************

!********************************************************************************
!>
!  Perform the 4D integration.

    subroutine integrate_4d (me, ans, ierr, err)

    implicit none

    class(integration_class_4d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    ! set the functions to the contained wrappers:
    me%iq%fun => f_of_q
    me%iz%fun => f_of_z
    me%iy%fun => f_of_y
    me%ix%fun => f_of_x

    ! call the low-level routine:
    call me%iq%integrate(ans, ierr, err)

    contains

        function f_of_x(ix,x) result(f)
            class(integration_class_1d),intent(inout) :: ix
            real(wp), intent(in) :: x
            real(wp)             :: f
            f = me%fxyzq(x,me%iy%val,me%iz%val,me%iq%val)
        end function f_of_x

        function f_of_y(iy,y) result(f)
            class(integration_class_1d),intent(inout)  :: iy
            real(wp), intent(in) :: y
            real(wp)             :: f
            iy%val = y
            call me%ix%integrate(f, ierr, err)
        end function f_of_y

        function f_of_z(iz,z) result(f)
            class(integration_class_1d),intent(inout)  :: iz
            real(wp), intent(in) :: z
            real(wp)             :: f
            iz%val = z
            call me%iy%integrate(f, ierr, err)
        end function f_of_z

        function f_of_q(iq,q) result(f)
            class(integration_class_1d),intent(inout)  :: iq
            real(wp), intent(in) :: q
            real(wp)             :: f
            iq%val = q
            call me%iz%integrate(f, ierr, err)
        end function f_of_q

    end subroutine integrate_4d
!********************************************************************************

!********************************************************************************
!>
!  Perform the 5D integration.

    subroutine integrate_5d (me, ans, ierr, err)

    implicit none

    class(integration_class_5d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    ! set the functions to the contained wrappers:
    me%ir%fun => f_of_r
    me%iq%fun => f_of_q
    me%iz%fun => f_of_z
    me%iy%fun => f_of_y
    me%ix%fun => f_of_x

    ! call the low-level routine:
    call me%ir%integrate(ans, ierr, err)

    contains

        function f_of_x(ix,x) result(f)
            class(integration_class_1d),intent(inout) :: ix
            real(wp), intent(in) :: x
            real(wp)             :: f
            f = me%fxyzqr(x,me%iy%val,me%iz%val,me%iq%val,me%ir%val)
        end function f_of_x

        function f_of_y(iy,y) result(f)
            class(integration_class_1d),intent(inout)  :: iy
            real(wp), intent(in) :: y
            real(wp)             :: f
            iy%val = y
            call me%ix%integrate(f, ierr, err)
        end function f_of_y

        function f_of_z(iz,z) result(f)
            class(integration_class_1d),intent(inout)  :: iz
            real(wp), intent(in) :: z
            real(wp)             :: f
            iz%val = z
            call me%iy%integrate(f, ierr, err)
        end function f_of_z

        function f_of_q(iq,q) result(f)
            class(integration_class_1d),intent(inout)  :: iq
            real(wp), intent(in) :: q
            real(wp)             :: f
            iq%val = q
            call me%iz%integrate(f, ierr, err)
        end function f_of_q

        function f_of_r(ir,r) result(f)
            class(integration_class_1d),intent(inout)  :: ir
            real(wp), intent(in) :: r
            real(wp)             :: f
            ir%val = r
            call me%iq%integrate(f, ierr, err)
        end function f_of_r

    end subroutine integrate_5d
!********************************************************************************

!********************************************************************************
!>
!  Perform the 6D integration.

    subroutine integrate_6d (me, ans, ierr, err)

    implicit none

    class(integration_class_6d),intent(inout)  :: me
    real(wp),intent(out)  :: ans
    integer,intent(out)   :: ierr
    real(wp),intent(out)  :: err

    ! set the functions to the contained wrappers:
    me%is%fun => f_of_s
    me%ir%fun => f_of_r
    me%iq%fun => f_of_q
    me%iz%fun => f_of_z
    me%iy%fun => f_of_y
    me%ix%fun => f_of_x

    ! call the low-level routine:
    call me%is%integrate(ans, ierr, err)

    contains

        function f_of_x(ix,x) result(f)
            class(integration_class_1d),intent(inout) :: ix
            real(wp), intent(in) :: x
            real(wp)             :: f
            f = me%fxyzqrs(x,me%iy%val,me%iz%val,me%iq%val,me%ir%val,me%is%val)
        end function f_of_x

        function f_of_y(iy,y) result(f)
            class(integration_class_1d),intent(inout)  :: iy
            real(wp), intent(in) :: y
            real(wp)             :: f
            iy%val = y
            call me%ix%integrate(f, ierr, err)
        end function f_of_y

        function f_of_z(iz,z) result(f)
            class(integration_class_1d),intent(inout)  :: iz
            real(wp), intent(in) :: z
            real(wp)             :: f
            iz%val = z
            call me%iy%integrate(f, ierr, err)
        end function f_of_z

        function f_of_q(iq,q) result(f)
            class(integration_class_1d),intent(inout)  :: iq
            real(wp), intent(in) :: q
            real(wp)             :: f
            iq%val = q
            call me%iz%integrate(f, ierr, err)
        end function f_of_q

        function f_of_r(ir,r) result(f)
            class(integration_class_1d),intent(inout)  :: ir
            real(wp), intent(in) :: r
            real(wp)             :: f
            ir%val = r
            call me%iq%integrate(f, ierr, err)
        end function f_of_r

        function f_of_s(is,s) result(f)
            class(integration_class_1d),intent(inout)  :: is
            real(wp), intent(in) :: s
            real(wp)             :: f
            is%val = s
            call me%ir%integrate(f, ierr, err)
        end function f_of_s

    end subroutine integrate_6d
!********************************************************************************

!********************************************************************************
!>
!  Integrate a real function of one variable over a finite
!  interval using the specified adaptive algorithm.
!  Intended primarily for high accuracy
!  integration or integration of smooth functions.
!
!### License
!  * SLATEC is public domain software: http://www.netlib.org/slatec/guide
!
!### See also
!  * Original sourcecode from: http://www.netlib.org/slatec/src/dgaus8.f
!
!### Author
!  * Jones, R. E., (SNLA) -- Original SLATEC code.
!  * Jacob Williams : 1/20/2020 : refactored to modern Fortran and generalized.
!
!@note This function is recursive.
!      [It can call itself indirectly during double integration]

    recursive subroutine dgauss_generic (me, lb, ub, error_tol, ans, ierr, err)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp),intent(in)   :: lb         !! lower bound of the integration
    real(wp),intent(in)   :: ub         !! upper bound of the integration
    real(wp),intent(in)   :: error_tol  !! is a requested pseudorelative error tolerance.  normally
                                        !! pick a value of abs(error_tol) so that
                                        !! dtol < abs(error_tol) <= 1.0e-3 where dtol is the larger
                                        !! of 1.0e-18 and the real(wp) unit roundoff d1mach(4).
                                        !! ans will normally have no more error than abs(error_tol)
                                        !! times the integral of the absolute value of fun(x).  usually,
                                        !! smaller values of error_tol yield more accuracy and require
                                        !! more function evaluations.
    real(wp),intent(out)  :: ans        !! computed value of integral
    integer,intent(out)   :: ierr       !! status code:
                                        !!
                                        !!  * normal codes:
                                        !!    * 1 : `ans` most likely meets requested error tolerance,
                                        !!      or `lb=ub`.
                                        !!    * -1 : `lb` and `ub` are too nearly equal to allow normal
                                        !!      integration. `ans` is set to zero.
                                        !!  * abnormal code:
                                        !!    * 2 : `ans` probably does not meet requested error tolerance.
    real(wp),intent(out)  :: err        !! an estimate of the absolute error in `ans`.
                                        !! the estimated error is solely for information to the user and
                                        !! should not be used as a correction to the computed integral.

    real(wp),parameter  :: sq2      = sqrt(two)
    real(wp),parameter  :: ln2      = log(two)
    integer,parameter   :: nlmn     = 1                   !! ??
    integer,parameter   :: kmx      = 5000                !! ??
    integer,parameter   :: kml      = 6                   !! ??
    real(wp),parameter  :: magic    = 0.30102000_wp       !! ??
    integer,parameter   :: iwork    = 60                  !! size of the work arrays. ?? Why 60 ??
    real(wp),parameter  :: bb       = radix(one)          !! machine constant
    real(wp),parameter  :: d1mach4  = bb**(1-digits(one)) !! machine constant
    real(wp),parameter  :: d1mach5  = log10(bb)           !! machine constant

    integer                   :: k,l,lmn,lmx,mxl,nbits,nib,nlmx
    real(wp)                  :: ae,anib,area,c,ee,ef,eps,est,gl,glr,tol
    real(wp),dimension(iwork) :: aa,hh,vl,gr
    integer,dimension(iwork)  :: lr

    ans = zero
    ierr = 1
    err = zero
    if (lb == ub) return
    aa = zero
    hh = zero
    vl = zero
    gr = zero
    lr = 0
    k = digits(one)
    anib = d1mach5*k/magic
    nbits = anib
    nlmx = min(60,(nbits*5)/8)         ! ... is this the same 60 as iwork???
    lmx = nlmx
    lmn = nlmn
    if (ub /= zero) then
        if (sign(one,ub)*lb > zero) then
            c = abs(one-lb/ub)
            if (c <= 0.1_wp) then
                if (c <= zero) return
                anib = one_half - log(c)/ln2
                nib = anib
                lmx = min(nlmx,nbits-nib-7)
                if (lmx < 1) then
                    ! lb and ub are too nearly equal to allow
                    ! normal integration [ans is set to zero]
                    ierr = -1
                    return
                end if
                lmn = min(lmn,lmx)
            end if
        end if
    end if
    tol = max(abs(error_tol),two**(5-nbits))/two
    if (error_tol == zero) tol = sqrt(d1mach4)
    eps = tol
    hh(1) = (ub-lb)/four
    aa(1) = lb
    lr(1) = 1
    l = 1
    est = me%g(aa(l)+two*hh(l),two*hh(l))
    k = 8
    area = abs(est)
    ef = one_half
    mxl = 0

    !compute refined estimates, estimate the error, etc.
    main : do

        gl = me%g(aa(l)+hh(l),hh(l))
        gr(l) = me%g(aa(l)+three*hh(l),hh(l))
        k = k + 16
        area = area + (abs(gl)+abs(gr(l))-abs(est))
        glr = gl + gr(l)
        ee = abs(est-glr)*ef
        ae = max(eps*area,tol*abs(glr))
        if (ee-ae > zero) then
            !consider the left half of this level
            if (k > kmx) lmx = kml
            if (l >= lmx) then
                mxl = 1
            else
                l = l + 1
                eps = eps*one_half
                ef = ef/sq2
                hh(l) = hh(l-1)*one_half
                lr(l) = -1
                aa(l) = aa(l-1)
                est = gl
                cycle main
            end if
        end if

        err = err + (est-glr)
        if (lr(l) > 0) then
            !return one level
            ans = glr
            do
                if (l <= 1) exit main ! finished
                l = l - 1
                eps = eps*two
                ef = ef*sq2
                if (lr(l) <= 0) then
                    vl(l) = vl(l+1) + ans
                    est = gr(l-1)
                    lr(l) = 1
                    aa(l) = aa(l) + four*hh(l)
                    cycle main
                end if
                ans = vl(l+1) + ans
            end do
        else
            !proceed to right half at this level
            vl(l) = glr
            est = gr(l-1)
            lr(l) = 1
            aa(l) = aa(l) + four*hh(l)
            cycle main
        end if

    end do main

    if ((mxl/=0) .and. (abs(err)>two*tol*area)) ierr = 2 ! ans is probably insufficiently accurate

    end subroutine dgauss_generic
!********************************************************************************

!************************************************************************************
!>
!  6-point method.
!
!### See also
!  * Coefficients from:
!    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

    function g6(me, x, h) result(f)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp), intent(in)                    :: x
    real(wp), intent(in)                    :: h
    real(wp)                                :: f

    !> abscissae:
    real(wp),dimension(3),parameter ::  a = [   0.6612093864662645136613995950199053470064485643&
                                                &951700708145267058521834966071431009442864037464&
                                                &614564298883716392751466795573467722253804381723&
                                                &198010093367423918538864300079016299442625145884&
                                                &9024557188219703863032236201173523213570221879361&
                                                &8906974301231555871064213101639896769013566165126&
                                                &1150514997832_wp,&
                                                &0.2386191860831969086305017216807119354186106301&
                                                &400213501813951645742749342756398422492244272573&
                                                &491316090722230970106872029554530350772051352628&
                                                &872175189982985139866216812636229030578298770859&
                                                &440976999298617585739469216136216592222334626416&
                                                &400139367778945327871453246721518889993399000945&
                                                &408150514997832_wp,&
                                                &0.9324695142031520278123015544939946091347657377&
                                                &122898248725496165266135008442001962762887399219&
                                                &259850478636797265728341065879713795116384041921&
                                                &786180750210169211578452038930846310372961174632&
                                                &524612619760497437974074226320896716211721783852&
                                                &305051047442772222093863676553669179038880252326&
                                                &771150514997832_wp ]
    !> weights:
    real(wp),dimension(3),parameter ::  w = [   0.36076157304813860756983351383771611166152189274&
                                                &674548228973924023714003783726171832096220198881&
                                                &934794311720914037079858987989027836432107077678&
                                                &721140858189221145027225257577711260007323688285&
                                                &916316028951118005174081368554707448247248610118&
                                                &325993144981721640242558677752676819993095031068&
                                                &73150514997832_wp,&
                                                0.46791393457269104738987034398955099481165560576&
                                                &921053531162531996391420162039812703111009258479&
                                                &198230476626878975479710092836255417350295459356&
                                                &355927338665933648259263825590180302812735635025&
                                                &362417046193182590009975698709590053347408007463&
                                                &437682443180817320636917410341626176534629278889&
                                                &17150514997832_wp,&
                                                0.17132449237917034504029614217273289352682250148&
                                                &404398239863543979894576054234015464792770542638&
                                                &866975211652206987440430919174716746217597462964&
                                                &922931803144845206713510916832108437179940676688&
                                                &721266924855699404815942932735702498405343382418&
                                                &236324411837461039120523911904421970357029774978&
                                                &12150514997832_wp ]

    f = h * ( w(1)*( me%fun(x-a(1)*h) + me%fun(x+a(1)*h) ) + &
              w(2)*( me%fun(x-a(2)*h) + me%fun(x+a(2)*h) ) + &
              w(3)*( me%fun(x-a(3)*h) + me%fun(x+a(3)*h) ) )

    end function g6
!************************************************************************************

!************************************************************************************
!>
!  This is the 8-point formula from the original SLATEC routine
!  [DGAUS8](http://www.netlib.org/slatec/src/dgaus8.f).
!
!@note replaced coefficients with high-precision ones from:
!      http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

    function g8(me, x, h) result(f)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp), intent(in) :: x
    real(wp), intent(in) :: h
    real(wp)             :: f

    !> abscissae:
    real(wp),parameter ::   x1 = 0.18343464249564980493947614236018398066675781291297378231718847&
                                                &369920447422154211411606822371112335374526765876&
                                                &428676660891960125238768656837885699951606635681&
                                                &044755516171385019663858107642055323708826547494&
                                                &928123149612477646193635627706457164566131594051&
                                                &34052985058171969174306064445289638150514997832_wp
    real(wp),parameter ::   x2 = 0.52553240991632898581773904918924634904196424312039285775085709&
                                                &927245482076856127252396140019363198206190968292&
                                                &482526085071087937666387799398053953036682536311&
                                                &190182730324023600607174700061279014795875767562&
                                                &412888953366196435283308256242634705401842246036&
                                                &88817537938539658502113876953598879150514997832_wp
    real(wp),parameter ::   x3 = 0.79666647741362673959155393647583043683717173161596483207017029&
                                                &503921730567647309214715192729572593901919745345&
                                                &309730926536564949170108596027725620746216896761&
                                                &539350162903423256455826342053015458560600957273&
                                                &426035574157612651404288519573419337108037227831&
                                                &36113628137267630651413319993338002150514997832_wp
    real(wp),parameter ::   x4 = 0.96028985649753623168356086856947299042823523430145203827163977&
                                                &737242489774341928443943895926331226831042439281&
                                                &729417621023895815521712854793736422049096997004&
                                                &339826183266373468087812635533469278673596634808&
                                                &705975425476039293185338665681328688426134748962&
                                                &8923208763998895240977248938732425615051499783203_wp

    !> weights:
    real(wp),parameter ::   w1 = 0.36268378337836198296515044927719561219414603989433054052482306&
                                                &756668673472390667732436604208482850955025876992&
                                                &629670655292582155698951738449955760078620768427&
                                                &783503828625463057710075533732697147148942683287&
                                                &804318227790778467229655355481996014024877675059&
                                                &28976560993309027632737537826127502150514997832_wp
    real(wp),parameter ::   w2 = 0.31370664587788728733796220198660131326032899900273493769026394&
                                                &507495627194217349696169807623392855604942757464&
                                                &107780861624724683226556160568906242764697589946&
                                                &225031187765625594632872220215204316264677947216&
                                                &038226012952768986525097231851579983531560624197&
                                                &51736972560423953923732838789657919150514997832_wp
    real(wp),parameter ::   w3 = 0.22238103445337447054435599442624088443013087005124956472590928&
                                                &929361681457044904085365314237719792784215926610&
                                                &121221812311143757985257224193818266745320905779&
                                                &086132895368404027893986488760043856972021574820&
                                                &632532471955902286315706513199655897335454406059&
                                                &52819880671616779621183704306688233150514997832_wp
    real(wp),parameter ::   w4 = 0.10122853629037625915253135430996219011539409105168495705900369&
                                                &806474017876347078486028273930404500655815438933&
                                                &141326670771549403089234876787319730411360735846&
                                                &905332088240507319763065757292054679614357794675&
                                                &524923287300550259929540899466768105108107294683&
                                                &66466585774650346143712142008566866150514997832_wp

    f = h * ( w1*( me%fun(x-x1*h) + me%fun(x+x1*h)) + &
              w2*( me%fun(x-x2*h) + me%fun(x+x2*h)) + &
              w3*( me%fun(x-x3*h) + me%fun(x+x3*h)) + &
              w4*( me%fun(x-x4*h) + me%fun(x+x4*h)) )

    end function g8
!************************************************************************************

!************************************************************************************
!>
!  10-point method.
!
!### See also
!  * Coefficients from:
!    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

    function g10(me, x, h) result(f)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp), intent(in) :: x
    real(wp), intent(in) :: h
    real(wp)             :: f

    !> abscissae:
    real(wp),dimension(5),parameter ::  a = [   0.14887433898163121088482600112971998461756485942&
                                                &069169570798925351590361735566852137117762979946&
                                                &369123003116080525533882610289018186437654023167&
                                                &619699680909130507378277203710590709424758594227&
                                                &432498371771742473462169148529029429290031934666&
                                                &590824338380943550759968335702300050038372806343&
                                                &51_wp,&
                                                0.43339539412924719079926594316578416220007183765&
                                                &624649650270151314376698907770350122510275795011&
                                                &772122368293504099893794727422475772324920512677&
                                                &410328220862009523192709334620320113283203876915&
                                                &840634111498011298231414887874432043247664144215&
                                                &76788807708483879452488118549797039287926964254222_wp,&
                                                0.67940956829902440623432736511487357576929471183&
                                                &480946766481718895255857539507492461507857357048&
                                                &037949983390204739931506083674084257663009076827&
                                                &417182029235431978528469774097183691437120135529&
                                                &628377331531086791269325449548547293413247272116&
                                                &80274268486617121011712030227181051010718804444161_wp,&
                                                0.86506336668898451073209668842349304852754301496&
                                                &533045252195973184537475513805556135679072894604&
                                                &577069440463108641176516867830016149345356373927&
                                                &293968909500115713496898930516120724357604809009&
                                                &797259233179237955357392905958797769568324277022&
                                                &36942765911483643714816923781701572597289139322313_wp,&
                                                0.97390652851717172007796401208445205342826994669&
                                                &238211923121206669659520323463615962572356495626&
                                                &855625823304251877421121502216860143447777992054&
                                                &095872599424367044136957648812587991466331435107&
                                                &587371198778752105670674524353687136830338609093&
                                                &88311646653581707125686970668737259229449284383797_wp ]

    !> weights:
    real(wp),dimension(5),parameter ::  w = [   0.29552422471475287017389299465133832942104671702&
                                                &685360135430802975599593821715232927035659579375&
                                                &421672271716440125255838681849078955200582600193&
                                                &634249418696660956271864888416804323130506153586&
                                                &740908305127066386528748390174687472659751595445&
                                                &0775158914556548308329986393605934912382356670244_wp,&
                                                0.26926671930999635509122692156946935285975993846&
                                                &088379580056327624215343231917927676422663670925&
                                                &276075559581145036869830869292346938114524155646&
                                                &588466344237116560144322599601417290445280303444&
                                                &112979029770671425375348062846083992765750069116&
                                                &86749842814086288868533208042150419508881916391898_wp,&
                                                0.21908636251598204399553493422816319245877187052&
                                                &267708988095654363519991065295128124268399317720&
                                                &219278659121687281288763476662690806694756883092&
                                                &118433166566771052699153220775367726528266710278&
                                                &782468510102088321733200642734832547562506684158&
                                                &85349420711613410227291565477768928313300688702802_wp,&
                                                0.14945134915058059314577633965769733240255663966&
                                                &942736783547726875323865472663001094594726463473&
                                                &195191400575256104543633823445170674549760147137&
                                                &160119371095287981348288651187709535664396393337&
                                                &739399092016902046490838156187791575225783003434&
                                                &27785361756927642128792412282970150172590842897331_wp,&
                                                0.06667134430868813759356880989333179285786483432&
                                                &015814512869488161341206408408710177678550968505&
                                                &887782109005471452041933148750712625440376213930&
                                                &498731699404163449536370640018701124231550439352&
                                                &624245062983271819871864748056604411786208647844&
                                                &9236378557180717569208295026105115288152794421677_wp ]

    f = h * ( w(1)*(  me%fun(x-a(1)*h)   +  me%fun(x+a(1)*h) ) + &
              w(2)*(  me%fun(x-a(2)*h)   +  me%fun(x+a(2)*h) ) + &
              w(3)*(  me%fun(x-a(3)*h)   +  me%fun(x+a(3)*h) ) + &
              w(4)*(  me%fun(x-a(4)*h)   +  me%fun(x+a(4)*h) ) + &
              w(5)*(  me%fun(x-a(5)*h)   +  me%fun(x+a(5)*h) )   )

    end function g10
!************************************************************************************

!************************************************************************************
!>
!  12-point method.
!
!### See also
!  * Coefficients from:
!    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

    function g12(me, x, h) result(f)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp), intent(in) :: x
    real(wp), intent(in) :: h
    real(wp)             :: f

    !> abscissae:
    real(wp),dimension(6),parameter ::  a = [   0.12523340851146891547244136946385312998339691630&
                                                &544427321292175474846205624138968874286829846949&
                                                &135959410459879132051097315159969664463407959720&
                                                &578930281363427149751877364610797786290401085851&
                                                &749803458163536009061915338533985792224380950454&
                                                &5097342064247739686883799517760948964137522919201_wp,&
                                                0.36783149899818019375269153664371756125636014133&
                                                &540962131179987950408992951678787387873442850054&
                                                &657723463312639597714521513515217932743935324199&
                                                &163774275382871320389664162274303718284470963188&
                                                &934547884841822611461227526979609371629600504639&
                                                &62319787423676668046033025242558536362617894366679_wp,&
                                                0.58731795428661744729670241894053428036909851404&
                                                &805248151027087966734069937589526243571076498874&
                                                &820190960155999292889267723106959108867175142499&
                                                &189843704151965799654931521792486834699342245746&
                                                &542270556959107871794349154143635139191674285545&
                                                &96877940491139756923177447689738849120865435563147_wp,&
                                                0.76990267419430468703689383321281807598492575001&
                                                &893163766441906424911654310847122401642499922342&
                                                &191061761754045422185620704016285265354759491942&
                                                &035158754711514435184626896570143367857869960707&
                                                &068262822102488760216156789235759062543109515384&
                                                &10899341797549230707021382467596975621464477134163_wp,&
                                                0.90411725637047485667846586611909619253759670921&
                                                &329754655407576068123479572923579048696942782373&
                                                &326781186038289641042234889971981954299601063524&
                                                &901258268291998347354448614206140899100247009682&
                                                &576258221693446448698746167580757842398074380920&
                                                &64065954540171679180850205196702894963912359448494_wp,&
                                                0.98156063424671925069054909014928082296015519981&
                                                &373151046268212180779324431825398222525726789045&
                                                &223578555649237284127318524545703044707716708276&
                                                &967488752886112565550184482662910041202137201539&
                                                &996961235882788466302337187351583920530374414763&
                                                &9383170419389543470920618543180673569225988370568_wp]

    !> weights:
    real(wp),dimension(6),parameter ::  w = [   0.24914704581340278500056243604295121083046090256&
                                                &961883139535100311627942745728804303115680061804&
                                                &235306483347611787718583305851107360364968803964&
                                                &210377008542941502737221091728257019684301646591&
                                                &924021619820796255207324340857766137885796625403&
                                                &29347837170742904111565650371846972323325015720931_wp,&
                                                0.23349253653835480876084989892487805625940997219&
                                                &975487473052349782149200007941167528067902650856&
                                                &369046673875643970886883389854278840891609661975&
                                                &038847380753533248145179488750388812162792803042&
                                                &489598308782293577290791644231030018795306547073&
                                                &15375809270840669989018891281956753131165193423269_wp,&
                                                0.20316742672306592174906445580979837650651814727&
                                                &459014639859456579764563251047284379514439506460&
                                                &523243116042933686325996496137135190210132907910&
                                                &420189599423685656890245260738280276852445703846&
                                                &681240064758134063899875305215261728059344541572&
                                                &2327927963339557545261423500783899286052850767594_wp,&
                                                0.16007832854334622633465252954335907187201173049&
                                                &086417790989954415795422517329115068165655263705&
                                                &773052707487709681280262724376386088264904467503&
                                                &100243409511213679869020659979278560098046375913&
                                                &998387244938872586336059061667742863824552984448&
                                                &70458396283884610940466728874776625823124924247387_wp,&
                                                0.10693932599531843096025471819399622421457017347&
                                                &032488000512604210281899362749757654053731809631&
                                                &645741357635933314114416117033051696355084484800&
                                                &865232691960050114390447649204829355515358576079&
                                                &107052492180710337954704248957128309309678064675&
                                                &98358517298903536137451828089012822811396588037254_wp,&
                                                0.04717533638651182719461596148501706031702907399&
                                                &484708956050534700380972115203871067082590707541&
                                                &453609661610167559673857966748040823913299673846&
                                                &365109909808575797967885849598965975687054894525&
                                                &799700269519193179311245399071070942125321236826&
                                                &63180160342232703368882666374567833050364187887189_wp]

    f = h * ( w(1)*(  me%fun(x-a(1)*h)   +  me%fun(x+a(1)*h) ) + &
              w(2)*(  me%fun(x-a(2)*h)   +  me%fun(x+a(2)*h) ) + &
              w(3)*(  me%fun(x-a(3)*h)   +  me%fun(x+a(3)*h) ) + &
              w(4)*(  me%fun(x-a(4)*h)   +  me%fun(x+a(4)*h) ) + &
              w(5)*(  me%fun(x-a(5)*h)   +  me%fun(x+a(5)*h) ) + &
              w(6)*(  me%fun(x-a(6)*h)   +  me%fun(x+a(6)*h) ) )

    end function g12
!************************************************************************************

!************************************************************************************
!>
!  14-point method.
!
!### See also
!  * Coefficients from:
!    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

    function g14(me, x, h) result(f)

    implicit none

    class(integration_class_1d),intent(inout)  :: me
    real(wp), intent(in) :: x
    real(wp), intent(in) :: h
    real(wp)             :: f

    !> abscissae:
    real(wp),dimension(7),parameter ::  a = [  0.10805494870734366206624465021983474761195160547&
                                               &4237557040821061308013529011730007130100688176689&
                                               &3672374502026424466474638099232632258191427567218&
                                               &1973150409752806137273842265069487944308775321508&
                                               &8445556391329819060204836416480024319739665907101&
                                               &2506161702814425014635643221773541001328892761_wp,&
                                               &0.31911236892788976043567182416847546683426120353&
                                               &3843956596650187257333440512792783164933705421346&
                                               &4131802793151826090394496145640578710017716508863&
                                               &2222396245608012120993128542172348808287716458637&
                                               &8479374239121304478425121768114783511643536777896&
                                               &2949997448460558214759676525644841351801594858_wp,&
                                               &0.51524863635815409196529071855118866230888528256&
                                               &9306036951504769092784951832055660452072020350772&
                                               &8923922907932905090138695274035571340047593918260&
                                               &5653057211011637652073200342580823038532041784020&
                                               &3436173906624491224801618641571038235567674745455&
                                               &3979637438627635490786064892912451481973721288_wp,&
                                               &0.68729290481168547014801980301933413753840121274&
                                               &7170675619266488628184896183133256947373070505211&
                                               &8384106603630216790054729627432715418501010682124&
                                               &6881727389082952662885443589912839338608106959371&
                                               &4595904926885388784713769175169784875289055161406&
                                               &7877996475717650653147982694804026342351254071_wp,&
                                               &0.82720131506976499318979474265039496103970110147&
                                               &5081181560709054241479830810028873570426390137889&
                                               &5453991241406273986535333275661226737816179582645&
                                               &1069907936808669317564778014567859855078251147291&
                                               &5830426696849656086721489336979443959282673643228&
                                               &6425172143208924251106624044295037127737490111_wp,&
                                               &0.92843488366357351733639113937787426447703921040&
                                               &9837618717962447482131093544359853111413905683657&
                                               &5176363551261559882603607008578010786539258018984&
                                               &5400440650494157888098179531161147719130825235345&
                                               &8596605653673043686690855550898698329741248613224&
                                               &5749388483890945436457404705549484348178721002_wp,&
                                               &0.98628380869681233884159726670405280167609140723&
                                               &9225881644070811777749554132491637910646239665151&
                                               &7527602612562941358578689852603067447974494119727&
                                               &0324710898207170072955675048180261687970555989447&
                                               &5396929426197069500447181272675429908986256542893&
                                               &3676463914802477677291745002965827767360741735_wp ]

    !> weights:
    real(wp),dimension(7),parameter ::  w = [   0.21526385346315779019587644331626003527499755805&
                                               &4128800219776392543618787353994604001024441410819&
                                               &5782372566723324367709929481659764649301890356019&
                                               &0805098142804175780269156508228762641736544919294&
                                               &6281203662033345376460522564310634412912654698349&
                                               &487266562730897512393716549425155133887783267_wp,&
                                               &0.20519846372129560396592406566121805571033906130&
                                               &9419451716897290283367144825249720339431839991890&
                                               &8957243692694424494287284534856133850644865918702&
                                               &3021403166714178733299347482783913811132568481282&
                                               &5439676020905052976535424973123755325146919285189&
                                               &8072394707049964721031773292256965337005468577_wp,&
                                               &0.18553839747793781374171659012515703624892260293&
                                               &7331659020034925069098350263525444425552731146712&
                                               &2229825611215057289188990778964974252160895085525&
                                               &2415283643607286404060027232379697141385075345609&
                                               &3331227890449938852384485366393922617921879824760&
                                               &6150274514935557012909889503067356410067833406_wp,&
                                               &0.15720316715819353456960193862384215660566803733&
                                               &7323374969317043874768176369608298513958093362418&
                                               &0762768531519990811885018854374920646576267489242&
                                               &9103726460198700102219564745910784232280561068611&
                                               &6907713218466935160138377442838502265889923868443&
                                               &9084685022864905124096570215866733146092008329_wp,&
                                               &0.12151857068790318468941480907247662595666934569&
                                               &0074672291075392543159743892526492318819906270375&
                                               &0071489155506530592569942811574313408868548096421&
                                               &2571445460802891854106154207862005646754562932960&
                                               &2540610239636717985405900755004972904989241013019&
                                               &1072357341821083329663867464821867539341968434_wp,&
                                               &0.08015808715976020980563327706285430958369778539&
                                               &4594765201399065489571474457287169863536190819137&
                                               &7559686225015908038847487953091382572604434376755&
                                               &1198447409477973877237005366105771785226539545491&
                                               &2313554662497115946457665357652160093748935412771&
                                               &0937535198838649279475628473516378736712929573_wp,&
                                               &0.03511946033175186303183287613819178061970560927&
                                               &7127276581499890196416322837808270537676796998646&
                                               &4636614217324764405511345585478510619843098677334&
                                               &0884595716394793248808744456729064741484147706750&
                                               &3186014306010893702617623540676052379390445897465&
                                               &9810087587180865408885105556219147609526200925_wp ]

    f = h * ( w(1)*(  me%fun(x-a(1)*h)   +  me%fun(x+a(1)*h) ) + &
              w(2)*(  me%fun(x-a(2)*h)   +  me%fun(x+a(2)*h) ) + &
              w(3)*(  me%fun(x-a(3)*h)   +  me%fun(x+a(3)*h) ) + &
              w(4)*(  me%fun(x-a(4)*h)   +  me%fun(x+a(4)*h) ) + &
              w(5)*(  me%fun(x-a(5)*h)   +  me%fun(x+a(5)*h) ) + &
              w(6)*(  me%fun(x-a(6)*h)   +  me%fun(x+a(6)*h) ) + &
              w(7)*(  me%fun(x-a(7)*h)   +  me%fun(x+a(7)*h) ) )

    end function g14
!************************************************************************************

!*******************************************************************************************************
    end module quadrature_module
!*******************************************************************************************************