!------------------------------------------------------------------------------
!        Colab+Atlantic, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : kernelNonlinearSolverScalar_mod
! URL           : http://www.mohid.com
! AFFILIATION   : Colab+Atlantic, Marine Modelling Group
! DATE          : May 2026
!> @author
!> Mohsen Shabani Email: shabani.mohsen@outlook.com
!
! DESCRIPTION:
!> General nonlinear solvers for scalar equations f(x)=0
!> Object-oriented version with type-bound procedures
!
!------------------------------------------------------------------------------

	module kernelNonlinearSolverScalar_mod

	use common_modules, only: prec
	implicit none
	private

	public :: kernelNonlinearSolverScalar_class

	abstract interface
	 function fun_scalar(x) result(fx)
	   import :: prec
	   real(prec), intent(in) :: x
	   real(prec) :: fx
	 end function fun_scalar
	end interface

	type :: kernelNonlinearSolverScalar_class
	contains
	 procedure :: NewtonSolver
	 procedure :: SecantSolver
	 procedure, private :: derivative_function
	end type kernelNonlinearSolverScalar_class

	contains

	  !--------------------------------------------------------------------------
	  ! Numerical derivative using central difference
	  !--------------------------------------------------------------------------
	  function derivative_function(self, f, x) result(dfdx)

		class(kernelNonlinearSolverScalar_class), intent(in) :: self
		procedure(fun_scalar) :: f
		real(prec), intent(in) :: x
		real(prec) :: dfdx
		real(prec) :: h

		h = epsilon(1.0_prec)**(1.0_prec/3.0_prec) * max(1.0_prec, abs(x))

		dfdx = ( f(x + h) - f(x - h) ) / (2.0_prec*h)

	  end function derivative_function

	  !--------------------------------------------------------------------------
	  ! Newton method with finite-difference derivative
	  !--------------------------------------------------------------------------
	  subroutine NewtonSolver(self, f, x0, root, converged, tol, maxit)

		class(kernelNonlinearSolverScalar_class), intent(in) :: self
		procedure(fun_scalar) :: f
		real(prec), intent(in) :: x0
		real(prec), intent(out) :: root
		logical, intent(out) :: converged
		real(prec), intent(in), optional :: tol
		integer, intent(in), optional :: maxit

		real(prec) :: x, xnew, fx, fxnew, dfx, mytol
		integer :: i, nmax

		mytol = 1.0e-10_prec
		if (present(tol)) mytol = tol

		nmax = 100
		if (present(maxit)) nmax = maxit

		x = x0
		converged = .false.

		do i = 1, nmax

		   fx  = f(x)
		   dfx = self%derivative_function(f, x)

		   if (abs(dfx) < 1.0e-14_prec) then
			  root = x
			  return
		   end if

		   xnew  = x - fx/dfx
		   fxnew = f(xnew)

		   if (abs(xnew - x) < mytol .or. abs(fxnew) < mytol) then
			  root = xnew
			  converged = .true.
			  return
		   end if

		   x = xnew

		end do

		root = x

	  end subroutine NewtonSolver

	  !--------------------------------------------------------------------------
	  ! Secant method
	  !--------------------------------------------------------------------------
	  subroutine SecantSolver(self, f, x0, x1, root, converged, tol, maxit)

		class(kernelNonlinearSolverScalar_class), intent(in) :: self
		procedure(fun_scalar) :: f
		real(prec), intent(in) :: x0, x1
		real(prec), intent(out) :: root
		logical, intent(out) :: converged
		real(prec), intent(in), optional :: tol
		integer, intent(in), optional :: maxit

		real(prec) :: xm1, x, xnew, fm1, fx, fxnew, mytol
		integer :: i, nmax

		mytol = 1.0e-10_prec
		if (present(tol)) mytol = tol

		nmax = 100
		if (present(maxit)) nmax = maxit

		xm1 = x0
		x   = x1
		converged = .false.

		do i = 1, nmax

		   fm1 = f(xm1)
		   fx  = f(x)

		   if (abs(fx - fm1) < 1.0e-14_prec) then
			  root = x
			  return
		   end if

		   xnew  = x - fx * (x - xm1) / (fx - fm1)
		   fxnew = f(xnew)

		   if (abs(xnew - x) < mytol .or. abs(fxnew) < mytol) then
			  root = xnew
			  converged = .true.
			  return
		   end if

		   xm1 = x
		   x   = xnew

		end do

		root = x

	  end subroutine SecantSolver

	end module kernelNonlinearSolverScalar_mod