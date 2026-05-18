!------------------------------------------------------------------------------
!        Colab+Atlantic, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : kernelNonlinearSolverVector_mod
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
	module kernelNonlinearSolverVector_mod

	  use common_modules, only: prec
	  implicit none
	  private

	  public :: kernelNonlinearSolverVector_class

	  abstract interface
		 function fun_array(x) result(fx)
		   import :: prec
		   real(prec), dimension(:), intent(in) :: x
		   real(prec), dimension(size(x)) :: fx
		 end function fun_array
	  end interface

	  type :: kernelNonlinearSolverVector_class
	  contains
		 procedure :: NewtonSolver
		 procedure :: SecantSolver
		 procedure, private :: derivative_function
	  end type kernelNonlinearSolverVector_class

	contains
	  !--------------------------------------------------------------------------
	  ! Numerical derivative using central difference
	  !--------------------------------------------------------------------------
	  function derivative_function(self, f, x) result(dfdx)
		class(kernelNonlinearSolverVector_class), intent(in) :: self
		procedure(fun_array) :: f
		real(prec), dimension(:), intent(in) :: x
		real(prec), dimension(size(x)) :: dfdx

		real(prec), dimension(size(x)) :: xp, xm, fp, fm, h
		integer :: i

		h = epsilon(1.0_prec)**(1.0_prec/3.0_prec) * max(1.0_prec, abs(x))
		dfdx = 0.0_prec

		do i = 1, size(x)
		   xp = x
		   xm = x

		   xp(i) = xp(i) + h(i)
		   xm(i) = xm(i) - h(i)

		   fp = f(xp)
		   fm = f(xm)

		   ! diagonal derivative: dF_i / dx_i
		   dfdx(i) = (fp(i) - fm(i)) / (2.0_prec*h(i))
		end do

	  end function derivative_function

	  !--------------------------------------------------------------------------
	  ! Newton method with finite-difference derivative
	  !--------------------------------------------------------------------------

	  subroutine NewtonSolver(self, f, x0, root, converged, tol, maxit)
		class(kernelNonlinearSolverVector_class), intent(in) :: self
		procedure(fun_array) :: f
		real(prec), dimension(:), intent(in) :: x0
		real(prec), dimension(:), intent(out) :: root
		logical, dimension(:), intent(out) :: converged
		real(prec), intent(in), optional :: tol
		integer, intent(in), optional :: maxit

		real(prec), dimension(size(x0)) :: x, xnew, fx, fxnew, dfx
		real(prec) :: mytol
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

		   where (abs(dfx) >= 1.0e-14_prec)
			  xnew = x - fx/dfx
		   elsewhere
			  xnew = x
		   end where

		   fxnew = f(xnew)

		   where (abs(xnew - x) < mytol .or. abs(fxnew) < mytol)
			  converged = .true.
		   end where

		   x = xnew

		   if (all(converged)) exit
		end do

		root = x
	  end subroutine NewtonSolver

	  !--------------------------------------------------------------------------
	  ! Secant method
	  !--------------------------------------------------------------------------
	  subroutine SecantSolver(self, f, x0, x1, root, converged, tol, maxit)
		class(kernelNonlinearSolverVector_class), intent(in) :: self
		procedure(fun_array) :: f
		real(prec), dimension(:), intent(in) :: x0, x1
		real(prec), dimension(:), intent(out) :: root
		logical, dimension(:), intent(out) :: converged
		real(prec), intent(in), optional :: tol
		integer, intent(in), optional :: maxit

		real(prec), dimension(size(x0)) :: xm1, x, xnew, fm1, fx, fxnew
		real(prec) :: mytol
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

		   where (abs(fx - fm1) >= 1.0e-14_prec)
			  xnew = x - fx * (x - xm1) / (fx - fm1)
		   elsewhere
			  xnew = x
		   end where

		   fxnew = f(xnew)

		   where (abs(xnew - x) < mytol .or. abs(fxnew) < mytol)
			  converged = .true.
		   end where

		   xm1 = x
		   x   = xnew

		   if (all(converged)) exit
		end do

		root = x
	  end subroutine SecantSolver

	end module kernelNonlinearSolverVector_mod