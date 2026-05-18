!------------------------------------------------------------------------------
!        Colab+Atlantic, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : kernelNonlinearSolverVectorCoupling_mod
! URL           : http://www.mohid.com
! AFFILIATION   : Colab+Atlantic, Marine Modelling Group
! DATE          : May 2026
!> @author
!> Mohsen Shabani Email: shabani.mohsen@outlook.com
!
! DESCRIPTION:
!> General nonlinear solvers for coupled vector equations F(x)=0
!> Robust object-oriented version with type-bound procedures
!
!------------------------------------------------------------------------------
module kernelNonlinearSolverVectorCoupling_mod

   use common_modules, only: prec
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public :: kernelNonlinearSolverVectorCoupling_class

   abstract interface
      function fun_array(x) result(fx)
         import :: prec
         real(prec), intent(in) :: x(:)
         real(prec)             :: fx(size(x))
      end function fun_array
   end interface

   type :: kernelNonlinearSolverVectorCoupling_class
   contains
      procedure :: NewtonSolver
      procedure :: SecantSolver
      procedure, private :: derivative_function
      procedure, private :: solve_linear_system
      procedure, private :: regularized_solve_linear_system
      procedure, private :: outer_product
   end type kernelNonlinearSolverVectorCoupling_class

contains

   !--------------------------------------------------------------------------
   ! Infinity norm of a vector
   !--------------------------------------------------------------------------
   pure function vec_norm_inf(v) result(nrm)
      real(prec), intent(in) :: v(:)
      real(prec) :: nrm

      if (size(v) == 0) then
         nrm = 0.0_prec
      else
         nrm = maxval(abs(v))
      end if
   end function vec_norm_inf

   !--------------------------------------------------------------------------
   ! Full numerical Jacobian using robust finite differences
   ! J(i,j) = dF_i / dx_j
   !
   ! Strategy:
   !   1) Try central difference
   !   2) If that fails, reduce step
   !   3) If still failing, try one-sided difference
   !   4) If all fail, return NaN matrix so caller can reject Jacobian
   !--------------------------------------------------------------------------
   function derivative_function(self, f, x) result(J)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      procedure(fun_array) :: f
      real(prec), intent(in) :: x(:)
      real(prec) :: J(size(x), size(x))

      real(prec) :: xp(size(x)), xm(size(x))
      real(prec) :: fp(size(x)), fm(size(x)), f0(size(x))
      real(prec) :: h(size(x)), hcol
      real(prec) :: nanv
      integer :: icol, itry
      logical :: got_col

      nanv = ieee_value(1.0_prec, ieee_quiet_nan)
      J = 0.0_prec

      if (.not. all(ieee_is_finite(x))) then
         J = nanv
         return
      end if

      f0 = f(x)
      if (.not. all(ieee_is_finite(f0))) then
         J = nanv
         return
      end if

      ! bounded perturbation size
      h = 1.0e-6_prec * max(1.0_prec, abs(x))
      h = min(h, 1.0e-4_prec)
      h = max(h, 1.0e-8_prec)

      do icol = 1, size(x)

         got_col = .false.
         hcol = h(icol)

         !---------------------------------------
         ! Try central difference with shrinking h
         !---------------------------------------
         do itry = 1, 5
            xp = x
            xm = x

            xp(icol) = xp(icol) + hcol
            xm(icol) = xm(icol) - hcol

            if (.not. all(ieee_is_finite(xp))) then
               hcol = 0.25_prec * hcol
               cycle
            end if

            if (.not. all(ieee_is_finite(xm))) then
               hcol = 0.25_prec * hcol
               cycle
            end if

            fp = f(xp)
            fm = f(xm)

            if (all(ieee_is_finite(fp)) .and. all(ieee_is_finite(fm))) then
               J(:,icol) = (fp - fm) / (2.0_prec * hcol)
               got_col = .true.
               exit
            end if

            hcol = 0.25_prec * hcol
         end do

         !---------------------------------------
         ! Fallback: one-sided forward difference
         !---------------------------------------
         if (.not. got_col) then
            hcol = h(icol)

            do itry = 1, 5
               xp = x
               xp(icol) = xp(icol) + hcol

               if (.not. all(ieee_is_finite(xp))) then
                  hcol = 0.25_prec * hcol
                  cycle
               end if

               fp = f(xp)

               if (all(ieee_is_finite(fp))) then
                  J(:,icol) = (fp - f0) / hcol
                  got_col = .true.
                  exit
               end if

               hcol = 0.25_prec * hcol
            end do
         end if

         if (.not. got_col) then
            J = nanv
            return
         end if

         if (.not. all(ieee_is_finite(J(:,icol)))) then
            J = nanv
            return
         end if

      end do

   end function derivative_function

   !--------------------------------------------------------------------------
   ! Outer product: A = u * v^T
   !--------------------------------------------------------------------------
   function outer_product(self, u, v) result(A)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      real(prec), intent(in) :: u(:), v(:)
      real(prec) :: A(size(u), size(v))
      integer :: j

      A = 0.0_prec
      do j = 1, size(v)
         A(:,j) = u * v(j)
      end do

   end function outer_product

   !--------------------------------------------------------------------------
   ! Solve linear system A*x = b using Gaussian elimination with partial pivoting
   !--------------------------------------------------------------------------
   subroutine solve_linear_system(self, A_in, b_in, x, ok)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      real(prec), intent(in)  :: A_in(:,:), b_in(:)
      real(prec), intent(out) :: x(size(b_in))
      logical, intent(out)    :: ok

      real(prec) :: A(size(A_in,1), size(A_in,2))
      real(prec) :: b(size(b_in))
      real(prec) :: factor, pivot_abs, scaleA
      real(prec) :: tmp
      real(prec) :: row_tmp(size(A_in,2))
      integer :: n
      integer :: i, k, p

      n = size(b_in)
      ok = .false.
      x  = 0.0_prec

      if (size(A_in,1) /= n .or. size(A_in,2) /= n) return
      if (.not. all(ieee_is_finite(A_in))) return
      if (.not. all(ieee_is_finite(b_in))) return

      A = A_in
      b = b_in

      scaleA = max(1.0_prec, maxval(abs(A)))

      do k = 1, n-1

         p = k
         pivot_abs = abs(A(k,k))

         do i = k+1, n
            if (abs(A(i,k)) > pivot_abs) then
               p = i
               pivot_abs = abs(A(i,k))
            end if
         end do

         if (pivot_abs < 1.0e-12_prec * scaleA) return

         if (p /= k) then
            row_tmp  = A(k,:)
            A(k,:)   = A(p,:)
            A(p,:)   = row_tmp

            tmp = b(k)
            b(k) = b(p)
            b(p) = tmp
         end if

         do i = k+1, n
            factor = A(i,k) / A(k,k)
            A(i,k:n) = A(i,k:n) - factor * A(k,k:n)
            b(i)     = b(i)     - factor * b(k)
         end do

         if (.not. all(ieee_is_finite(A))) return
         if (.not. all(ieee_is_finite(b))) return

      end do

      if (abs(A(n,n)) < 1.0e-12_prec * scaleA) return

      x(n) = b(n) / A(n,n)

      do i = n-1, 1, -1
         x(i) = (b(i) - sum(A(i,i+1:n) * x(i+1:n))) / A(i,i)
      end do

      if (.not. all(ieee_is_finite(x))) return

      ok = .true.

   end subroutine solve_linear_system

   !--------------------------------------------------------------------------
   ! Solve linear system with small diagonal regularization if needed
   !--------------------------------------------------------------------------
   subroutine regularized_solve_linear_system(self, A_in, b_in, x, ok)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      real(prec), intent(in)  :: A_in(:,:), b_in(:)
      real(prec), intent(out) :: x(size(b_in))
      logical, intent(out)    :: ok

      real(prec) :: A(size(A_in,1), size(A_in,2))
      real(prec) :: scaleA, lambda
      integer :: i, itry

      ok = .false.
      x  = 0.0_prec

      if (.not. all(ieee_is_finite(A_in))) return
      if (.not. all(ieee_is_finite(b_in))) return

      scaleA = max(1.0_prec, maxval(abs(A_in)))

      do itry = 0, 6
         A = A_in

         if (itry > 0) then
            lambda = (10.0_prec**(itry-1)) * 1.0e-10_prec * scaleA
            do i = 1, size(A,1)
               A(i,i) = A(i,i) + lambda
            end do
         end if

         call self%solve_linear_system(A, b_in, x, ok)

         if (ok) then
            if (all(ieee_is_finite(x))) return
            ok = .false.
         end if
      end do

   end subroutine regularized_solve_linear_system

   !--------------------------------------------------------------------------
   ! Coupled Newton method with robust finite-difference Jacobian
   !
   ! NOTE:
   ! converged(:) is kept as an array only for backward compatibility.
   ! All entries are set equal to the same global convergence flag.
   !--------------------------------------------------------------------------
   subroutine NewtonSolver(self, f, x0, root, converged, tol, maxit)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      procedure(fun_array) :: f
      real(prec), intent(in)  :: x0(:)
      real(prec), intent(out) :: root(:)
      logical, intent(out)    :: converged(:)
      real(prec), intent(in), optional :: tol
      integer, intent(in), optional :: maxit

      real(prec) :: x(size(x0)), xnew(size(x0)), xtrial(size(x0))
      real(prec) :: fx(size(x0)), fxnew(size(x0)), fxtrial(size(x0))
      real(prec) :: dx(size(x0))
      real(prec) :: J(size(x0), size(x0))
      real(prec) :: mytol, alpha, fnorm, fnorm_trial, max_step
      integer :: it, nmax
      logical :: ok, accepted

      mytol = 1.0e-10_prec
      if (present(tol)) mytol = tol

      nmax = 100
      if (present(maxit)) nmax = maxit

      max_step = 10.0_prec

      x = x0
      root = x0
      ok = .false.
      converged = .false.

      if (.not. all(ieee_is_finite(x))) return

      fx = f(x)
      if (.not. all(ieee_is_finite(fx))) return

      do it = 1, nmax

         fnorm = vec_norm_inf(fx)

         if (fnorm < mytol) then
            ok = .true.
            exit
         end if

         J = self%derivative_function(f, x)
         if (.not. all(ieee_is_finite(J))) then
            write(*,*) 'BAD J IN NEWTON'
            root = x
            converged = .false.
            return
         end if

         call self%regularized_solve_linear_system(J, -fx, dx, ok)
         if (.not. ok) then
            root = x
            converged = .false.
            return
         end if

         if (.not. all(ieee_is_finite(dx))) then
            write(*,*) 'BAD DX IN NEWTON', dx
            root = x
            converged = .false.
            return
         end if

         ! limit excessively large steps
         if (vec_norm_inf(dx) > max_step) then
            dx = dx * (max_step / vec_norm_inf(dx))
         end if

         accepted = .false.
         alpha = 1.0_prec

         do while (alpha >= 1.0e-6_prec)

            xtrial = x + alpha * dx

            if (all(ieee_is_finite(xtrial))) then
               fxtrial = f(xtrial)

               if (all(ieee_is_finite(fxtrial))) then
                  fnorm_trial = vec_norm_inf(fxtrial)

                  if (fnorm_trial <= fnorm) then
                     xnew = xtrial
                     fxnew = fxtrial
                     accepted = .true.
                     exit
                  end if
               end if
            end if

            alpha = 0.5_prec * alpha
         end do

         if (.not. accepted) then
            root = x
            converged = .false.
            return
         end if

         if (.not. all(ieee_is_finite(xnew))) then
            write(*,*) 'BAD XNEW IN NEWTON', x, dx, xnew
            root = x
            converged = .false.
            return
         end if

         x = xnew
         fx = fxnew

         if (vec_norm_inf(alpha*dx) < mytol .or. vec_norm_inf(fx) < mytol) then
            ok = .true.
            exit
         end if

      end do

      root = x
      converged = ok

   end subroutine NewtonSolver

   !--------------------------------------------------------------------------
   ! Coupled "SecantSolver" implemented as robust Broyden's method
   !
   ! Keeps the same external interface as your old SecantSolver.
   !--------------------------------------------------------------------------
   subroutine SecantSolver(self, f, x0, x1, root, converged, tol, maxit)
      class(kernelNonlinearSolverVectorCoupling_class), intent(in) :: self
      procedure(fun_array) :: f
      real(prec), intent(in)  :: x0(:), x1(:)
      real(prec), intent(out) :: root(:)
      logical, intent(out)    :: converged(:)
      real(prec), intent(in), optional :: tol
      integer, intent(in), optional :: maxit

      real(prec) :: x(size(x0)), xnew(size(x0)), xtrial(size(x0))
      real(prec) :: fx(size(x0)), fxnew(size(x0)), fxtrial(size(x0)), fm1(size(x0))
      real(prec) :: dx(size(x0))
      real(prec) :: s(size(x0)), y(size(x0)), Js(size(x0))
      real(prec) :: J(size(x0), size(x0))
      real(prec) :: denom, mytol, alpha, fnorm, fnorm_trial, max_step
      integer :: it, nmax
      logical :: ok, accepted

      mytol = 1.0e-10_prec
      if (present(tol)) mytol = tol

      nmax = 100
      if (present(maxit)) nmax = maxit

      max_step = 10.0_prec

      root = x1
      converged = .false.
      ok = .false.

      if (.not. all(ieee_is_finite(x0))) return
      if (.not. all(ieee_is_finite(x1))) return

      x   = x1
      fx  = f(x)
      fm1 = f(x0)

      if (.not. all(ieee_is_finite(fx)))  return
      if (.not. all(ieee_is_finite(fm1))) return

      J = self%derivative_function(f, x)
      if (.not. all(ieee_is_finite(J))) return

      s = x1 - x0
      y = fx - fm1
      denom = dot_product(s, s)

      if (denom > 1.0e-20_prec) then
         Js = matmul(J, s)
         if (all(ieee_is_finite(Js))) then
            J = J + self%outer_product(y - Js, s) / denom
         end if
      end if

      do it = 1, nmax

         fnorm = vec_norm_inf(fx)

         if (fnorm < mytol) then
            ok = .true.
            exit
         end if

         call self%regularized_solve_linear_system(J, -fx, dx, ok)

         if (.not. ok) then
            J = self%derivative_function(f, x)
            if (.not. all(ieee_is_finite(J))) then
               root = x
               converged = .false.
               return
            end if

            call self%regularized_solve_linear_system(J, -fx, dx, ok)
            if (.not. ok) then
               root = x
               converged = .false.
               return
            end if
         end if

         if (.not. all(ieee_is_finite(dx))) then
            root = x
            converged = .false.
            return
         end if

         if (vec_norm_inf(dx) > max_step) then
            dx = dx * (max_step / vec_norm_inf(dx))
         end if

         accepted = .false.
         alpha = 1.0_prec

         do while (alpha >= 1.0e-6_prec)

            xtrial = x + alpha * dx

            if (all(ieee_is_finite(xtrial))) then
               fxtrial = f(xtrial)

               if (all(ieee_is_finite(fxtrial))) then
                  fnorm_trial = vec_norm_inf(fxtrial)

                  if (fnorm_trial <= fnorm) then
                     xnew = xtrial
                     fxnew = fxtrial
                     accepted = .true.
                     exit
                  end if
               end if
            end if

            alpha = 0.5_prec * alpha
         end do

         if (.not. accepted) then
            root = x
            converged = .false.
            return
         end if

         if (vec_norm_inf(alpha*dx) < mytol .or. vec_norm_inf(fxnew) < mytol) then
            x = xnew
            ok = .true.
            exit
         end if

         s = xnew - x
         y = fxnew - fx
         denom = dot_product(s, s)

         if (denom > 1.0e-20_prec) then
            Js = matmul(J, s)
            if (all(ieee_is_finite(Js))) then
               J = J + self%outer_product(y - Js, s) / denom
            else
               J = self%derivative_function(f, xnew)
            end if
         else
            J = self%derivative_function(f, xnew)
         end if

         if (.not. all(ieee_is_finite(J))) then
            J = self%derivative_function(f, xnew)
            if (.not. all(ieee_is_finite(J))) then
               root = x
               converged = .false.
               return
            end if
         end if

         x  = xnew
         fx = fxnew

      end do

      root = x
      converged = ok

   end subroutine SecantSolver

end module kernelNonlinearSolverVectorCoupling_mod