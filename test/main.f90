module test_stats
  use stats
  use stdlib_kinds, only: int32, dp
  use stdlib_stats, only: mean, var
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_stats

  contains

    !> Collect all exported unit tests
    subroutine collect_stats(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("mean", test_mean), new_unittest("variance", test_var)]
    end subroutine collect_stats

    !> Check accumulator mean
    subroutine test_mean(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type (acc)     :: t
      real(dp)       :: x(1000), avg1, avg2

      call random_init(.false., .false.)
      call random_number(x)
      call t%push(x)

      avg1 = t%mean()
      avg2 = mean(x, 1)

      ! call check(error, avg2, avg1, thr=real(1.e-8, dp))
      call check(error, avg2, avg1, thr=epsilon(avg2)*100)
      ! call check(error, avg2, avg1)
      if (allocated(error)) return

    end subroutine test_mean

    !> Check accumulator variance
    subroutine test_var(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type (acc)     :: t
      real(dp)       :: x(1000), var1, var2

      call random_init(.false., .false.)
      call random_number(x)
      call t%push(x)

      var1 = t%vars()
      var2 = var(x, 1)

      ! call check(error, avg2, avg1, thr=real(1.e-8, dp))
      call check(error, var2, var1, thr=epsilon(var2)*100)
      ! call check(error, avg2, avg1)
      if (allocated(error)) return

    end subroutine test_var

end module test_stats

program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite
  use test_stats, only : collect_stats
  implicit none
  integer :: stat

  stat = 0
  call run_testsuite(collect_stats, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester
