module stats
  !> real64 statistics accumulator
  !> based on the nim std/stats module: https://nim-lang.org/docs/stats.html

  use stdlib_kinds, only: int32, dp

  implicit none
  private

  type, public                     :: acc
    integer(int32)                 :: n    = 0
    real(dp)                       :: min  = 0
    real(dp)                       :: max  = 0
    real(dp)                       :: sum  = 0
    real(dp)                       :: mom1 = 0
    real(dp)                       :: mom2 = 0
    real(dp)                       :: mom3 = 0
    real(dp)                       :: mom4 = 0

    contains

      procedure                    :: push_real, push_reala
      generic                      :: push => push_real, push_reala
      procedure                    :: clear, mean, vars, varp, stds, stdp, skews, skewp, kurts, kurtp, show

  end type

  contains

    subroutine push_real(s, x)
      !> pushes a real value x into the accumulator.
      class(acc), intent(inout)  :: s
      real(dp), intent(in)       :: x
      real(dp)                   :: m, p, delta, delta_n, delta_n2, term1
      if (s%n == 0) then
        s%min = x
        s%max = x
      else
        if (s%min > x) then
          s%min = x
        end if
        if (s%max < x) then
          s%max = x
        end if
      end if
      ! see Knuth TAOCP vol 2, 3rd edition, page 232
      s%n      = s%n+1
      s%sum    = s%sum+x
      m        = real(s%n, dp)
      p        = real(s%n-1, dp)
      delta    = x-s%mom1
      delta_n  = delta/m
      delta_n2 = delta_n*delta_n
      term1    = delta*delta_n*p
      s%mom4   = s%mom4+term1*delta_n2*(m*m-3*m+3)+6*delta_n2*s%mom2-4*delta_n*s%mom3
      s%mom3   = s%mom3+term1*delta_n*(m-2)-3*delta_n*s%mom2
      s%mom2   = s%mom2+term1
      s%mom1   = s%mom1+delta_n
    end subroutine push_real

    subroutine push_reala(s, a)
      !> pushes a real array into the accumulator.
      class(acc), intent(inout)  :: s
      real(dp), intent(in)       :: a(:)
      integer(int32)             :: i, n
      n = size(a)
      do i = 1, n
        call s%push(a(i))
      end do
    end subroutine push_reala

    subroutine clear(s)
      !> rets accumulator.
      class(acc), intent(inout)    :: s
      s%n    = 0
      s%min  = 0
      s%max  = 0
      s%sum  = 0
      s%mom1 = 0
      s%mom2 = 0
      s%mom3 = 0
      s%mom4 = 0
    end subroutine clear

    function mean(s) result(r)
      !> computes the current mean of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = s%mom1
    end function mean

    function varp(s) result(r)
      !> computes the current population variance of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      if (s%n > 1) then
        r = s%mom2/real(s%n, dp)
      else
        r = 0
      end if
    end function varp

    function vars(s) result(r)
      !> computes the current sample variance of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      if (s%n > 1) then
        r = s%mom2/real(s%n-1, dp)
      else
        r = 0
      end if
    end function vars

    function stdp(s) result(r)
      !> computes the current population standard deviation of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = sqrt(varp(s))
    end function stdp

    function stds(s) result(r)
      !> computes the current sample standard deviation of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = sqrt(vars(s))
    end function stds

    function skewp(s) result(r)
      !> computes the current population skewness of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = sqrt(real(s%n, dp))*s%mom3/s%mom2**1.5
    end function skewp

    function skews(s) result(r)
      !> computes the current sample skewness of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r, s2
      s2 = skewp(s)
      if (s%n > 2) then
        r = sqrt(real(s%n*(s%n-1), dp))*s2/real(s%n-2, dp)
      else
        r = 0
      end if
    end function skews

    function kurtp(s) result(r)
      !> computes the current population kurtosis of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = real(s%n, dp)*s%mom4/(s%mom2*s%mom2)-3.0_dp
    end function kurtp

    function kurts(s) result(r)
      !> computes the current sample kurtosis of the accumulator.
      class(acc), intent(in)       :: s
      real(dp)                     :: r
      r = real(s%n-1, dp)/real((s%n-2)*(s%n-3), dp)*(real(s%n+1,dp)*kurtp(s)+6)
    end function kurts

    subroutine show(s)
      !> print stats.
      class(acc), intent(in)       :: s
      print *            , ""
      print "(a25)"      , "running stats"
      print "(25a)"      , repeat("-", 25)
      print "(a12, i13)"  , "n:"          , s%n
      print "(a12, f13.3)", "mean:"       , s%mean()
      print "(a12, f13.3)", "min:"        , s%min
      print "(a12, f13.3)", "max:"        , s%max
      print "(a12, f13.3)", "sum:"        , s%sum
      print "(a12, f13.3)", "stdev, p:"   , s%stdp()
      print "(a12, f13.3)", "stdev, s:"   , s%stds()
      print "(a12, f13.3)", "skewness, p:", s%skewp()
      print "(a12, f13.3)", "skewness, s:", s%skews()
      print "(a12, f13.3)", "kurtosis, p:", s%kurtp()
      print "(a12, f13.3)", "kurtosis, s:", s%kurts()
      print "(25a)"      , repeat("-", 25)
    end subroutine show

end module stats

