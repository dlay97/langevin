module normal
  contains

    function r8_normal_ab ( a, b )

    !*****************************************************************************80
    !
    !! R8_NORMAL_AB returns a scaled pseudonormal R8.
    !
    !  Discussion:
    !
    !    The normal probability distribution function (PDF) is sampled,
    !    with mean A and standard deviation B.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = rk ) A, the mean of the PDF.
    !
    !    Input, real ( kind = rk ) B, the standard deviation of the PDF.
    !
    !    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
    !
      implicit none

      integer, parameter :: rk = kind ( 1.0D+00 )

      real ( kind = rk ) a
      real ( kind = rk ) b
      real ( kind = rk ) r1
      real ( kind = rk ) r2
      real ( kind = rk ) r8_normal_ab
      real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = rk ) x

      call random_number ( harvest = r1 )
      call random_number ( harvest = r2 )
      x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

      r8_normal_ab = a + b * x

      return
    end
    
end module normal