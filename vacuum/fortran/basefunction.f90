!> \file
!> \brief Polynomials evaluation


!> \brief Get the Zernike polynomials \f$\hat{R}^{m}_{l}\f$ with zeroth, first derivatives
!>
!> The original Zernike polynomial is defined by
!> The Zernike polynomials take the form
!> \f{eqnarray*}{
!>    Z^{-m}_{l}(s,\theta) &= R^{m}_{l}(s) \sin m\theta,\\
!>    Z^{m}_{l}(s,\theta) &= R^{m}_{l}(s) \cos m\theta,
!> \f}
!> where \f$R^{m}_{l}(s)\f$ is a \f$l\f$-th order polynomial given by
!> \f{eqnarray*}{
!>    R_l^m(s) = \sum^{\frac{l-m}{2}}_{k=0} \frac{(-1)^k (l-k)!}
!>    {k!\left[\frac{1}{2}(l+m) - k \right]! \left[\frac{1}{2}(l-m) - k \right]!} s^{l-2k},
!> \f}
!> and is only non-zero for \f$l \ge m\f$ and even \f$l-m\f$.
!>
!> In this subroutine, \f$R^{m}_{l}(s) \f$ is computed using the iterative relationship
!> \f{eqnarray*}{
!> R_l^m(s) = \frac{2(l-1)(2l(l-2)s^2 -m^2 -l(l-2))R_{l-2}^m(s) - l (l+m-2) (l-m-2)R_{l-4}^m(s) }{(l+m)(l-m)(l-2)}
!> \f}
!>
!> For \f$ m=0 \f$ and \f$ m=1 \f$, a basis recombination method is used by defining new radial basis functions as
!> \f{eqnarray*}{
!>    \hat{R}_{0}^{0} &= 1,
!>    \hat{R}^{0}_{l} &= \frac{1}{l+1}R^{0}_{l} - \frac{(-1)^{l/2}}{l+1},
!>    \\
!>    \hat{R}_{1}^{1} &= s,
!>    \hat{R}^{1}_{l} &= \frac{1}{l+1}R^{1}_{l} - \frac{(-1)^{(l-1)/2}}{2} s.
!> \f}
!> so that the basis scales as \f$s^{m+2}\f$ except for \f$\hat{R}_{0}^{0}\f$ and \f$\hat{R}_{1}^{1}\f$, which are excluded from the representation of \f${A}_{\theta,m,n}\f$.
!> For \f$m\ge2\f$, the radial basis functions are only rescaled as
!> \f[
!>   \hat{R}^{m}_{l} = \frac{1}{l+1}R^{m}_{l}.
!> \f]
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value, first derivative of Zernike polynomial
subroutine get_zernike(r, lrad, mpol, zernike)

  implicit none

  real(kind=8),intent(in) :: r
  integer, intent(in) :: lrad, mpol
  real(kind=8), intent(out) :: zernike(0:lrad,0:mpol,0:1)

  real(kind=8), parameter :: zero = 0.0
  real(kind=8), parameter :: one  = 1.0
  real(kind=8), parameter :: two  = 2.0

  real(kind=8) ::    rm, rm1  ! r to the power of m'th and m-1'th
  real(kind=8) ::    factor1, factor2, factor3, factor4
  integer :: m, n  ! Zernike R^m_n

  rm = one  ! r to the power of m'th
  rm1 = zero ! r to the power of m-1'th
  zernike(:,:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m,0:1) = (/ rm, real(m)*rm1 /)
    endif

    if (lrad >= m+2) then
      zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
      zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
      zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
    enddo

    rm1 = rm
    rm = rm * r

  enddo

  do n = 2, lrad, 2
    zernike(n,0,0) = zernike(n,0,0) - (-1)**(n/2)
  enddo

  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1,0) = zernike(n,1,0) - (-1)**((n-1)/2) * real((n+1)/2) * r
      zernike(n,1,1) = zernike(n,1,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m,:) = zernike(n,m,:) / real(n+1)
    end do
  end do
end subroutine get_zernike

!> \brief Get the Zernike polynomials  \f$\hat{R}^{m}_{l}\f$ with zeroth, first, second derivatives
!>
!> See get_zernike for more detail.
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value, first/second derivative of Zernike polynomial
subroutine get_zernike_d2(r, lrad, mpol, zernike)

  implicit none

  real(kind=8),intent(in) :: r
  integer, intent(in) :: lrad, mpol
  real(kind=8), intent(out) :: zernike(0:lrad,0:mpol,0:2)

  real(kind=8), parameter :: zero = 0.0
  real(kind=8), parameter :: one  = 1.0
  real(kind=8), parameter :: two  = 2.0

  real(kind=8) ::    rm, rm1, rm2  ! r to the power of m'th, m-1'th and m-2'th
  real(kind=8) ::    factor1, factor2, factor3, factor4
  integer :: m, n  ! Zernike R^m_n

  rm = one  ! r to the power of m'th
  rm1 = zero ! r to the power of m-1'th
  rm2 = zero ! r to the power of m-2'th
  zernike(:,:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m,0:2) = (/ rm, real(m)*rm1, real(m*(m-1))*rm2 /)
      !write(0, *) m, m, r, zernike(m,m,:)
    endif

    if (lrad >= m+2) then
      zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
      zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1
      zernike(m+2,m,2) = real((m+2)**2*(m+1))*rm - real((m+1)*m*(m-1))*rm2
      !write(0, *) m+2, m, r, zernike(m+2,m,:)
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
      zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
      zernike(n, m, 2) = factor1 * (two*factor2*(two*r*zernike(n-2,m,1) + zernike(n-2,m,0)) &
                        +(factor2*r**2 - factor3)*zernike(n-2,m,2) - factor4*zernike(n-4,m,2))
      !write(0, *) n, m, r, zernike(n,m,:)
    enddo

    rm2 = rm1
    rm1 = rm
    rm = rm * r

  enddo
  do n = 2, lrad, 2
    zernike(n,0,0) = zernike(n,0,0) - (-1)**(n/2)
  enddo
  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1,0) = zernike(n,1,0) - (-1)**((n-1)/2) * real((n+1)/2) * r
      zernike(n,1,1) = zernike(n,1,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m,:) = zernike(n,m,:) / real(n+1)
    end do
  end do
end subroutine get_zernike_d2

!> \brief Get the Zernike polynomials \f$\hat{R}^{m}_{l}/r^m\f$
!>
!> See get_zernike for more detail.
!>
!> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
!> @param[in] lrad radial resolution
!> @param[in] mpol poloidal resolution
!> @param[out] zernike the value
subroutine get_zernike_rm(r, lrad, mpol, zernike)

  implicit none

  real(kind=8),intent(in) :: r
  integer, intent(in) :: lrad, mpol
  real(kind=8), intent(out) :: zernike(0:lrad,0:mpol)

  real(kind=8), parameter :: zero = 0.0
  real(kind=8), parameter :: one  = 1.0
  real(kind=8), parameter :: two  = 2.0

  real(kind=8) ::    factor1, factor2, factor3, factor4
  integer :: m, n  ! Zernike R^m_n

  zernike(:,:) = zero
  do m = 0, mpol
    if (lrad >= m) then
      zernike(m,m) = one
    endif

    if (lrad >= m+2) then
      zernike(m+2,m) = real(m+2)*r**2 - real(m+1)
    endif

    do n = m+4, lrad, 2
      factor1 = real(n)/real(n**2 - m**2)
      factor2 = real(4 * (n-1))
      factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
      factor4 = real((n-2)**2-m**2) / real(n-2)

      zernike(n, m) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m) - factor4*zernike(n-4,m))
    enddo

  enddo
  do n = 2, lrad, 2
    zernike(n,0) = zernike(n,0) - (-1)**(n/2)
  enddo
  if (mpol >= 1) then
    do n = 3, lrad, 2
      zernike(n,1) = zernike(n,1) - (-1)**((n-1)/2) * real((n+1)/2)
    enddo
  end if

  do m = 0, mpol
    do n = m, lrad, 2
      zernike(n,m) = zernike(n,m) / real(n+1)
    end do
  end do
end subroutine get_zernike_rm
