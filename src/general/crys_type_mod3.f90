! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module crys_type_mod3

  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none

  real(rk), parameter, private :: sq2_i = 1.0d0/dsqrt(2.0d0)
  real(rk), parameter, private :: sq6_i = 1.0d0/dsqrt(6.0d0)
  real(rk), parameter, private :: twosq6_i = 2.0d0*dsqrt(6.0d0)

  public

contains

  ! this guy should be in matrix_operations_mod, but there was circular dependencies
  ! when I did it (rq 03/2023)
  subroutine tensor3ddecompose(mat, dev, skw, sph)

    ! Decompose matrix into deviatoric, skew, or spherical parts.

    ! Note that the three components of the decomposition are orthogonal, but
    !   the underlying basis is not orthonormal due to scaling.

    ! Note: dev(2) = sqrt(3/2)*mat(3,3), when mat is deviatoric.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! mat: Array of 3x3 matrices
    ! dev: Array of 5-vec representing symmetric portion (mpsim convention)
    ! skw: Array of 3-vec representing skew portion (mpsim convention)
    ! sph: Array of scalars representing one third of the trace

    real(rk), intent(in) :: mat(:, :, :)
    real(rk), intent(out), optional :: dev(:, :)
    real(rk), intent(out), optional :: skw(:, :)
    real(rk), intent(out), optional :: sph(:)

    !---------------------------------------------------------------------------

    if (present(dev)) then
      dev(1, :) = (mat(1, 1, :) - mat(2, 2, :))*sq2_i
      dev(2, :) = (mat(3, 3, :) + mat(3, 3, :) - mat(1, 1, :) - mat(2, 2, :))* &
          & sq6_i

      dev(3, :) = sq2_i*(mat(1, 2, :) + mat(2, 1, :))
      dev(4, :) = sq2_i*(mat(1, 3, :) + mat(3, 1, :))
      dev(5, :) = sq2_i*(mat(2, 3, :) + mat(3, 2, :))
    end if

    if (present(skw)) then
      skw(1, :) = 0.5d0*(mat(2, 1, :) - mat(1, 2, :))
      skw(2, :) = -0.5d0*(mat(1, 3, :) - mat(3, 1, :))
      skw(3, :) = 0.5d0*(mat(3, 2, :) - mat(2, 3, :))
    end if

    if (present(sph)) then
      sph = (mat(1, 1, :) + mat(2, 2, :) + mat(3, 3, :))*(1.0d0/3.0d0)
    end if

  end subroutine tensor3ddecompose

end module crys_type_mod3
