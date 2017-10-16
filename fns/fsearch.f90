! MIT License
!
! Copyright (c) 2017 Lars Andersen Bratholm
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

subroutine fmanhattan_distance(A, B, D)

    implicit none

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(:,:), intent(in) :: B
    double precision, dimension(:,:), intent(inout) :: D

    integer :: na, nb 
    integer :: i, j

    na = size(A, dim=2)
    nb = size(B, dim=2)

!$OMP PARALLEL DO
    do i = 1, nb
        do j = 1, na
            D(j,i) = sum(abs(a(:,j) - b(:,i)))
        enddo
    enddo
!$OMP END PARALLEL DO

end subroutine fmanhattan_distance
