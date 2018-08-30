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

module fastmath

    implicit none

contains

function fast_erf(x2, x, n, m, k) result(val)

    double precision, dimension(n), intent(in) :: x2
    double precision, dimension(n), intent(in) :: x
    integer, intent(in) :: m
    integer, intent(in) :: n
    double precision, intent(in) :: k

    double precision, dimension(n) :: val

    val = -0.6956521739130435d0 * x2 - 1.1283791670955126d0 * x
    val = 1.0d0 - fast_exp(val, n, m, k)

    ! second order
    !val = 1.0d0 - 0.7897872340425532 * fast_exp(-0.9069444444444444 * x2-0.749958629819626*x,n,m,k) &
    !    & -0.21021276595744676 * fast_exp(-0.9069444444444444 * x2-2.5501372990470221*x,n,m,k)

    ! exact
    !val = erf(x)

end function fast_erf


function fast_exp(x, n, m, k) result(val)
    double precision, dimension(n), intent(in) :: x
    integer, intent(in) :: n
    integer, intent(in) :: m
    double precision, intent(in) :: k

    integer :: i
    double precision, dimension(n) :: val

    val = (1.0d0+k*x)

    ! Makes sure the expression converges
    where (val < 0.0d0) val = 0.0d0

    do i=2, m
        val = val * val
    enddo

    ! exact
    !val = exp(x)


end function fast_exp

end module fastmath


!! Algorithm 2 but subtracting the mean of each feature
!! in the preprocessing
!subroutine fffn(features, seed, nmax, nmem, final_ordering)
!
!    implicit none
!
!    double precision, dimension(:,:), intent(in) :: features
!    integer, intent(in) :: seed
!    integer, intent(in) :: nmax
!    integer, intent(in) :: nmem
!    integer, dimension(nmax), intent(out) :: final_ordering
!
!    double precision, allocatable, dimension(:) :: means
!    double precision, allocatable, dimension(:) :: feature_means
!    double precision, allocatable, dimension(:) :: std
!    integer, allocatable, dimension(:) :: ordering
!
!    integer :: partition_size, remainder, nsamples, nfeatures
!    integer :: i, j, k, idx, n, m, p, l
!    double precision :: maxdist, ub, d, huge_double, lb, tmp
!
!    nsamples = size(features,dim=1)
!    nfeatures = size(features, dim=2)
!
!    ! Allocate temporary
!    allocate(means(nsamples))
!    allocate(feature_means(nfeatures))
!    allocate(std(nsamples))
!    allocate(ordering(nsamples))
!
!    huge_double = huge(1.0d0)
!
!    ! Begin preprocessing
!    means = 0.0d0
!    std = 0.0d0
!
!    feature_means = sum(features, dim=1) / nsamples
!
!    !$OMP PARALLEL DO
!    do i = 1, nsamples
!        means(i) = sum(features(i,:) - feature_means) / nfeatures
!        std(i) = sqrt(sum((features(i,:) - feature_means - means(i))**2) / nfeatures)
!    enddo
!    !$OMP END PARALLEL DO
!
!    ! Initialize the ordering
!    ordering = (/ (i, i = 1, nsamples) /)
!
!    ! The first item will be the seed
!    ordering(seed) = ordering(1)
!    ordering(1) = seed
!
!
!    ! Do the actual algorithm
!    do i = 2, nmax
!        p = ordering(i)
!        maxdist = 0.0d0
!        idx = i
!        do j = i, nsamples
!            n = ordering(j)
!            d = huge_double
!            !$OMP PARALLEL DO PRIVATE(m, lb, ub, tmp) REDUCTION(min:d)
!            do k = max(1, i - nmem), i-1
!                if (d <= maxdist) then
!                    cycle
!                endif
!
!                m = ordering(k)
!                tmp = (means(n) - means(m))**2 + std(n)**2 + std(m)**2
!                ub = nfeatures * (tmp + 2*std(n)*std(m))
!
!                if (ub <= maxdist) then
!                    d = 0.0d0
!                    cycle
!                endif
!
!                lb = nfeatures * (tmp - 2*std(n)*std(m))
!
!                if (lb >= d) then
!                    cycle
!                endif
!
!                d = min(d,sum((features(n,:) - features(m,:))**2))
!
!            enddo
!            !$OMP END PARALLEL DO
!
!            if (d > maxdist) then
!                idx = j
!                maxdist = d
!            endif
!        enddo
!
!        ordering(i) = ordering(idx)
!        ordering(idx) = p
!    enddo
!
!
!    ! Generate output in python format
!    final_ordering = ordering(:nmax) - 1
!
!    ! Deallocate temporary
!    deallocate(means)
!    deallocate(feature_means)
!    deallocate(std)
!    deallocate(ordering)
!
!end subroutine fffn

! Algorithm 4 but subtracting the mean of each feature
! in the preprocessing
subroutine fiffn(features, seed, npartitions, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: npartitions
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:,:) :: means
    double precision, allocatable, dimension(:) :: feature_means
    double precision, allocatable, dimension(:) :: tmp
    double precision, allocatable, dimension(:,:) :: std
    integer, allocatable, dimension(:) :: ordering
    integer, allocatable, dimension(:) :: weights

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist, ub, d, huge_double, lb

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Repartition the features
    partition_size = int(nfeatures/npartitions)
    remainder = mod(nfeatures, partition_size)

    if (remainder == 0) then
        remainder = partition_size
    else
        partition_size = partition_size + 1
        remainder = mod(nfeatures, partition_size)
    endif

    ! Allocate temporary
    allocate(means(nsamples, npartitions))
    allocate(feature_means(nfeatures))
    allocate(std(nsamples, npartitions))
    allocate(weights(npartitions))
    allocate(tmp(npartitions))
    allocate(ordering(nsamples))

    huge_double = huge(1.0d0)

    ! Begin preprocessing
    means = 0.0d0
    std = 0.0d0

    feature_means = sum(features, dim=1) / nsamples

    !$OMP PARALLEL DO PRIVATE(idx)
    do i = 1, nsamples
        do j = 1, npartitions-1
            idx = (j-1) * partition_size
            means(i,j) = sum(features(i, idx+1:idx+partition_size) - feature_means(idx+1:idx+partition_size)) / partition_size
            std(i,j) = sqrt(sum((features(i, idx+1:idx+partition_size) - &
                & feature_means(idx+1:idx+partition_size) - means(i,j))**2) / partition_size)
        enddo

        means(i,npartitions) = sum(features(i, (npartitions-1)*partition_size+1:) - &
            & feature_means((npartitions-1)*partition_size+1:)) / remainder
        std(i,npartitions) = sqrt(sum((features(i, (npartitions-1)*partition_size+1:) &
            & - feature_means((npartitions-1)*partition_size+1:) - means(i,npartitions))**2) / remainder)
    enddo
    !$OMP END PARALLEL DO

    weights = partition_size
    weights(npartitions) = remainder
    ! end preprocessing

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m, ub, lb, tmp)
            do k = max(1, i - nmem), i-1
                if (d <= maxdist) then
                    cycle
                endif

                m = ordering(k)
                tmp = (means(n,:) - means(m,:))**2 + std(n,:)**2 + std(m,:)**2
                ub = sum(weights * (tmp + 2*std(n,:)*std(m,:)))

                if (ub <= maxdist) then
                    d = 0.0d0
                    cycle
                endif

                lb = sum(weights * (tmp - 2*std(n,:)*std(m,:)))

                if (lb >= d) then
                    cycle
                endif

                d = min(d,sum((features(n,:) - features(m,:))**2))

            enddo
            !$OMP END PARALLEL DO

            if (d > maxdist) then
                idx = j
                maxdist = d
            endif

        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo



    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(means)
    deallocate(feature_means)
    deallocate(std)
    deallocate(tmp)
    deallocate(ordering)
    deallocate(weights)

end subroutine fiffn

! Brute force approach (Algorithm 5 but for furthest neighbours)
subroutine fobf(features, seed, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    integer, allocatable, dimension(:) :: ordering

    integer :: nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist,  d, huge_double

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(ordering(nsamples))
    huge_double = huge(1.0d0)

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m)
            do k = max(1, i - nmem), i-1
                if (d <= maxdist) then
                    cycle
                endif
                m = ordering(k)
                d = min(sum((features(n,:) - features(m,:))**2),d)
            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(ordering)

end subroutine fobf

! Variation that uses an approximation to the l2 distance
subroutine fifafn(features, seed, npartitions, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: npartitions
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:,:) :: means
    double precision, allocatable, dimension(:) :: feature_means
    double precision, allocatable, dimension(:,:) :: var
    integer, allocatable, dimension(:) :: ordering
    integer, allocatable, dimension(:) :: weights

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist, d, huge_double

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Repartition the features
    partition_size = int(nfeatures/npartitions)
    remainder = mod(nfeatures, partition_size)

    if (remainder == 0) then
        remainder = partition_size
    else
        partition_size = partition_size + 1
        remainder = mod(nfeatures, partition_size)
    endif

    ! Allocate temporary
    allocate(means(nsamples, npartitions))
    allocate(feature_means(nfeatures))
    allocate(var(nsamples, npartitions))
    allocate(weights(npartitions))
    allocate(ordering(nsamples))

    huge_double = huge(1.0d0)

    ! Begin preprocessing
    means = 0.0d0
    var = 0.0d0

    feature_means = sum(features, dim=1) / nsamples

    !$OMP PARALLEL DO PRIVATE(idx)
    do i = 1, nsamples
        do j = 1, npartitions-1
            idx = (j-1) * partition_size
            means(i,j) = sum(features(i, idx+1:idx+partition_size) - feature_means(idx+1:idx+partition_size)) / partition_size
            var(i,j) = (sum((features(i, idx+1:idx+partition_size) - &
                & feature_means(idx+1:idx+partition_size) - means(i,j))**2) / partition_size)
        enddo

        means(i,npartitions) = sum(features(i, (npartitions-1)*partition_size+1:) - &
            & feature_means((npartitions-1)*partition_size+1:)) / remainder
        var(i,npartitions) = (sum((features(i, (npartitions-1)*partition_size+1:) &
            & - feature_means((npartitions-1)*partition_size+1:) - means(i,npartitions))**2) / remainder)
    enddo
    !$OMP END PARALLEL DO

    weights = partition_size
    weights(npartitions) = remainder
    ! end preprocessing

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        ! Redundant extra calculation, but might be better to start with a realistic guess
        ! than 0 depending on the parallelisation
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m)
            do k = max(1, i - nmem), i-1
                if (d <= maxdist) then
                    cycle
                endif

                m = ordering(k)
                d = sum(weights * ((means(n,:) - means(m,:))**2 + var(n,:) + var(m,:)))

            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            ! and allows for easy parallelisation
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(means)
    deallocate(feature_means)
    deallocate(var)
    deallocate(ordering)
    deallocate(weights)

end subroutine fifafn

! l1-distance brute force approach, with distances stored in memory
subroutine fobf_l1_inmem(features, seed, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    integer, allocatable, dimension(:) :: ordering
    real, allocatable, dimension(:) :: distances

    integer :: nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p,r,s
    double precision :: maxdist,  d, huge_double

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(ordering(nsamples))
    allocate(distances((nsamples*(nsamples-1))/2))
    huge_double = huge(1.0d0)

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed

    !$OMP PARALLEL DO SCHEDULE(dynamic)
    do i = 1, nsamples-1
        do j = i+1, nsamples
            distances(((2*i-2)*nsamples+2*j-i**2-i)/2) = sum(abs(features(i,:)-features(j,:)))
        enddo
    enddo
    !$OMP END PARALLEL DO

    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m)
            do k = max(1, i - nmem), i-1
                m = ordering(k)
                if (n > m) then
                    d = min(d, distances(((2*m-2)*nsamples+2*n-m**2-m)/2))
                else
                    d = min(d, distances(((2*n-2)*nsamples+2*m-n**2-n)/2))
                endif

            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(ordering)

end subroutine fobf_l1_inmem

! l2-distance brute force approach, with distances stored in memory
subroutine fobf_l2_inmem(features, seed, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    integer, allocatable, dimension(:) :: ordering
    real, allocatable, dimension(:) :: distances

    integer :: nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p,r,s
    double precision :: maxdist,  d, huge_double

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(ordering(nsamples))
    allocate(distances((nsamples*(nsamples-1))/2))
    huge_double = huge(1.0d0)

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed

    !$OMP PARALLEL DO SCHEDULE(dynamic)
    do i = 1, nsamples-1
        do j = i+1, nsamples
            distances(((2*i-2)*nsamples+2*j-i**2-i)/2) = sum((features(i,:)-features(j,:))**2)
        enddo
    enddo
    !$OMP END PARALLEL DO

    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m)
            do k = max(1, i - nmem), i-1
                m = ordering(k)
                if (n > m) then
                    d = min(d, distances(((2*m-2)*nsamples+2*n-m**2-m)/2))
                else
                    d = min(d, distances(((2*n-2)*nsamples+2*m-n**2-n)/2))
                endif
            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(ordering)

end subroutine fobf_l2_inmem

! l1-distance brute force approach
subroutine fobf_l1(features, seed, nmax, nmem, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    integer, allocatable, dimension(:) :: ordering

    integer :: nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist,  d, huge_double

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(ordering(nsamples))
    huge_double = huge(1.0d0)

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m)
            do k = max(1, i - nmem), i-1
                if (d <= maxdist) then
                    cycle
                endif
                m = ordering(k)
                d = min(sum(abs(features(n,:) - features(m,:))),d)
            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(ordering)

end subroutine fobf_l1

! Variation that uses an approximation to the l1 distance
! Slightly more complicated than the l2 case but fast implementations
! of erf and exp should keep the speed reasonable
subroutine fifafn_l1(features, seed, npartitions, nmax, exp_approx_factor, nmem, final_ordering)

    use fastmath, only: fast_erf, fast_exp

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: npartitions
    integer, intent(in) :: nmax
    integer, intent(in) :: exp_approx_factor
    integer, intent(in) :: nmem
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:,:) :: means
    double precision, allocatable, dimension(:,:) :: var
    double precision, allocatable, dimension(:) :: sigma, mu, expon, sigma2, feature_means
    integer, allocatable, dimension(:) :: ordering, weights

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist, d, pi, sqrt_2_over_pi, one_over_sqrt_2, c1, huge_double, dd

    pi = 4.0d0 * atan(1.0d0)
    huge_double = huge(1.0d0)

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Repartition the features
    partition_size = int(nfeatures/npartitions)
    remainder = mod(nfeatures, partition_size)

    if (remainder == 0) then
        remainder = partition_size
    else
        partition_size = partition_size + 1
        remainder = mod(nfeatures, partition_size)
    endif

    ! Allocate temporary
    allocate(means(nsamples, npartitions))
    allocate(var(nsamples, npartitions))
    allocate(ordering(nsamples))
    allocate(sigma(npartitions))
    allocate(sigma2(npartitions))
    allocate(mu(npartitions))
    allocate(feature_means(nfeatures))
    allocate(expon(npartitions))
    allocate(weights(npartitions))

    ! Begin preprocessing
    means = 0.0d0
    var = 0.0d0
    weights = partition_size
    weights(npartitions) = remainder

    feature_means = sum(features, dim=1) / nsamples

    !$OMP PARALLEL DO PRIVATE(idx)
    do i = 1, nsamples
        do j = 1, npartitions-1
            idx = (j-1) * partition_size
            means(i,j) = sum(features(i, idx+1:idx+partition_size) - feature_means(idx+1:idx+partition_size)) / partition_size
            var(i,j) = (sum((features(i, idx+1:idx+partition_size) - &
                & feature_means(idx+1:idx+partition_size) - means(i,j))**2) / partition_size)
        enddo

        means(i,npartitions) = sum(features(i, (npartitions-1)*partition_size+1:) - &
            & feature_means((npartitions-1)*partition_size+1:)) / remainder
        var(i,npartitions) = (sum((features(i, (npartitions-1)*partition_size+1:) &
            & - feature_means((npartitions-1)*partition_size+1:) - means(i,npartitions))**2) / remainder)
    enddo
    !$OMP END PARALLEL DO

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    sqrt_2_over_pi = sqrt(2.0d0/pi)
    one_over_sqrt_2 = 1.0d0 / sqrt(2.0d0)

    c1 = 0.5d0**exp_approx_factor

    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = 0.0d0
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = huge_double
            !$OMP PARALLEL DO REDUCTION(min:d) PRIVATE(m, mu, sigma2, sigma, expon, dd)
            do k = max(1, i - nmem), i-1
                if (d <= maxdist) then
                    cycle
                endif
                m = ordering(k)
                mu = abs(means(n,:) - means(m,:))
                sigma2 = var(n,:) + var(m,:)
                sigma = sqrt(sigma2)
                expon = 0.5d0 * mu**2 / sigma2
                d = sum(weights * (sqrt_2_over_pi * sigma * fast_exp(-expon, npartitions, exp_approx_factor, c1) + &
                    & abs(mu) * fast_erf(expon, sqrt(expon), npartitions, exp_approx_factor, c1)))
                !d = sum(weights * (sqrt_2_over_pi * sigma * exp(-expon) + &
                !    & abs(mu) * erf(sqrt(expon))))
                !write(*,*) d, dd,sum(abs(features(n,:) - features(m,:)))

            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo

    ! Generate output in python format
    final_ordering = ordering(:nmax) -1

    ! Deallocate temporary
    deallocate(means)
    deallocate(feature_means)
    deallocate(var)
    deallocate(ordering)
    deallocate(sigma)
    deallocate(sigma2)
    deallocate(mu)
    deallocate(expon)
    deallocate(weights)

end subroutine fifafn_l1
