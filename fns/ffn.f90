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

subroutine fffn(features, seed, nmax, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:) :: means
    double precision, allocatable, dimension(:) :: std
    integer, allocatable, dimension(:) :: ordering

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p, l
    double precision :: maxdist, ub, d

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(means(nsamples))
    allocate(std(nsamples))
    allocate(ordering(nsamples))

    ! Begin preprocessing
    means = 0.0d0
    std = 0.0d0

    means = sum(features, dim=1) / nsamples
    !$OMP PARALLEL DO
    do i = 1, nsamples
        std(i) = sqrt(sum(features(i,:) - means(i)))
    enddo
    !$OMP END PARALLEL DO

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = sum((features(p,:) - features(seed,:))**2)
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = 0.0d0
            !$OMP PARALLEL DO PRIVATE(m, ub) REDUCTION(max:d)
            do k = 1, i-1
                m = ordering(k)
                ub = nfeatures * ((means(n) - means(m))**2 + (std(n) + std(m))**2)

                if (ub <= maxdist) then
                    continue
                endif

                d = max(d,sum((features(n,:) - features(m,:))**2))

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



    ! Generate output
    final_ordering = ordering(:nmax)

    ! Deallocate temporary
    deallocate(means)
    deallocate(std)
    deallocate(ordering)

end subroutine fffn

subroutine fiffn(features, seed, npartitions, nmax, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: npartitions
    integer, intent(in) :: nmax
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:,:) :: means
    double precision, allocatable, dimension(:,:) :: std
    integer, allocatable, dimension(:) :: ordering

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist, ub, d

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
    allocate(std(nsamples, npartitions))
    allocate(ordering(nsamples))

    ! Begin preprocessing
    means = 0.0d0
    std = 0.0d0

    !$OMP PARALLEL DO PRIVATE(idx)
    do i = 1, nsamples
        do j = 1, npartitions-1
            idx = (j-1) * partition_size
            means(i,j) = sum(features(i, idx+1:idx+partition_size)) / partition_size
            std(i,j) = sqrt(sum((features(i, idx+1:idx+partition_size)-means(i,j))**2) / partition_size)
        enddo

        means(i,npartitions) = sum(features(i, (npartitions-1)*partition_size+1:)) / remainder
        std(i,npartitions) = sqrt(sum((features(i, (npartitions-1)*partition_size+1:))**2) / remainder)
    enddo
    !$OMP END PARALLEL DO

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = sum((features(p,:) - features(seed,:))**2)
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = 0.0d0
            !$OMP PARALLEL DO REDUCTION(max:d) PRIVATE(m, ub)
            do k = 1, i-1
                m = ordering(k)
                ub = partition_size * sum((means(n,:npartitions-1) - means(m,:npartitions-1))**2 + &
					& (std(n,:npartitions-1) + std(m,:npartitions-1))**2)
                ub = ub + remainder * ((means(n,npartitions) - means(m,npartitions))**2 + &
                    & (std(n,npartitions) + std(m, npartitions))**2)

                if (ub <= maxdist) then
                    continue
                endif

                d = max(d,sum((features(n,:) - features(m,:))**2))
            if (d > maxdist) then
                idx = j
                maxdist = d
            endif

            enddo
            !$OMP END PARALLEL DO

            ! Having this outside the loop is somewhat faster
        enddo

        ordering(i) = ordering(idx)
        ordering(idx) = p
    enddo



    ! Generate output
    final_ordering = ordering(:nmax)

    ! Deallocate temporary
    deallocate(means)
    deallocate(std)
    deallocate(ordering)

end subroutine fiffn

subroutine fifafn(features, seed, npartitions, nmax, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: npartitions
    integer, intent(in) :: nmax
    integer, dimension(nmax), intent(out) :: final_ordering

    double precision, allocatable, dimension(:,:) :: means
    double precision, allocatable, dimension(:,:) :: var
    integer, allocatable, dimension(:) :: ordering

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p
    double precision :: maxdist, d

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

    ! Begin preprocessing
    means = 0.0d0
    var = 0.0d0

    !$OMP PARALLEL DO PRIVATE(idx)
    do i = 1, nsamples
        do j = 1, npartitions-1
            idx = (j-1) * partition_size
            means(i,j) = sum(features(i, idx+1:idx+partition_size)) / partition_size
            var(i,j) = sum((features(i, idx+1:idx+partition_size)-means(i,j))**2) / partition_size
        enddo

        means(i,npartitions) = sum(features(i, (npartitions-1)*partition_size+1:)) / remainder
        var(i,npartitions) = sum((features(i, (npartitions-1)*partition_size+1:))**2) / remainder
    enddo
    !$OMP END PARALLEL DO

    ! Initialize the ordering
    ordering = (/ (i, i = 1, nsamples) /)

    ! The first item will be the seed
    ordering(seed) = ordering(1)
    ordering(1) = seed


    ! TODO: parallelise in some way, ie. threadlocking maxdist
    ! Do the actual algorithm
    do i = 2, nmax
        p = ordering(i)
        maxdist = sum((features(p,:) - features(seed,:))**2)
        idx = i
        do j = i, nsamples
            n = ordering(j)
            d = 0.0d0
            !$OMP PARALLEL DO REDUCTION(max:d) PRIVATE(m)
            do k = 1, i-1
                m = ordering(k)
                d = max(d,partition_size * sum((means(n,:npartitions-1) - means(m,:npartitions-1))**2 + &
					& var(n,:npartitions-1) + var(m,:npartitions-1)) + &
                    & remainder * ((means(n,npartitions) - means(m,npartitions))**2 + &
                    & var(n,npartitions) + var(m,npartitions)))

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

    ! Generate output
    final_ordering = ordering(:nmax)

    ! Deallocate temporary
    deallocate(means)
    deallocate(var)
    deallocate(ordering)

end subroutine fifafn

subroutine fobf(features, seed, nmax, final_ordering)

    implicit none

    double precision, dimension(:,:), intent(in) :: features
    integer, intent(in) :: seed
    integer, intent(in) :: nmax
    integer, dimension(nmax), intent(out) :: final_ordering

    integer, allocatable, dimension(:) :: ordering

    integer :: partition_size, remainder, nsamples, nfeatures
    integer :: i, j, k, idx, n, m, p, l
    double precision :: maxdist, ub, d

    nsamples = size(features,dim=1)
    nfeatures = size(features, dim=2)

    ! Allocate temporary
    allocate(ordering(nsamples))

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
            d = 0.0d0
            !$OMP PARALLEL DO REDUCTION(max:d) PRIVATE(m)
            do k = 1, i-1
                m = ordering(k)
                d = max(sum((features(n,:) - features(m,:))**2),d)

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

    ! Generate output
    final_ordering = ordering(:nmax)

    ! Deallocate temporary
    deallocate(ordering)

end subroutine fobf
