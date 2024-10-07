program cluster_probability
    implicit none
    ! Declare variables
    integer, parameter :: max_size = 20 ! Adjust based on expected cluster sizes
    real(kind=8), dimension(:), allocatable :: cluster_time,hist
    real(kind=8) :: total_time, time,t
    integer :: cluster_size, i, n
    character(len=256) :: line

    ! Initialize variables
    total_time = 0.0
    allocate(cluster_time(1:max_size))
    allocate(hist(1:max_size))
    cluster_time = 0.0

    ! Open the input file
    open(unit=10, file='cluster.txt', status='old', action='read')
    open(unit=12, file='hist.txt', status='unknown')
    ! Read data from file and accumulate time for each cluster size
    do
        read(10, '(A)', iostat=n) line
        if (n /= 0) exit  ! Exit loop if end of file is reached

        ! Parse the line to extract cluster time and size
        read(line, *) t, t, time, cluster_size  ! Read time and cluster size
        cluster_time(cluster_size) = cluster_time(cluster_size) + time
        total_time = total_time + time
    end do

    ! Close the file
    close(10)

    ! Print the probabilities of each cluster size
    !print *, 'Cluster Size', 'Probability'
    do i = 1, max_size
        if (cluster_time(i) > 0.0) then
            !print *, i, cluster_time(i) / total_time
            hist(i)=cluster_time(i) / total_time
            write(12, '(I3, F10.4)') i,hist(i) 
        end if
    end do

    ! Deallocate memory
    deallocate(cluster_time)

end program cluster_probability

