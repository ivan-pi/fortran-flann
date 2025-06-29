program test_flann

    use flann_c
    use iso_c_binding

    implicit none

    type(FLANNParameters) :: p

    integer, parameter :: rows = 9000
    integer, parameter :: cols = 128
    integer, parameter :: tcount = 1000

    integer :: nn

    real, allocatable :: dataset(:,:), testset(:,:)

    integer, allocatable :: result(:,:)
    real, allocatable :: dists(:,:)

    real :: speedup, t1, t2, t3
    type(c_ptr) :: index_id
    integer, pointer :: index_int => null()
    integer :: ires

    allocate(dataset(cols,rows),testset(cols,tcount))

    write(*,*) "Reading input data file."
    call read_points("dataset.dat",rows,cols,dataset)

    write(*,*) "Reading test data file."
    call read_points("testset.dat",tcount,cols,testset)

    nn = 3

    p%algorithm = FLANN_INDEX_KDTREE_SINGLE
    p%trees = 1
    p%log_level = FLANN_LOG_INFO
    p%checks = 128
    p%random_seed = 12345678

    write(*,*) "Computing index."

    call cpu_time(t1)
    index_id = flann_build_index_float(dataset,rows,cols,speedup,p)
    call cpu_time(t2)
    print *, "Building tree ", t2 - t1, " s"

    call c_f_pointer(index_id,index_int)
    print*, "pointer to integer = ", index_int

    allocate(result(nn,tcount))
    allocate(dists(nn,tcount))

    call cpu_time(t2)
    ires = flann_find_nearest_neighbors_index(index_id,testset,tcount,result,dists,nn,p)
    call cpu_time(t3)
    print *, "Query time ", t3 - t2, " s"
    
    print *, "Find neighbors ", ires

    call print_points("results_f.dat",result,dists)


    ! deallocate(index_int)


    ires = flann_free_index(index_id,p)
    deallocate(dataset,testset,result,dists)

    print*, "Clean memory ", ires
    ! call print_params(p)

contains

    subroutine read_points(fname,rows,cols,points)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: rows, cols
        real, intent(out) :: points(cols,rows)
        integer :: i,j,unit

        open(newunit=unit,file=fname,status='old')

        do i = 1, rows
            read(unit,*) (points(j,i),j=1,cols)
        end do

        close(unit)
    end subroutine

    subroutine print_points(fname,points, dists)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: points(:,:)
        real, intent(in), optional :: dists(:,:)

        integer :: i, j, unit

        open(newunit=unit,file=fname,action="write")
        print *, "Writing data to file "//fname

        if (present(dists)) then
            do i = 1, size(points,2)
                write(unit,*) (points(j,i),j=1,size(points,1)), (dists(j,i),j=1,size(dists,1))
            end do
        else
            do i = 1, size(points,2)
                write(unit,'(*(i0,:," "))') (points(j,i),j=1,size(points,1))
            end do
        end if

        close(unit)
    end subroutine

    subroutine print_params(p)
        type(FLANNParameters), intent(in) :: p

        print '("algorithm ",(i0))',p%algorithm
        print '("checks ",(i0))',p%checks
        print '("eps ",(g0))',p%eps
        print '("sorted ",(i0))',p%sorted
        print '("max_neighbors ",(i0))',p%max_neighbors
        print '("cores ",(i0))',p%cores
        print '("trees ",(i0))',p%trees
        print '("leaf_max_size ",(i0))',p%leaf_max_size
        print '("branching ",(i0))',p%branching
        print '("iterations ",(i0))',p%iterations
        print '("centers_init ",(i0))',p%centers_init
        print '("cb_index ",(g0))',p%cb_index
        print '("target_precision ",(g0))',p%target_precision
        print '("build_weight ",(g0))',p%build_weight
        print '("memory_weight ",(g0))',p%memory_weight
        print '("sample_fraction ",(g0))',p%sample_fraction
        print '("table_number_ ",(i0))',p%table_number_
        print '("key_size_ ",(i0))',p%key_size_
        print '("multi_probe_level_ ",(i0))',p%multi_probe_level_
        print '("log_level ",(i0))',p%log_level
        print '("random_seed ",(i0))',p%random_seed        

    end subroutine

end program