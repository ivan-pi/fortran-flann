program kdtree_benchmark_flann

  use iso_c_binding
  use flann_c
  implicit none

  real(c_float), allocatable :: data(:,:) ! n by m
  real(c_float), allocatable :: queries(:,:) ! r by m
  
  type(FLANNParameters) :: params
  type(c_ptr) :: index_id = c_null_ptr
  real(c_float) :: speedup
  integer(c_int) :: ires

  real :: t0, t1, t2
  integer(c_int), parameter :: nn = 1
  integer, parameter :: m(3) = [3, 8, 16]
  integer, parameter :: n = 10000, r = 1000
  integer(c_int), allocatable :: idxs(:,:)
  real(c_float), allocatable :: dists(:,:)
  integer :: i, s
  integer, allocatable :: seed(:)

  call random_seed(size=s)
  allocate(seed(s))
  seed = 1234
  call random_seed(put=seed)

  ! params comes initialized by default
  params%algorithm = FLANN_INDEX_KDTREE_SINGLE
  params%trees = 1
  params%log_level = FLANN_LOG_INFO

  allocate(idxs(nn,r))
  ! idxs = 0
  allocate(dists(nn,r))
  ! dists = 0

  do i = 1, 3

    write(*,*) "Dimension = ", m(i)

    ! FLANN expects a pointer to a rows x cols matrix stored in 
    ! row-major order (one feature/data-point on each row).
    ! Since Fortran follows column major ordering, we must
    ! swap the dimensions.

    allocate(data(m(i),n))
    call random_number(data)

    allocate(queries(m(i),r))
    call random_number(queries)

    ! Populate tree
    call cpu_time(t0)
    ! call c_f_pointer(index_id,index_int)
    index_id = flann_build_index_float(data,n,m(i),speedup,params)
    call cpu_time(t1)

    ! Query random vectors
    ires = flann_find_nearest_neighbors_index_float(index_id,queries,r,idxs,dists,nn,params)
    call cpu_time(t2)

    write(*,*) "Tree build (s): ", t1 - t0
    write(*,*) "Tree query (s): ", t2 - t1

    ires = flann_free_index_float(index_id,params)
    deallocate(data)
    deallocate(queries)

  end do

end program
