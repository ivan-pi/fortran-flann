
# Fortran FLANN binding

Fortran bindings to the [FLANN](https://github.com/mariusmuja/flann
) library for performing fast approximate nearest neighbor searches in high dimensional spaces. 

> :warning:
> WARNING: Development of the Fortran FLANN bindings is ongoing. This means the API is not yet stable and might be subject to changes. Several functions remain untested.
> :warning:

## Minimal usage example

```fortran
use flann
implicit none

integer, parameter :: nn = 5, d = 128, ndata = 10000, ntest = 1000

real, allocatable, dimension(:,:) :: dataset, testset, dists
integer, allocatable :: indexes(:,:)
type(flann_index) :: idx

! Allocate data and result arrays
allocate(dataset(d,n),testset(d,ntest))
allocate(indexes(nn,ntest),dists(nn,ntest))

! Insert some random values
call random_number(dataset)
call random_number(testset)

! Create a FLANN index instance
idx = flann_index(dataset,kdtree_index_params(trees=8))
call idx%build_index()

! Perform K-nearest neighbor search
call idx%knn_search(testset,indexes,dists,nn,search_params(checks=128))

end
```


## Installing FLANN

On Linux (Ubuntu) you can install FLANN using the command

```
sudo apt install libflann-doc libflann-dev libflann1.9
```

Assuming default install paths, you can use the command
```
<PDF Viewer> /usr/share/doc/flann/manual.pdf
```
where the `<PDF Viewer>` is a program like Atril, Okular, Evince, or others.

Windows users can follow the instructions provided in documentation of the [original FLANN project](https://github.com/mariusmuja/flann).

## Using FLANN in Fortran

We recommend trying the brand new Fortran package manager - [`fpm`](). To integrate FLANN in your project add the following lines to the appropriate sections of your manifest file:

```toml

```

A basic Makefile is also provided in case you don't want to use `fpm`.



