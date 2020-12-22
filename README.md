
> :warning:
> WARNING: Development of the Fortran FLANN binding is ongoing. This means the API is not yet stable and might be subject to changes. Several functions remain untested. The instructions given below might not be complete.
> :warning:

# Fortran FLANN binding

Fortran bindings to the [FLANN](https://github.com/mariusmuja/flann
) library for performing fast approximate nearest neighbor searches in high dimensional spaces. 

* [Minimal usage example](#minimal-usage-example)
* [Installing FLANN](#installing-flann)
* [Using FLANN with fpm](#using-flann-with-fpm)
* [Learn more](#learn-more)
* [Contributing](#contributing)

## Minimal usage example

The example below shows how to use a `flann_index` instance to find the 5 nearest neighbors for 1000 random test points from a 128-dimensional data set containing 10000 points:

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
idx = flann_index(dataset,kdtree_index_params(trees=8),'Euclidean')
call idx%build_index()

! Perform K-nearest neighbor search
call idx%knn_search(testset,indexes,dists,nn,search_params(checks=128))

end
```

## Installing FLANN

On Linux you can install FLANN using the command

```
sudo apt install libflann-doc libflann-dev libflann1.9
```

Windows users can follow the instructions provided in documentation of the [original FLANN project](https://github.com/mariusmuja/flann).

## Using FLANN with `fpm`

To use FLANN in your project we recommed trying the new Fortran package manager - [`fpm`](https://github.com/fortran-lang/fpm). To integrate FLANN in your project add the following lines to the `[dependencies]` section of your TOML manifest file:

```toml
fortran-flann = { git = "https://github.com/ivan-pi/fortran-flann" }
```

## Learn more

On Linux systems, assuming you installed the target `libflann-doc` under default path settings, you can view the original project documentation with the command:
```
<PDF Viewer> /usr/share/doc/flann/manual.pdf
```
where the `<PDF Viewer>` is a program like Atril, Okular, Evince, or others. Since this package only contains Fortran bindings, for the most part the API follows the original project in C++, with the important difference we decided to use "snake case" across all methods.

A more complete explanation of the algorithms available in FLANN can be found in the paper: 

> Muja, M., & Lowe, D. G. (2009). Fast approximate nearest neighbors with automatic algorithm configuration. *Proceedings of the Fourth International Conference on Computer Vision Theory and Applications - Volume 1: VISAPP, (VISIGRAPP 2009)*, pages 331-340. DOI:[10.5220/0001787803310340](https://doi.org/10.5220/0001787803310340)

PDF versions can be found easily with your favorite search engine. Some working links at the time of writing (2020/12/22) include: [(link 1)](https://lear.inrialpes.fr/~douze/enseignement/2014-2015/presentation_papers/muja_flann.pdf), [(link 2)](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.160.1721&rep=rep1&type=pdf), [(link 3)](http://image.ntua.gr/iva/files/MujaLowe_ICCVTA2009%20-%20Fast%20Approximate%20Nearest%20Neighbors%20with%20Automatic%20Algorithm%20Configuration.pdf)

## Contributing

Feel welcome to submit bug reports or suggest changes to the Fortran bindings by opening a new [issue](https://github.com/ivan-pi/fortran-flann/issues).

If you think you are facing an issue with the underlying FLANN library, you might be able to find an answer in the [list of open/closed issues](https://github.com/mariusmuja/flann/issues) of the original FLANN project. Unfortunately, the original project has gone stale and doesn't seem to be supported anymore.

Since the Fortran bindings provided here are only a wrapper of the C interface exported in the original FLANN project, we are limited to a subset of the original FLANN functionality.




