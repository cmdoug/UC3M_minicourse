# This is a script that can be used to build FreeFEM locally on Macs

# First update homebrew and install necessary packages/dependencies
```sh
brew update
brew upgrade
brew install gfortran
brew install m4 bison hdf5 autoconf automake
```
# Set target directory for FreeFEM and PETSc (from Chris' forks!) and clone
```sh
cd $HOME/local
export LOCAL_DIR=${PWD}
git clone https://github.com/cmdoug/FreeFem-sources
```

# Build FreeFEM
```sh
cd ${LOCAL_DIR}/FreeFem-sources
git checkout develop
autoreconf -i
./configure --enable-download --enable-optim --prefix=${LOCAL_DIR}/FreeFem-sources
./3rdparty/getall -a
cd 3rdparty/getall && make petsc-slepc && cd -
./reconfigure
make -j4
make install
# if necessary, add `sudo` before the line above
```