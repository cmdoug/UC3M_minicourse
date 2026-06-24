# Bifurcation analysis of Navier–Stokes system

This example reproduces the vortex shedding phenomenon behind a cylinder in incompressible flow using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*}
\frac{\partial \vec{u}}{\partial t} + \vec{u}\cdot\nabla \vec{u} + \nabla p - \frac{1}{Re}\nabla^2\vec{u} &= 0\\\\
\nabla\cdot \vec{u} &= 0
\end{align*}
$$

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega$. Solutions $`\left[\vec{u},p\right]^T`$ are then sought, in the appropriate space, such that for all test functions $`\left[\check{\vec{u}},\check{p}\right]^T`$,

$$
\begin{align*}
\int_\Omega\left\lbrace\check{\vec{u}}\cdot\left[\frac{\partial\vec{u}}{\partial t} + \left(\vec{u}\cdot\nabla\vec{u}\right)\right] - \left(\nabla\cdot\check{\vec{u}}\right)p + \frac{1}{Re}\nabla\check{\vec{u}}:\nabla\vec{u} - \check{p}\left(\nabla\cdot\vec{u}\right)\right\rbrace\mathrm{d}\vec{x}=0
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_NS.idp](./eqns_NS.idp).


## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Define working directory and number of processors:
```sh
export workdir=~/your/path/to/Navier-Stokes_problem/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf ~/your/path/to/Navier-Stokes_problem/eqns_NS.idp eqns.idp
ln -sf ~/your/path/to/Navier-Stokes_problem/settings_NS.idp settings.idp
````

## Build initial meshes

#### Build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/Navier-Stokes_problem/cylinder.md -mo $workdir/cylinder
```

## Perform parallel computations using `ff-bifbox`
### Continue base state along the parameter $Da$ from trivial solution

```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cylinder.msh -fo test -Re 1 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi test.base -mo test2 -fo test2 -Re 1 -pv 1

ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi test2.base -fo Reswp -pv 1 -param Re

ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi test2.base -fo Reswp_adapt -pv 1 -param Re -mo Reswp_adapt -maxcount 20

ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi test2.base -fo Reswp_adapt -pv 1 -param Re -mo Reswp_adapt -maxcount 20 -thetamax 0.01

ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.7i -eps_nev 10 -fi Reswp_adapt_15.base -so spectrum -fo eigenmodes -pv 1 -sym 0
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.7i -eps_nev 10 -fi Reswp_adapt_15.base -so spectrum -fo eigenmodes -pv 1 -sym 1
```


