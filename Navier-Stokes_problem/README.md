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
ln -sf $workdir/../eqns_NS.idp eqns.idp
ln -sf $workdir/../settings_NS.idp settings.idp
````

## Build initial meshes

#### Build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 $workdir/../cylinder.md -mo $workdir/cylinder
```

## Perform parallel computations using `ff-bifbox`
### Compute an initial base state at $Re=1$ from null initial guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cylinder.msh -fo test -Re 1 -pv 1
```

### Adapt mesh to initial solution at $Re=1$ and recompute state, exporting solution to ParaView file
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi test.base -mo test2 -fo test2 -Re 1 -pv 1
```

### Continue initial base state along $Re$ parameter (without adaptation)
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi test2.base -fo Reswp -pv 1 -param Re -maxcount 20
```

### Continue initial base state along $Re$ parameter (with adaptation at every step)
Note: the `-thetamax` parameter helps to prevent "coarsening" of the corners along the round cylinder surface during adaptation. 
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi test2.base -fo Reswp_adapt -pv 1 -param Re -mo Reswp_adapt -maxcount 20 -scount 1 -thetamax 0.01
```
### Compute leading 10 eigenvalues about a base state assuming both symmetric and anti-symmetric perturbations
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.7i -eps_nev 10 -fi Reswp_adapt_15.base -so spectrum -fo sym_modes -pv 1 -sym 0
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.7i -eps_nev 10 -fi Reswp_adapt_15.base -so spectrum -fo asym_modes -pv 1 -sym 1
```

### Scan across a line slightly above neutral stability threshold
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.1i -ntarget 10 -targetf 0.1+1i -eps_nev 20 -fi Reswp_adapt_15.base -so spectrum_swp -pv 1 -sym 0
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -eps_target 0.1+0.1i -ntarget 10 -targetf 0.1+1i -eps_nev 20 -fi Reswp_adapt_15.base -so spectrum_swp -pv 1 -sym 1
```

### Perform DNS
First, we must reflect the mesh
```sh
FreeFem++-mpi $workdir/../reflectmesh.md -dir $workdir -mi test2.msh -mo reflectedmesh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi reflectedmesh.msh -fo reflectedbase -Re 50 -pv 1
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -mi reflectedmesh.msh -fi asym_modes_5.mode -bfi reflectedbase.base -fo dns_solution -pv 1 -scount 5 -maxcount 1000 -amp 10
```

### Resolvent analysis
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -mi test2.msh -fi test2.base -sym 1 -eps_nev 5
```