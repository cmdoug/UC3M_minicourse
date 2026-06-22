# FreeFEM + PETSc "Hello world!" example

Solves Poisson's equation with homogeneous Dirichlet BCs on L-corner domain.

### Load PETSc and hpddm:
```freefem
load "PETSc"
include "macro_ddm.idp"
NewMacro def(u)u EndMacro
NewMacro init(u)u EndMacro
```

### Define distributed and global meshes:
```freefem
int[int] lab = [1,1,1,1];
int n = 20;
mesh Th = square(n, n, label = lab);
Th = trunc(Th, x < 0.5 | y < 0.5, label = 1); // subtractive construction
/*
mesh Th = square(n, n/2, [x, 0.5*y], label = lab);
{
	mesh Th2 = square(n/2, n/2, [0.5*x, 0.5+0.5*y], label = lab);
	Th = Th + Th2; // additive construction
}
Th = change(Th, rmInternalEdges = true);
*/
mesh Thg = Th; // Define the global mesh
DmeshCreate(Th); // Partition and distribute mesh
```

### Define variational formulation
```freefem
func f = 1; // f can be any function in $L^2(\Omega)$ 
varf Poisson(q, v) = int2d(Th) (dx(q)*dx(v) + dy(q)*dy(v)) // Defines a(q, v)
			 	   + int2d(Th) (f*v) // Defines b(f, v)
				   + on(1, q = 0); // Restricts to $H_0^1(\Omega)$
```

### Define discrete finite element space and functions
```freefem
func Pk = P1;
fespace Vh(Th, Pk); // must at least span $H^1(\Omega)$
```

### Build algebraic structures
```freefem
Mat Ah;
MatCreate(Th, Ah, Pk);
Ah = Poisson(Vh, Vh); // stiffness matrix
real[int] bh = Poisson(0, Vh); // equal to Mh*fh
```

### Solve linear system
```freefem
Vh qh; // Define discrete solution vector in the space Vh
qh[] = Ah^-1*bh;
plotMPI(Th, qh, Pk, def, real, cmm = "Global solution");
```
