# FreeFEM "Hello world!" example

Solves Poisson's equation with homogeneous Dirichlet BCs on L-corner domain.

### Define mesh:
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
```

### Define variational formulation
```freefem
func f = 1;
varf Poisson(q, v) = int2d(Th) (dx(q)*dx(v) + dy(q)*dy(v)) // Defines a(q, v)
			 	   + int2d(Th) (f*v) // Defines b(v)
				   + on(1, q = 0); // Restricts to $H_0^1(\Omega)$
```

### Define finite element space and functions
```freefem
func Pk = P1;
fespace Vh(Th, Pk); // must at least span $H^1(\Omega)$
```

### Build algebraic structures
```freefem
matrix Ah = Poisson(Vh, Vh); // stiffness matrix
real[int] bh = Poisson(0, Vh); // equal to Mh*fh
```

### Solve linear system
```freefem
Vh qh; // Define discrete solution vector in the space Vh
qh[] = Ah^-1*bh;
plot(Th, qh, cmm = "Global solution");
```
