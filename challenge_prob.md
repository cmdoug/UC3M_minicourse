# FreeFEM UC3M minicourse
Challenge problem Day 1

```freefem
int n = 100;
mesh Th = square(n, n, [10*x, -1+2*y]);
```

### Define variational formulation
```freefem
real alpha = 0.1;
varf advdiff(q, v) = int2d(Th) (alpha*(dx(q)*dx(v) + dy(q)*dy(v)) + (1-y^2)*dx(q)*v) // bilinear form
				   + on(4, q = y^3 - 3*y + 2); // homogeneous Dirichlet BCs
```
### Define finite element space and functions
```freefem
func Pk = P1;
fespace Vh(Th, Pk);
```

### Build algebraic structures
```freefem
matrix Ah = advdiff(Vh, Vh);
real[int] bh = advdiff(0, Vh);
```

### Solve and plot
```freefem
Vh qh;
qh[] = Ah^-1*bh;
plot(qh);
```