# FreeFEM Newton iteration example

Solves Liouville-Bratu-Gelfand equation with homogeneous Dirichlet BCs on L-corner domain.

Our nonlinear equation is:

$$
\dot{q}-\nabla^2q-\lambda \exp(q) = 0
$$

with homogeneous Dirichlet BCs. 

### Define mesh:
```freefem
int[int] lab = [1,1,1,1];
int n = 20;
mesh Th = square(n, n, label = lab);
Th = trunc(Th, x < 0.5 | y < 0.5, label = 1); // subtractive construction
```

### Define finite element space and functions
```freefem
func Pk = P1;
fespace Vh(Th, Pk); // must at least span $H^1(\Omega)$
```

### Define variational formulation

The weak form of the residual is, find $q\in{}V$ such that:

$$
\forall v\in V,\qquad{}R(q;\lambda) = \int_{\Omega}\nabla v\cdot\nabla q - \lambda\exp(q)v\mathrm{d}\boldsymbol{x}
$$

The weak form of the Jacobian is, find $\delta{}q\in{}V$ such that:

$$
\forall v\in V,\qquad{}J(q;\lambda)\delta{}q = \int_{\Omega}\nabla v\cdot\nabla (\delta q) - \lambda\exp(q)\delta{}qv\mathrm{d}\boldsymbol{x}
$$

```freefem
Vh qh; // Define discrete solution vector in the space Vh
real lambda = 20;
varf Res(dq, v) = int2d(Th) (dx(qh)*dx(v) + dy(qh)*dy(v) - lambda*exp(qh)*v)
				+ on(1, dq = 0); // Restricts to $H_0^1(\Omega)$

varf Jac(dq, v) = int2d(Th) (dx(dq)*dx(v) + dy(dq)*dy(v) - lambda*exp(qh)*dq*v) // Defines a(q, v)
				+ on(1, dq = 0); // Restricts to $H_0^1(\Omega)$
```


### Solve nonlinear system using Newton iteration
```freefem
real tol = 1.0e-14;
real[int] R = Res(0, Vh);
real err = R.l2;
int iter = 0;
cout << " Iteration = " + iter + "; ||R|| = " + err << endl;
while (err > tol){
	matrix J = Jac(Vh, Vh);
	real[int] dq = J^-1*R; // here we solve for -dq
	qh[] -= dq;
	R = Res(0, Vh);
	err = R.l2;
	iter++;
	cout << " Iteration = " + iter + "; ||R|| = " + err + "; ||q|| = " + qh[].l2 + "; ||dq|| = " + dq.l2 << endl;
}
```

### Plot solution
```freefem
plot(Th, qh, cmm = "Global solution");
```
