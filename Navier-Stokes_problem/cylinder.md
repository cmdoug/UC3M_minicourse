# cylinder.md

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"

int n = 20;

// mesh filename, if not provided, defaults to "cylinder.msh"
string meshout = getARGV("-mo","cylinder.msh"); 
if(meshout.rfind(".msh") <0) meshout = meshout + ".msh"; // add extension if not provided


//	o------------------------------------o
//  |									 |
//  |									 |
//  |									 |
//  |			.--.					 |
//  o----1-----o    o--------------------o


border UpAxis(t = -20, -0.5) { x = t; y = 0; label = BCaxis; }
border Cyl(t = pi, 0) { x = 0.5*cos(t); y = 0.5*sin(t); label = BCcyl; }
border DownAxis(t = 0.5, 50) { x = t; y = 0; label = BCaxis; }
border DownStream(t = 0, 20) { x = 50.0; y = t; label = BCout; }
border Lateral(t = 50, -20) { x = t; y = 20.0; label = BClat; }
border UpStream(t = 20, 0) { x = -20.0; y = t; label = BCin; }

mesh Thg = buildmesh(UpAxis(n) + Cyl(2*n) + DownAxis(n) + DownStream(n) + Lateral(n) + UpStream(n));

int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
plot(Thg);
```