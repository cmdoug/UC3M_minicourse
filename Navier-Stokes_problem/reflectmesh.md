# Reflect mesh across centerline

```freefem
assert(mpirank == 0);
include "settings.idp"
include "macros_bifbox.idp"
string meshin = getARGV("-mi", ""); // input meshfile with extension
string meshout = getARGV("-mo", ""); // output meshfile without extension
string meshroot, meshext = parsefilename(meshin, meshroot);
parsefilename(meshout, meshout); // trim extension from output mesh, if given
meshout = meshout + "." + meshext;
Th = readmeshN(workdir + meshin); // load half mesh
plot(Th);
mesh Th2 = movemesh(Th, [x, -y]); // create reflected copy of mesh
plot(Th2);
Th = Th + Th2; // unify both halves of the mesh
plot(Th);
Th = change(Th, rmInternalEdges = true); // remove the internal edge along the centerline
    cout << "  Saving reflected mesh '" + meshout + "' in '" + workdir + "'." << endl;
savemesh(Th, workdir + meshout); // save full mesh
plot(Th);
```