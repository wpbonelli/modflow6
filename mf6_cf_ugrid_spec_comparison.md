### Comparison

| Aspect                     | UGRID-1.0                                  | MODFLOW 6                                              |
|:---------------------------|:-------------------------------------------|:-------------------------------------------------------|
| Scope                      | Unstructured mesh topology                 | Application-specific                                   |
| Compliance                 | Extends CF                                 | Extends CF + UGRID                                     |
| Grid Types                 | 1D/2D/3D unstructured                      | Structured (DIS), extends UGRID w/ layers (DIS, DISV)  |
| Metadata                   | Extends CF: cf_role, etc                   | Extends CF/UGRID: grid, model, varname, aux, etc.      |
| Dimensions                 | Extends CF                                 | x, y, z, time (DIS), defers to UGRID (DIS, DISV)       |
| Topology                   | Nodes, edges, faces, volumes               | Nodes, lay/col/row (DIS), defers to UGRID (DIS, DISV)  |
