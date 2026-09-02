# Mesh formats and generator portability

Status: feasibility note. Records what was measured while evaluating Netgen as a
possible alternative to Gmsh, so the expensive half of that work can be scoped
honestly. The motivation is licensing (see below), not mesh quality. No
migration is proposed here; nothing in the codebase depends on Netgen.

## The portability property we rely on

A configuration binds a logical name to a set of raw mesh attribute ids:

```json
{"name": "LoopDomain", "dim": 2, "attribute_ids": [2]}
```

`PhysicsSolver::MarkerFromGroup()` resolves that name using the group's own
`AttributeIds` and nothing else. The mesh's `attribute_sets` /
`bdr_attribute_sets` are deliberately not consulted, and the solver no longer
writes into them.

This matters because named sets are **not** a universal mesh feature. Of MFEM's
readers:

| Reader | Populates `attribute_sets`? | Source of names |
|---|---|---|
| MFEM v1.3+ | yes | explicit `attribute_sets` block |
| Gmsh (`.msh`) | yes | `$PhysicalNames` |
| Cubit (`.e`) | yes | block/sideset names |
| **Netgen** (`areamesh2`, 3D neutral) | **no** | none; bare integer attributes |
| VTK | no | none |

Had we resolved names through the mesh, an identical config would work on Gmsh
and silently resolve to nothing on Netgen. Keeping the config authoritative
makes group resolution independent of the generator.

## Measured: the read path is generator-independent

`test/test_solvers.cpp` writes the current-loop geometry twice, from one shared
`CurrentLoopGrid`, so the two files differ only in serialization:

- `CreateCurrentLoopMesh()` -> `MFEM mesh v1.0`
- `CreateCurrentLoopNetgenMesh()` -> `areamesh2` (Netgen 2D)

The test `"Magnetostatic loop inductance is mesh-format independent (Netgen)"`
runs the *same* config over the Netgen file and asserts the analytic ring
inductance. Result:

| Mesh format | Inductance [H] |
|---|---|
| MFEM v1.0 | 6.023785e-07 |
| Netgen `areamesh2` | 6.023785e-07 |

Bit-identical, and the test additionally asserts that the Netgen reader produced
**zero** named sets. Reading Netgen meshes therefore requires no code changes at
all. The test is tagged `[netgen]`.

### `areamesh2` format

1-based indices; note that boundary elements precede elements, the reverse of
the MFEM format:

```
areamesh2
<n_bdr>     then per line: attr v1 v2
<n_elem>    then per line: attr nverts v1 ... vn
<n_vert>    then per line: x y
```

## Measured: geometry emission works, but id provenance is the trap

Netgen 6.x (`C:\netgen`, Python 3.14 bindings) was then installed and the loop
was rebuilt as an actual **geometry** definition -- `SplineGeometry`, two
materials, a `SetDomainMaxH` on the conductor -- and meshed by Netgen.
`tools/netgen_current_loop.py` does this and emits `areamesh2`.

Result on a coarse unstructured mesh (510 vertices, 948 triangles):

| Source | Inductance [H] | vs analytic |
|---|---|---|
| analytic (GMD ring) | 6.0272e-07 | - |
| structured, MFEM + Netgen format | 6.023785e-07 | -0.06% |
| **Netgen-meshed geometry** | **6.019738e-07** | **-0.12%** |

So the end-to-end path is feasible: geometry in, correct physics out, with no
solver changes.

### The id-provenance trap, concretely

The predicted risk was real and took three attempts to get right. In Netgen's
2D Python API:

1. **Domain ids follow `SetMaterial` order.** These behaved as expected
   (`air` -> 1, `loop` -> 2).
2. **`bc="name"` numbers boundaries by order of first appearance**, not by any
   declared value. Declaring `updown` first gave it id 1 where 3 was wanted.
3. **`bc=<int>` does not control `Elements1D().index`.** Passing `bc=3` does
   *not* make segments come back as index 3. `index` is the **spline
   declaration index**, so an 8-spline geometry yields ids 1..8 regardless.
4. **`Append()` returns a 0-based spline index while `Elements1D().index` is
   1-based.** A clean off-by-one between the two numbering systems.

The working approach is to ignore `bc=` for numbering entirely: capture
`Append()`'s return value at declaration time, `+1` it, and map spline index ->
attribute yourself. The geometry declares the *role*; the code derives the
*number*. `tools/netgen_current_loop.py` asserts the resulting contract after
meshing, so any future drift fails loudly instead of producing a plausible
solve.

Note also that Netgen emits **interior** segments (the air/conductor interface)
as boundary elements, so MFEM sees a bdr attribute that no boundary condition
references. It must be given an id so it cannot collide with an exterior one.

### Implication for configs

`attribute_ids` are hand-maintained integers, and none of the four hazards above
produce an error -- they produce a *different mesh with the same shape*. A wrong
boundary id silently moves a Dirichlet condition. If Netgen becomes real, the
generator should emit a sidecar mapping names to the ids it actually assigned,
and `entity_groups` should be generated from that rather than typed by hand.

Netgen has **no `areamesh2` writer** of its own (`Export()` lists Gmsh2,
Neutral, Abaqus and others, but not this one), so the writer in
`tools/netgen_current_loop.py` is ours. The obvious alternative -- exporting
Gmsh2 and letting MFEM rebuild names from `$PhysicalNames` -- was tested and
does not work; see the next section.

## Measured: Netgen's Gmsh2 export does NOT restore named groups

The obvious way to avoid a sidecar is to have Netgen `Export()` Gmsh2 and let
MFEM's Gmsh reader rebuild `attribute_sets` from `$PhysicalNames`. This was
tried directly. **It does not work**, for three independent reasons:

1. **No `$PhysicalNames` section is written at all.** `mesh.SetMaterial()` and
   `mesh.SetBCName()` both succeed, but the exporter discards the names. The
   file contains only `$MeshFormat`, `$Nodes` and `$Elements`. Element tags
   carry the raw integer ids (`tags=(1,1)` and `tags=(2,2)`), so the *numbers*
   survive and the *names* do not -- which is exactly the situation the sidecar
   exists to fix.
2. **Boundary elements are dropped.** The writer announces itself as
   "Write Gmsh v2.xx Surface Mesh" and emits 948 triangles and **zero** line
   elements, losing all 86 boundary segments. MFEM consequently sees
   `bdr_attributes` max of 1, and our own validator rejects the config:
   `entity_groups[2].attribute_ids: Attribute 2 is out of range [1, 1]`. Without
   boundary attributes there is nowhere to hang the far-field Dirichlet
   condition, so this path cannot express the problem at all.
3. **The version header is unreadable by MFEM.** Netgen writes
   `2.000000 0 8`; MFEM's `ReadGmshMesh` aborts with "Gmsh file version < 2.2".
   Reason 3 alone is a one-line patch, but it masks reasons 1 and 2 until fixed.

So Gmsh2-from-Netgen is strictly worse than `areamesh2`: it loses the boundary
information that `areamesh2` preserves, and gains no names. The custom
`areamesh2` writer in `tools/netgen_current_loop.py` is the right path, and the
sidecar remains the way to carry names for generators that have no place to put
them.

## Why this was investigated: licensing

The motivation is not technical. Gmsh is GPL, this project is MIT, and Gmsh is
therefore used only as a **separately installed executable** invoked by
`examples/generate_meshes.sh`. It is not a build dependency: `CMakeLists.txt`
fetches only MFEM, Eigen, nlohmann_json and Catch2. Generated `.mesh` files are
committed so the examples run without Gmsh present.

**This arrangement is deliberate and load-bearing.** Linking a GPL mesher's
library into this codebase would impose GPL on it. The practical consequence is
that in-process meshing -- meshing directly from a geometry definition, remeshing
during AMR, parametric sweeps without subprocess round-trips -- is closed off
while Gmsh is the mesher. That constraint, not mesh quality, is what makes
alternatives worth evaluating.

Netgen is **LGPL**, which permits use from an MIT project via *dynamic* linking
with the relink freedom preserved. Any adoption of `nglib` would have to be
dynamically linked and LGPL-compliant; static linking is not an option.

### Permissively licensed alternatives

Surveyed while looking for something that would avoid the LGPL obligations
entirely:

| Mesher | License | Fit |
|---|---|---|
| TetGen | AGPL | no |
| CGAL mesh generation | GPL / commercial | no |
| Triangle (Shewchuk) | non-commercial use only | no |
| MMG | LGPL | same category as Netgen |
| Salome SMESH | LGPL | same category as Netgen |
| fTetWild | MPL-2.0 | permissive, but 3D tetrahedra only |
| poly2tri | BSD | 2D constrained Delaunay only; no sizing or grading |

No permissively licensed general-purpose mesher covers the 2D axisymmetric and
3D cases this project needs. The realistic options are therefore GPL Gmsh at
arm's length (the status quo) or LGPL Netgen dynamically linked.

## Gmsh's id model is better than Netgen's

Worth recording explicitly, because it is the opposite of what the Netgen work
above might suggest. Gmsh's `.geo` DSL binds a name and an id together in the
geometry source:

```
Physical Surface("Loop", 2) = {loop_surfaces()};
```

and entities can be selected *geometrically* rather than by id, so the binding
survives renumbering after boolean operations:

```
loop_surfaces() = Surface In BoundingBox{ ... };
```

`examples/current_loop/loop.geo` uses exactly this to stay stable across
`BooleanDifference`. Netgen's 2D API offers no equivalent -- ids come from
declaration order, which is what produced the four hazards above.

So the sidecar concept is a **Netgen tax**, not a general requirement. It exists
to replace naming machinery that Gmsh already provides and Netgen lacks. Any
migration should count it as a cost of moving, not as a gap in the current setup.

## Not measured

- 3D. Only the 2D `SplineGeometry` / `areamesh2` path was exercised. The 3D
  neutral format and CSG geometry are untested here.
- Curved geometry. The loop is all straight splines; Netgen's curved elements
  and MFEM's `curved_areamesh2` variant were not tried.
- Whether a **newer** Netgen or its OCC/`nggui` export path writes richer Gmsh
  output. Only the bundled `Export(..., "Gmsh2 Format")` of the installed
  6.x build was tested. Gmsh itself can also read Netgen neutral files and
  re-export them with physical groups, which would be a third route.
- **Linking `nglib`.** The decisive question for the licensing motivation is
  whether Netgen can be driven in-process as a dynamically linked LGPL library,
  since that is the capability GPL Gmsh denies. Only the file-level and Python
  paths were exercised here; `nglib.lib` was not linked or called.
- Mesh quality on realistic geometry. The current-loop example is a toy: one
  rectangle in air. Nothing here says how Netgen compares to Gmsh on actual
  winding geometry, which is where meshing robustness would actually be decided.

## Related

- `docs/open_boundary.md` - far-field truncation in the same current-loop model.
