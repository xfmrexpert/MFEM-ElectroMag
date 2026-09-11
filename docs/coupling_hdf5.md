# Coupling matrices in HDF5

`analysis_type: "coupling_matrix"` writes HDF5 instead of CSV. Console tables
remain available. Field-output flags do not control coupling output.

Files are written under `simulation.results_path`, or next to the mesh if that
setting is omitted. Missing output directories are created. Each save replaces
the file for that physics type; it does not append results from previous runs.

| Physics | File | Matrix groups under `/coupling` |
|---------|------|--------------------------------|
| Electrostatics | `coupling_electrostatics.h5` | `Capacitance` |
| Magnetostatics | `coupling_magnetostatics.h5` | `Inductance` |
| Magnetoquasistatics | `coupling_magnetoquasistatics.h5` | `Inductance`, `Resistance` |

## Schema version 1

The file root has an integer `schema_version` attribute equal to `1`.
The `/coupling` group has string attributes `physics_type` and `geometry_type`,
using the same values as the configuration.

`/coupling/terminal_names` is a one-dimensional array of variable-length UTF-8
strings. Its order is the solver's name-sorted terminal order and labels **both**
matrix axes. Rows are measured terminals; columns are driven terminals.
Terminal names are data, not HDF5 paths, so punctuation in names is preserved.

Each matrix group contains a `values` dataset of 64-bit floating-point numbers.
For static solvers its shape is `[N, N]`, where `N` is the terminal count:

- `/coupling/Capacitance/values`
- `/coupling/Inductance/values`

The `values` dataset has a string `units` attribute. Units are `F`, `H`, or `Ohm`
for axisymmetric models and `F/m`, `H/m`, or `Ohm/m` for planar models. Values
are stored without rounding to console precision.

## MQS frequency axis

`/coupling/frequency_hz` is a one-dimensional array of 64-bit floating-point
frequencies, with `units = "Hz"`. It contains unique, positive frequencies in
ascending numerical order. Duplicate frequencies across scenarios or within a
sweep are solved and stored once. Equality is exact numeric equality, not a
formatted frequency-string comparison. Scenario names are not dictionary keys.

Both `/coupling/Inductance/values` and `/coupling/Resistance/values` have shape
`[K, N, N]`, where `K` is the frequency count. The first dimension always exists,
even for a single-frequency solve. Slice `i` in each dataset belongs to
`frequency_hz[i]`. These are the extracted L and R matrices; an impedance matrix
can be constructed as `Z(f) = R(f) + j * 2 * pi * f * L(f)`.

## Reading into a C# dictionary

Use an HDF5 library to read `frequency_hz` as `double[]`, `terminal_names` as
`string[]`, and each `values` dataset as either a rectangular `double[,,]` or a
flat `double[]`, depending on the library. Check `schema_version` and the dataset
shapes before consuming the arrays.

For flat arrays, HDF5's last dimension varies fastest. The element for frequency
index `i`, row `r`, and column `c` is at `(i * N + r) * N + c`. Do not transpose
based on MFEM's internal column-major storage: the writer has already converted
it to the documented order.

After reading flat arrays named `frequencies`, `inductanceValues`, and
`resistanceValues`, and string array `terminalNames`, the dictionary assembly is
independent of the HDF5 library:

```csharp
var matrices = new Dictionary<double, (double[,] Inductance, double[,] Resistance)>();
int n = terminalNames.Length;
for (int i = 0; i < frequencies.Length; ++i)
{
	var inductance = new double[n, n];
	var resistance = new double[n, n];
	for (int r = 0; r < n; ++r)
	{
		for (int c = 0; c < n; ++c)
		{
			int offset = (i * n + r) * n + c;
			inductance[r, c] = inductanceValues[offset];
			resistance[r, c] = resistanceValues[offset];
		}
	}
	matrices.Add(frequencies[i], (inductance, resistance));
}
```

Use the stored frequency values as keys, rather than reconstructing them from
rounded display labels. Iterate the frequency array when sorted order matters;
do not rely on dictionary enumeration order.
