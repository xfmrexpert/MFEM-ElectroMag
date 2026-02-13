#!/usr/bin/env python3
"""Analyze electromagnetic simulation results and compare to analytical solutions"""

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import sys

def read_vtu_file(filename):
    """Read a VTU file and extract field data"""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()

    points = vtk_to_numpy(output.GetPoints().GetData())

    # Extract fields
    fields = {}
    point_data = output.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        array_name = point_data.GetArrayName(i)
        array = vtk_to_numpy(point_data.GetArray(i))
        fields[array_name] = array

    return points, fields

def analyze_capacitor():
    """Analyze parallel plate capacitor results"""
    print("\n" + "="*60)
    print("CAPACITOR ANALYSIS")
    print("="*60)

    try:
        points, fields = read_vtu_file("build/results_electrostatic/simple_capacitor_000000.vtu")
    except:
        print("ERROR: Could not read capacitor results")
        return

    V = fields.get('V', None)
    E = fields.get('E', None)

    if V is None:
        print("ERROR: No voltage field found")
        return

    print(f"\nVoltage (V):")
    print(f"  Min: {np.min(V):.2f} V")
    print(f"  Max: {np.max(V):.2f} V")
    print(f"  Mean: {np.mean(V):.2f} V")

    if E is not None:
        E_mag = np.linalg.norm(E, axis=1)
        print(f"\nElectric Field (E):")
        print(f"  Min magnitude: {np.min(E_mag):.2f} V/m")
        print(f"  Max magnitude: {np.max(E_mag):.2f} V/m")
        print(f"  Mean magnitude: {np.mean(E_mag):.2f} V/m")

        # Expected E ~ 1000V / 0.01m = 100 kV/m = 100,000 V/m
        # (for parallel plates with 1000V across 10mm gap)
        expected_E = 100000  # V/m
        print(f"\nExpected E-field (parallel plates): ~{expected_E/1000:.0f} kV/m")
        print(f"Simulated mean E-field: {np.mean(E_mag)/1000:.2f} kV/m")

        # Check if we're in the right ballpark (within an order of magnitude)
        if np.mean(E_mag) > expected_E/10 and np.mean(E_mag) < expected_E*10:
            print("✅ E-field magnitude is reasonable")
        else:
            print("⚠️  E-field magnitude may be incorrect")

def analyze_solenoid():
    """Analyze solenoid results"""
    print("\n" + "="*60)
    print("SOLENOID ANALYSIS")
    print("="*60)

    try:
        points, fields = read_vtu_file("build/results_magnetostatic/solenoid_000000.vtu")
    except:
        print("ERROR: Could not read solenoid results")
        return

    A = fields.get('A', None)
    B = fields.get('B', None)

    if A is not None:
        print(f"\nVector Potential (A):")
        print(f"  Min: {np.min(A):.6f} Wb/m")
        print(f"  Max: {np.max(A):.6f} Wb/m")

    if B is not None:
        B_mag = np.linalg.norm(B, axis=1)
        print(f"\nMagnetic Field (B):")
        print(f"  Min magnitude: {np.min(B_mag):.6f} T")
        print(f"  Max magnitude: {np.max(B_mag):.6f} T")
        print(f"  Mean magnitude: {np.mean(B_mag):.6f} T")

        # Expected: B = μ₀ * n * I
        # From config: J = 166667 A/m^2, coil thickness = 0.03m
        # n*I ≈ J * thickness = 166667 * 0.03 = 5000 A/m
        # B = 4π×10^-7 * 5000 = 0.00628 T = 6.28 mT
        expected_B = 0.00628  # Tesla
        print(f"\nExpected B-field (inside solenoid): ~{expected_B*1000:.2f} mT")
        print(f"Simulated max B-field: {np.max(B_mag)*1000:.2f} mT")

        # Check if we're in the right ballpark
        if np.max(B_mag) > expected_B*0.5 and np.max(B_mag) < expected_B*2:
            print("✅ B-field magnitude is reasonable")
        else:
            print("⚠️  B-field magnitude may be incorrect")

def analyze_eddy_current():
    """Analyze eddy current results"""
    print("\n" + "="*60)
    print("EDDY CURRENT ANALYSIS")
    print("="*60)

    try:
        points, fields = read_vtu_file("build/results_magnetoquasistatic/eddy_current_000000.vtu")
    except:
        print("ERROR: Could not read eddy current results")
        return

    A = fields.get('A', None)
    B = fields.get('B', None)
    J = fields.get('J', None)

    if A is not None:
        print(f"\nVector Potential (A):")
        print(f"  Min: {np.min(A):.6e} Wb/m")
        print(f"  Max: {np.max(A):.6e} Wb/m")

    if B is not None:
        B_mag = np.linalg.norm(B, axis=1)
        print(f"\nMagnetic Field (B):")
        print(f"  Min magnitude: {np.min(B_mag):.6e} T")
        print(f"  Max magnitude: {np.max(B_mag):.6e} T")
        print(f"  Mean magnitude: {np.mean(B_mag):.6e} T")

    if J is not None:
        J_mag = np.linalg.norm(J, axis=1)
        print(f"\nEddy Current Density (J):")
        print(f"  Min magnitude: {np.min(J_mag):.2f} A/m²")
        print(f"  Max magnitude: {np.max(J_mag):.2f} A/m²")
        print(f"  Mean magnitude: {np.mean(J_mag):.2f} A/m²")

if __name__ == "__main__":
    analyze_capacitor()
    analyze_solenoid()
    analyze_eddy_current()

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
