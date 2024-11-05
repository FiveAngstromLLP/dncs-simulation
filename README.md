# DNCS Configuration and Setup Guide

## Prerequisites

- Python 3.8 or higher
- Git
- Rust toolchain (for building Rust components)
- just command runner (installation: https://github.com/casey/just)
  ```bash
  cargo install just
  ```

## Configuration File (dncs.toml)

The DNCS configuration file contains critical simulation parameters:

```toml
[simulation]
# Name of the molecule to simulate
moleculename = "YGGFM"

# Amino acid sequence
sequence = "AAAAA"

# Simulation interface to use (options: openmm)
interface = "openmm"

# Number of simulation samples to generate
n_samples = 10

# Temperature in Kelvin
temp = 300.0

# List of forcefield files
forcefield = ["amber99sb.xml", "amber14/tip3pfb.xml"]

# Computation device (options: CPU, CUDA, OpenCL)
device = "CPU"

# Solvent box size in angstroms
solvent = 100

# Number of simulation steps
steps = 5000

# Friction coefficient for Langevin dynamics
gamma = 1.0

# Integration timestep in picoseconds
dt = 0.002

# Number of molecular dynamics steps
md_steps = 5000

# Grid size for spatial decomposition
grid = 5
```

## Configuration Setup

1. Create a dncs.toml file in your working directory
2. Configure the simulation parameters based on your requirements
3. Save the file - DNCS will automatically load these settings at runtime

## Using the Justfile

The justfile provides convenient commands for managing DNCS:

### Installation
```bash
just install
```
This command:
- Installs required Python packages from requirements.txt
- Installs maturin build tool
- Builds and installs Rust components using maturin develop in release mode

### Running DNCS
```bash
just run
```
This command:
- Displays current working directory
- Executes the main DNCS Python script (src/main.py)

Note: The justfile uses DNCS_FOLDER environment variable to maintain correct paths regardless of execution location.
