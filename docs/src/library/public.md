# Public Documentation

Documentation for `AerostructuralDynamics.jl`'s public interface.

See the Models section of the manual for documentation covering the individual models.

## Contents

```@contents
Pages = ["public.md"]
```

## Index

```@index
Pages = ["public.md"]
```

## Public Interface

### Getting Model Properties

The following methods may be used to determine the properties of a model.

```@docs
number_of_states
number_of_inputs
number_of_parameters
```

### Getting Model Indices

The following methods may be used to find the state, input, and parameter indices associated with each model in a coupled model.

```@docs
state_indices
input_indices
parameter_indices
```

### Getting Model Inputs

The following methods may be used to calculate the value of the inputs for a coupled model.

```@docs
get_inputs
get_inputs!
```

### Getting State Rates

The following methods may be used to calculate the value of the state rates for a model.

```@docs
get_rates
get_rates!
```

### Getting Mass Matrices

The following methods may be used to calculate the mass matrix for a model.

```@docs
get_mass_matrix
get_mass_matrix!
```

### Getting Jacobians

The following methods may be used to calculate the jacobian of the state rates for a model.

```@docs
get_state_jacobian
get_state_jacobian!
```

### Performing a Stability Analysis

The following method may be used to perform a stability analysis for a model.

```@docs
get_eigen
```

### Interfacing with DifferentialEquations

The following method may be used to construct a function for use with DifferentialEquations.

```@docs
get_ode
```
