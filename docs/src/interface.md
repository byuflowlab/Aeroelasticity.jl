# Interface

## User Interface

```@docs
AbstractModel
number_of_states
number_of_inputs
number_of_parameters
state_indices
input_indices
parameter_indices
get_rates
get_rates!
get_mass_matrix
get_mass_matrix!
get_state_jacobian
get_state_jacobian!
get_inputs
get_inputs!
get_eigen
get_ode
```

# Developer Interface

```@docs
AerostructuralDynamics.NoStateModel
AerostructuralDynamics.inplaceness
AerostructuralDynamics.mass_matrix_type
AerostructuralDynamics.state_jacobian_type
AerostructuralDynamics.input_jacobian_type
AerostructuralDynamics.get_input_jacobian
AerostructuralDynamics.get_input_mass_matrix
AerostructuralDynamics.get_input_mass_matrix!
AerostructuralDynamics.get_input_state_jacobian
AerostructuralDynamics.get_input_state_jacobian!
```
