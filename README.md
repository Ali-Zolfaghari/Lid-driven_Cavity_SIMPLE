# Lid-driven_Cavity_SIMPLE
The lid-driven cavity is a well-known benchmark problem for viscous incompressible fluid flow. We are dealing with a square cavity consisting of three rigid walls with no-slip conditions and a lid moving with a tangential unit velocity. The lower left corner has a reference static pressure of 0. In computational fluid dynamics (CFD), the SIMPLE algorithm is a widely used numerical procedure to solve the Navier–Stokes equations. SIMPLE is an acronym for Semi-Implicit Method for Pressure Linked Equations.
- The algorithm is iterative. The basic steps in the solution update are as follows:  
- 1. Set the boundary conditions. 
- 2. Compute the gradients of velocity and pressure. 
- 3. Solve the discretized momentum equation to compute the intermediate velocity field. 
- 4. Compute the uncorrected mass fluxes at faces. 
- 5. Solve the pressure correction equation to produce cell values of the pressure correction. 
- 6. Update the pressure field and use the under-relaxation factor for pressure. 
- 7. Update the boundary pressure corrections. 
- 8. Correct the face mass fluxes. 
- 9. Correct the cell velocities by the gradient of the pressure corrections and the vector of central coefficients for the discretized linear system representing the velocity equation and Vol is the cell volume. 
- 10. Update density due to pressure changes.
