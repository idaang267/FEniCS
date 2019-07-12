X = VectorFunctionSpace(mesh,'DG',1)

# Within Loop
while (steps < tot_steps):
    FluxVector = project(Flux(u,mu),X)
    
