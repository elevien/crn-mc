from pylab import *
from events import *


def path_coupled(N,J,level,w,t_max):
    Np = 1000.
    Nt = 1000.
    clock = zeros(Nt)

    mesh_fine = int(N/(pow(J,level)))     # fine grid
    mesh_coarse  = int(N/(pow(J,level+1))) # coarse grid
    refinement = J


    X0 = zeros((Nt,mesh_fine))
    X1 = zeros((Nt,mesh_coarse))

    Y0 = zeros((Nt,mesh_fine))
    Y1 = zeros((Nt,mesh_coarse))

    X0[0] = Np*ones(mesh_fine)*(pow(J,level))
    Y0[0] = Np*ones(mesh_fine)*(pow(J,level))
    X1[0] = Np*ones(mesh_coarse)*(pow(J,level+1))
    Y1[0] = Np*ones(mesh_coarse)*(pow(J,level+1))
    x0 = Np*ones(mesh_fine)*(pow(J,level))
    y0 = Np*ones(mesh_fine)*(pow(J,level))
    x1 = Np*ones(mesh_coarse)*(pow(J,level+1))
    y1 = Np*ones(mesh_coarse)*(pow(J,level+1))
    system_state= [x0,y0,x1,y1]
    # compute initial reaction rates and times
    # reaction kinetics
    events = []
    order =2
    # create events
    for i in range(mesh_coarse):
        voxels_fine = []
        for m in range(J):
            voxel_fine = i*J+m
            voxels_fine.append(voxel_fine)
            voxel_coarse = i

            # common
            stoichiometric_coeffs =\
            [\
            -identity(mesh_fine)[voxel_fine],\
            -identity(mesh_fine)[voxel_fine],\
            -identity(mesh_coarse)[voxel_coarse],\
            -identity(mesh_coarse)[voxel_coarse]\
            ]

            reaction = Reaction_SplitCommon(mesh_fine,mesh_coarse,refinement,
            voxel_fine,voxel_coarse,2,
            stoichiometric_coeffs,system_state)
            events.append(reaction)

            # fine
            stoichiometric_coeffs =\
            [\
            -identity(mesh_fine)[voxel_fine],\
            -identity(mesh_fine)[voxel_fine],\
            0,0]\


            reaction = Reaction_SplitFine(mesh_fine,mesh_coarse,refinement,
            voxel_fine,voxel_coarse,2,
            stoichiometric_coeffs,system_state)
            events.append(reaction)

        # coarse
        stoichiometric_coeffs =\
        [0,0,\
        -identity(mesh_coarse)[voxel_coarse],\
        -identity(mesh_coarse)[voxel_coarse]\
        ]

        reaction =  Reaction_SplitCoarse(mesh_fine,mesh_coarse,refinement,voxels_fine,voxel_coarse,order,stoichiometric_coeffs,system_state)
        events.append(reaction)

    # diffusion
    for i in range(mesh_coarse-1):
        # left diffusion common
        voxel_out_fine = (i+1)*J+1
        voxel_in_fine = (i+1)*J
        voxel_out_coarse = i+1
        voxel_in_coarse = i

        reaction = Diffusion_SplitCommon(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)

        reaction = Diffusion_SplitFine(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)

        reaction = Diffusion_SplitCoarse(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)

        # right diffusion common
        voxel_out_fine  = (i+1)*J
        voxel_in_fine = (i+1)*J+1
        voxel_out_coarse = i
        voxel_in_coarse = i+1

        reaction = Diffusion_SplitCommon(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)

        reaction = Diffusion_SplitFine(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)

        reaction = Diffusion_SplitCoarse(mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,0,system_state)
        events.append(reaction)




        for m in range(J-1):

            # left diffusion
            left_diffusion = Diffusion(mesh_fine,i*J+m,i*J+m+1,0,system_state)
            events.append(left_diffusion)

            # right diffusion
            right_diffusion = Diffusion(mesh_fine,i*J+m+1,i*J+m,0,system_state)
            events.append(right_diffusion)

    clock[0] = 0.
    k = 1

    while k < Nt:

        x0[:] = X0[k-1]
        y0[:] = Y0[k-1]
        x1[:] = X1[k-1]
        y1[:] = Y1[k-1]
        system_state=[x0,y0,x1,y1]

        # find event with minimum absolute time
        firing_event = min(events, key=lambda e: e.wait_absolute)
        m = events.index(firing_event)
        # get information needed from event
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        # update events
        firing_event.fire(system_state,delta)
        # update all other events
        events.pop(m)
        for e in events:
            e.no_fire(system_state,delta)
        events.append(firing_event)
        # update system
        clock[k] = clock[k-1]+delta
        X0[k][:] = x0 + stoichiometric_coeffs[0]
        Y0[k][:] = y0 + stoichiometric_coeffs[1]
        X1[k][:] = x1 + stoichiometric_coeffs[2]
        Y1[k][:] = y1 + stoichiometric_coeffs[3]


        k = k+1
    return X0,Y0,X1,Y1,clock




def path_uncoupled(N,w,t_max):
    Np = 100.
    Nt = 3000.
    clock = zeros(Nt)

    X0 = zeros((Nt,N))
    Y0 = zeros((Nt,N))

    X0[0] = Np*ones(N)
    Y0[0] = Np*ones(N)
    x0 = Np*ones(N)
    y0 = Np*ones(N)
    system_state = [x0,y0]
    events = []

    # reactions
    for i in range(N):
        stoichiometric_coeffs =[-identity(N)[i],-identity(N)[i]]
        reaction = Reaction(N,i,2,stoichiometric_coeffs,system_state)
        events.append(reaction)

    # diffusions
    for i in range(N-1):

        # left diffusion
        left_diffusion = Diffusion(N,i,i+1,0,system_state)
        events.append(left_diffusion)

        # right diffusion
        right_diffusion = Diffusion(N,i+1,i,0,system_state)
        events.append(right_diffusion)

    clock[0] = 0.
    k = 1

    while k < Nt:

        x0[:] = X0[k-1]
        y0[:] = Y0[k-1]
        system_state=[x0,y0]

        # find event with minimum absolute time
        firing_event = min(events, key=lambda e: e.wait_absolute)
        m = events.index(firing_event)
        # get information needed from event
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        # update events
        firing_event.fire(system_state,delta)
        # update all other events
        events.pop(m)
        for e in events:
            e.no_fire(system_state,delta)
        events.append(firing_event)
        # update system
        clock[k] = clock[k-1]+delta
        X0[k][:] = x0 + stoichiometric_coeffs[0]
        Y0[k][:] = y0 + stoichiometric_coeffs[1]

        k = k+1
    return X0,Y0,clock
