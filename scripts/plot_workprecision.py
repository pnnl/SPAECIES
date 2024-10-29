import numpy as np 
import matplotlib.pyplot as plt 
import netCDF4
import glob
import itertools

#################################### User Options ########################################
# LOAD_PATH = location of results to load in
# REFERENCE = reference solution filepath
# SAVE_PATH = save filepath
# SIM_NAMES = collections of simulations to load in. see SIMULATION_NAME in workprecision_runscript.sh
# MRI_FLAG  = indicates whether the collection of simulations is MRI or not...
# LABELS    = corresponding plot labels for each collection
# SOL_FIELD = solution field for which to compute work-precision diagram. options: nr, qr, T, q
# MAX_ERR   = threshold above which relative error is assumed infinite and not plotted 
# MIN_ERR   = threshold below which relative error is not plotted (assumed stagnated)
# USE_TEX   = flag for toggling LaTeX in plots (currently false since TeX not installed on tux box)
##########################################################################################
LOAD_PATH = "/home/dong9/SPAECIES-feature-arg-parser/build/rainshaft/results"
REFERENCE = LOAD_PATH + "/rainshaft_reference_notables.nc" 
SAVE_PATH = LOAD_PATH    
SIM_NAMES = ["explicit", "mri"] 
MRI_FLAG = [False, True]
LABELS = ["ERK", "MRI"] 
SOL_FIELD = "nr" 
MAX_ERR = 1.e1   
MIN_ERR = 1.e-13 
USE_TEX = False

# dictionaries for storing L2 error, time step size, and wall clock time
errors_sim = {}
timesteps_sim = {}
walltimes_sim = {}
for SIM_NAME in SIM_NAMES:
    errors_sim[SIM_NAME] = {}
    timesteps_sim[SIM_NAME] = {}
    walltimes_sim[SIM_NAME] = {}
    
## load reference solution
ref_data = netCDF4.Dataset(REFERENCE)
Z = ref_data["z_int"][-1,:]
dZ = [np.abs(z_right - z_left) for z_left, z_right in zip(Z, Z[1:])]
sol_ref = ref_data[SOL_FIELD][-1,:]
L2_ref = np.sqrt(np.sum(dZ * sol_ref**2))

## for each collection of simulations, compute L2 errors and extract 
## wall clock time and time step info
for idx, SIM_NAME in enumerate(SIM_NAMES):
    filenames = glob.glob(LOAD_PATH + '/rainshaft_' + SIM_NAME + '*.nc')
    
    # find orders for the collection and slow time steps if MRI
    orders = []
    if MRI_FLAG[idx] == True:
        slow_timesteps = []
        
    for filename in filenames:
        order = (filename.split("order", 1)[1])[0] # NOTE: This will break if order >=10...
        orders.append(order) if order not in orders else orders
        
        if MRI_FLAG[idx] == True:
            slow_timestep = (filename.split("dtslow", 1)[1]).split("_", 1)[0] 
            slow_timesteps.append(slow_timestep) if slow_timestep not in slow_timesteps else slow_timesteps
        
    # organize data according to order
    for order in orders:
        if MRI_FLAG[idx] == True:
            filenames_order = glob.glob(LOAD_PATH + '/rainshaft_' + SIM_NAME + '_order' + order + '*.nc')

            timesteps_order = {}
            L2_errors_order = {}
            walltimes_order = {}
            
            for dt_slow in slow_timesteps:
                filenames_dtslow = [filename for filename in filenames_order if dt_slow in filename]
                print(filenames_dtslow)
                
                timesteps_order_dtslow = []
                L2_errors_order_dtslow = []
                walltimes_order_dtslow = []
                
                for filename in filenames_dtslow:
                    # extract fast time step
                    dt = (filename.split("dt")[2])[0:-3]
                    timesteps_order_dtslow.append(float(dt)) if dt not in timesteps_order_dtslow else timesteps_order_dtslow
                    
                    # compute relative L2 error
                    data = netCDF4.Dataset(filename)
                    Z = data["z_int"][-1,:]
                    sol = data[SOL_FIELD][-1,:]
                    L2 = np.sqrt(np.sum(dZ * (sol - sol_ref)**2)) / L2_ref

                    # filter out large errors (assumed unstable) and very small errors
                    # (assumed stagnated)
                    if L2 > MAX_ERR:
                        L2_errors_order_dtslow.append(np.nan)
                    elif L2 < MIN_ERR:
                        L2_errors_order_dtslow.append(np.nan)
                    else:
                        L2_errors_order_dtslow.append(np.sqrt(np.sum(dZ * (sol - sol_ref)**2)) / L2_ref)
                    
                    # extract wall clock time (in seconds)
                    walltimes_order_dtslow.append(data['walltime_ms'][:]/1000.0)
                    
                # sort and save data for each slow time step
                sort_index = np.argsort(timesteps_order_dtslow)
                L2_errors_order[dt_slow] = np.array(L2_errors_order_dtslow)[sort_index]
                timesteps_order[dt_slow] = np.array(timesteps_order_dtslow)[sort_index]
                walltimes_order[dt_slow] = np.array(walltimes_order_dtslow)[sort_index]
            
            errors_sim[SIM_NAME][order] = L2_errors_order
            timesteps_sim[SIM_NAME][order] = timesteps_order
            walltimes_sim[SIM_NAME][order] = walltimes_order
        else:
            filenames_order = glob.glob(LOAD_PATH + '/rainshaft_' + SIM_NAME + '_order' + order + '*.nc')

            timesteps_order = []
            L2_errors_order = []
            walltimes_order = []
            for filename in filenames_order:        
                # extract time step
                dt = (filename.split("dt", 1)[1])[0:-3]
                timesteps_order.append(float(dt)) if dt not in timesteps_order else timesteps_order
            
                # compute relative L2 error
                data = netCDF4.Dataset(filename)
                Z = data["z_int"][-1,:]
                sol = data[SOL_FIELD][-1,:]
                L2 = np.sqrt(np.sum(dZ * (sol - sol_ref)**2)) / L2_ref
                
                # filter out large errors (assumed unstable) and very small errors
                # (assumed stagnated)
                if L2 > MAX_ERR:
                    L2_errors_order.append(np.nan)
                elif L2 < MIN_ERR:
                    L2_errors_order.append(np.nan)
                else:
                    L2_errors_order.append(np.sqrt(np.sum(dZ * (sol - sol_ref)**2)) / L2_ref)
                
                # extract wall clock time (in seconds)
                walltimes_order.append(data['walltime_ms'][:]/1000.0)
                
            # sort and save
            sort_index = np.argsort(timesteps_order)
            errors_sim[SIM_NAME][order] = np.array(L2_errors_order)[sort_index]
            timesteps_sim[SIM_NAME][order] = np.array(timesteps_order)[sort_index]
            walltimes_sim[SIM_NAME][order] = np.array(walltimes_order)[sort_index]


## plot work-precision diagram
plt.rcParams.update({'font.size': 23})
plt.rcParams['text.usetex'] = USE_TEX

# markers to cycle through
marker1 = itertools.cycle(('o', 'v', 's', 'd')) 
marker2 = itertools.cycle(('o', 'v', 's', 'd')) 

# colors to cycle through
color1 = itertools.cycle(('tab:blue', 'tab:orange', 'tab:green', 'tab:purple'))
color2 = itertools.cycle(('tab:blue', 'tab:orange', 'tab:green', 'tab:purple'))

for idx, SIM_NAME in enumerate(SIM_NAMES):
    errors_sim_i = errors_sim[SIM_NAME]
    timesteps_sim_i = timesteps_sim[SIM_NAME]
    walltimes_sim_i = walltimes_sim[SIM_NAME]
    
    fig, ax = plt.subplots(1, 2, figsize=(20,8))
    
    if MRI_FLAG[idx] == True:
        for order in errors_sim_i:
            color_i = next(color1)
            for dt_slow in errors_sim_i[order]:
                marker_i = next(marker1)
                
                dt_slow_float = float(dt_slow)
                dt_slow_print = str(dt_slow_float) if dt_slow_float>0 else str(np.abs(dt_slow_float)) + "$\Delta t_{fast}$"
                
                ax[0].loglog(timesteps_sim_i[order][dt_slow], errors_sim_i[order][dt_slow], linestyle='-', marker=marker_i, label=LABELS[idx] + order + ", $\Delta t_{slow}=$" + dt_slow_print, color=color_i)
                ax[1].loglog(walltimes_sim_i[order][dt_slow], errors_sim_i[order][dt_slow], linestyle='-', marker=marker_i, label=LABELS[idx] + order + ", $\Delta t_{slow}=$" + dt_slow_print, color=color_i)
    else:
        for order in errors_sim_i:
            ax[0].loglog(timesteps_sim_i[order], errors_sim_i[order], 's-', label=LABELS[idx] + order)
            ax[1].loglog(walltimes_sim_i[order], errors_sim_i[order], 's-', label=LABELS[idx] + order)
            
    ax[0].set_xlabel("Time step $\Delta t$ (s)")
    ax[0].set_ylabel("Relative $\ell^{2}$ error")
    ax[1].set_xlabel("Wall clock time (s)")
    ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.suptitle("Work-precision diagram")
    
    ax[0].grid()
    ax[1].grid()
    plt.tight_layout()
    
    plt.savefig(SAVE_PATH + "/workprecision_" + SIM_NAME + ".png")