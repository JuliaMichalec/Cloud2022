import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time

# methane- declaration of mixtures by gri30 mechanism
gas = ct.Solution('gri30.xml', 'gri30_mix')
T0 = np.linspace(300, 400, 5)
P0 = [i*ct.one_atm for i in range(1, 3)]
species_dict = {'CH4': 0.9, 'O2': 2, 'N2': 7.52}  # giving the mole fractions of the elements in mixture

# calculating the temperature, the flame speed and the CO_2 variation along the length of the domain
for Pi in P0:
    x_plt = []
    T_plt = []
    v_plt = []
    conc_CO2_plt = []
    X_CO2_plt = []
    for Ti in T0:
        gas.TPX = Ti, Pi, species_dict
        f = ct.FreeFlame(gas, width=0.05)
        f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
        f.solve(loglevel=1, refine_grid=True, auto=False)
        x = f.grid
        T = f.T
        v = f.velocity
        x_plt.append(x)
        T_plt.append(T)
        v_plt.append(v)
        print('===========================================================================================')
        print('Mixture-averaged flamespeed = %.6f m/s' % (v[0]))
        print('===========================================================================================')
        conc_CO2 = f.concentrations[gas.species_index('CO2'), :]
        X_CO2 = f.X[gas.species_index('CO2'), :]

        conc_CO2_plt.append(conc_CO2)
        X_CO2_plt.append(X_CO2)

    # printing the plot of the temperature and velocity along the length of the domain
    plt.subplots()
    plt.subplot(2, 1, 1)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], T_plt[i], linewidth=1)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.ylabel('Temperature [K]', fontsize=12)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    plt.title('Temperature and velocity along the length of the domain for CH_4')
    plt.subplot(2, 1, 2)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], v_plt[i], linewidth=1)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.xlabel('Distance along the grid [m]', fontsize=12)
    plt.ylabel('Velocity [m/s]', fontsize=12)
    plt.tight_layout()

    # printing the plot of CO_2 variation and mole fraction along the length of the domain
    plt.subplots()
    plt.subplot(2, 1, 1)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], conc_CO2_plt[i], linewidth=1)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.ylabel('Concentration [kmol/m$^3$]', fontsize=9)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    plt.title('CO_2 variation along the length of the domain for CH_4')
    plt.subplot(2, 1, 2)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], X_CO2_plt[i], linewidth=1)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.xlabel('Distance along the grid [m]', fontsize=12)
    plt.ylabel('Mole fraction', fontsize=10)
    plt.tight_layout()

    plt.show()

# hydrogen

# giving the mole fractions of the elements in mixture
species_dict = {'H2': 1.8, 'O2': 1, 'N2': 3.76}

# calculating the temperature, the flame speed and the CO_2 variation along the length of the domain
for Pi in P0:
    x_plt = []
    T_plt = []
    v_plt = []
    conc_CO2_plt = []
    X_CO2_plt = []
    for Ti in T0:
        gas.TPX = Ti, Pi, species_dict
        f = ct.FreeFlame(gas, width=0.05)
        f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
        f.solve(loglevel=1, refine_grid=True, auto=False)
        x = f.grid
        T = f.T
        v = f.velocity
        x_plt.append(x)
        T_plt.append(T)
        v_plt.append(v)
        print('===========================================================================================')
        print('Mixture-averaged flamespeed = %.6f m/s' % (v[0]))
        print('===========================================================================================')
        conc_CO2 = f.concentrations[gas.species_index('CO2'), :]
        X_CO2 = f.X[gas.species_index('CO2'), :]

        conc_CO2_plt.append(conc_CO2)
        X_CO2_plt.append(X_CO2)

    # printing the plot of the temperature and velocity along the length of the domain
    plt.subplots()
    plt.subplot(2, 1, 1)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], T_plt[i], linewidth=1)
    # plt.legend(T0, loc="lower left")
    plt.grid('both', linestyle='--', linewidth=1)
    plt.ylabel('Temperature [K]', fontsize=12)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    plt.title('Temperature and velocity along the length of the domain for H_2')
    plt.subplot(2, 1, 2)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], v_plt[i], linewidth=1)
    # plt.legend(T0)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.xlabel('Distance along the grid [m]', fontsize=12)
    plt.ylabel('Velocity [m/s]', fontsize=12)
    plt.tight_layout()

    # printing the plot of CO_2 variation and mole fraction along the length of the domain
    plt.subplots()
    plt.subplot(2, 1, 1)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], conc_CO2_plt[i], linewidth=1)
    # plt.legend(T0)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.ylabel('Concentration [kmol/m$^3$]', fontsize=9)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    plt.title('CO_2 variation along the length of the domain for H_2')
    plt.subplot(2, 1, 2)
    for i in range(len(x_plt)):
        plt.plot(x_plt[i], X_CO2_plt[i], linewidth=1)
    # plt.legend(T0)
    plt.grid('both', linestyle='--', linewidth=1)
    plt.xlabel('Distance along the grid [m]', fontsize=12)
    plt.ylabel('Mole fraction', fontsize=10)
    plt.tight_layout()

    plt.show()

print("elapsed time:", time.process_time())
