# Sam Alvis, January 2022

import time
from common_funcs import *
import copy


def main():
    # define body
    Mars = Body()
    Mars.mu = 4.2828E4
    Mars.R = 3396.2
    Mars.J2 = 1.9555E-3
    Mars.omega = 2 * math.pi / 88619.6
    Mars.max_atm_height = 125.0
    Mars.load_density_info2("Mars_Atm_Data/0.txt")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    atms = np.zeros((2*5, 401))
    for idx in range(1, 6):
        body = copy.deepcopy(Mars)
        body.load_density_info2(f'Mars_Atm_Data/{idx}.txt')
        line, = ax.plot(body.alt_density_table[1, :], body.alt_density_table[0, :], 'b.')

    line, = ax.plot(Mars.alt_density_table[1, :], Mars.alt_density_table[0, :], 'k')

    ax.set_xscale('log')
    ax.set_xlabel("Log of Density (kg/m^3)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Sample Atmosphere Variations (Mars)")
    plt.show()

if __name__ == "__main__":
    main()
