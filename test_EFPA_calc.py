# Sam Alvis, January 2022
import math
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

    # define craft
    SC = Craft()
    SC.beta_i = 30.0
    SC.beta_f = 135.0
    SC.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    a_f_targ = Mars.R + 500.0  # km

    # define x_0
    x_0 = [Mars.R + 125.0, math.radians(18.44), math.radians(77.45), 6, 0, math.radians(45)]  # NOTE: 0 will be replaced by optimal EFPA

    # define target altitude and run function, then plot results
    tic = time.perf_counter()
    solution = EFPA_calc(a_f_targ, x_0, Mars, SC)
    toc = time.perf_counter()
    print(f"EFPA = {math.degrees(solution.craft_obj.x_0[4])}; t_P_max = {solution.t[np.argmax(solution.P_dyn)]}; t_eject = {solution.craft_obj.t_eject} [took {toc-tic:0.4f} seconds]")


if __name__ == "__main__":
    main()
