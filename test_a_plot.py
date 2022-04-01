# Sam Alvis, January 2022

import time

import matplotlib.pyplot as plt

from common_funcs import *
from f_analysis_2 import *


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

    # define target apoapsis
    a_f_targ = 500 + Mars.R

    # define x_0
    SC.x_0 = [Mars.R + 125.0, math.radians(18.44), math.radians(77.45), 6, 0, math.radians(45)]
    # SC.x_0[4] = EFPA_calc(a_f_targ, x_0, Mars, SC).craft_obj.x_0[4]
    SC.x_0[4] = math.radians(-9.875)

    # load linear distributions
    # lin_sols = gen_lin_atm(0.01, 0.75, 1.25, "test5", a_f_targ, Mars, SC)
    lin_sols = load_lin_atm("test1_fix")

    # run through linear variations and plot a/a_dot
    a_arr = []
    a_dot_arr = []
    for lin_sol in lin_sols:
        print(lin_sol.craft_obj.t_eject, lin_sol.body_obj.scale_factor)
        t_step_star = lin_sol.t[1] - lin_sol.t[0]
        a_arr.append(lin_sol.a_t[int(math.floor(lin_sol.craft_obj.t_eject / t_step_star))])
        a_dot_arr.append(lin_sol.a_dot_t[int(lin_sol.craft_obj.t_eject)])

    print(a_arr)
    print(a_dot_arr)
    plt.scatter(a_arr, a_dot_arr)
    for idx in range(len(a_arr)):
        plt.annotate(f"{idx}", (a_arr[idx]+2, a_dot_arr[idx]-0.02))

    plt.xlabel("a")
    plt.ylabel("a_dot")
    plt.show()

if __name__ == "__main__":
    main()
