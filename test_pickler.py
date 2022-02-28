# Sam Alvis, January 2022

import time
from common_funcs import *
import copy


def main():
    # define body
    Earth_0 = Body()
    Earth_0.mu = 3.9860E5
    Earth_0.R = 6371.0
    Earth_0.J2 = 1.9555E-3
    Earth_0.omega = 2 * math.pi / 88619.6
    Earth_0.max_atm_height = 200.0
    Earth_0.load_density_info("Atm_Data/0_ref.txt")

    # define craft
    SC_0 = Craft()
    SC_0.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC_0.beta_f = 150.0  # 140 by default?
    SC_0.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define params
    a_f_targ = Earth_0.R + 1000.0  # km

    # define x_0
    x_0 = [Earth_0.R + 200.0, 0, 0, 8, math.radians(-4.5), math.pi/2]

    sol_list = []
    for idx in range(11):
        # create objects (distinct slightly)
        Earth = copy.copy(Earth_0)
        SC = copy.copy(SC_0)
        Earth.mu = Earth_0.mu * (1.01 - idx / 100)

        # run sims (distinct slightly)
        solution = orbit_sim(x_0, Earth, SC, run_time=1024, time_res=1)
        sol_list.append(solution)

    # save objects
    save_outputs(sol_list, "abc4.pickle")


def main2():
    # load objects
    tic = time.perf_counter()
    loaded_list = load_outputs("abc.pickle")
    toc = time.perf_counter()
    print(toc-tic)

    # display objects
    for sol in loaded_list:
        plot_case(sol)


if __name__ == "__main__":
    main2()
