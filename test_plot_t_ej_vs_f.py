# Sam Alvis, January 2022

from common_funcs import *

def main():
    """Description"""

    # define body
    Mars = Body()
    Mars.mu = 4.2828E4
    Mars.R = 3396.2
    Mars.J2 = 1.9555E-3
    Mars.omega = 2 * math.pi / 88619.6
    Mars.max_atm_height = 125.0
    Mars.load_density_info("mars_atm_data.csv")

    # define craft
    SC = Craft()
    SC.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC.beta_f = 150.0  # 140 by default?
    SC.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    a_f_targ = Mars.R + 1000.0  # km

    plot_t_eject_vs_f(a_f_targ, Mars, SC)


if __name__ == "__main__":
    main()
