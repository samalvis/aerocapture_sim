# Sam Alvis, January 2022

import math
import numpy as np
from scipy import integrate
from scipy import optimize as opt
import matplotlib.pyplot as plt
import time
import pickle
import copy


# CLASSES ##############################################################################################################

class Body:
    """Description"""
    mu = 0
    R = 0
    J2 = 0
    omega = 0  # TODO: rads/s
    max_atm_height = 0
    scale_height = 0
    scale_density = 0
    scale_factor = 1
    alt_density_table = None

    def density_interp(self, r):
        """Description"""
        alt = np.linalg.norm(r) - self.R
        if alt > self.max_atm_height:
            return 0

        alt_table = self.alt_density_table[0, :]
        density_table = self.alt_density_table[1, :]

        idx = 0
        while alt_table[idx + 1] < alt and idx + 1 < len(alt_table) - 1:
            idx = idx + 1
        alt_low = alt_table[idx]
        alt_high = alt_table[idx + 1]
        density_low = density_table[idx]
        density_high = density_table[idx + 1]

        density = density_low - (alt - alt_low) * (density_high - density_low) / (alt_low - alt_high)
        return self.scale_factor * density

    def density_exp(self, r):
        """Description"""
        alt = np.linalg.norm(r) - self.R
        if alt > self.max_atm_height:
            return 0

        density = self.scale_density * math.exp(-alt / self.scale_height)
        return self.scale_factor * density

    def load_density_info(self, file_name):
        """Description"""
        data = np.loadtxt(file_name, skiprows=1, delimiter=' ')

        alt_table = data[:, 1]
        density_table = data[:, 3]
        self.alt_density_table = np.vstack((alt_table, density_table))

        self.scale_height = np.mean(data[:, 3])  # NOTE: is this reasonable?
        self.scale_density = density_table[0]

        return

    def load_density_info2(self, file_name):
        """Description"""
        data = np.loadtxt(file_name, skiprows=1, delimiter=' ')

        alt_table = data[:, 1]
        density_table = np.multiply(data[:, 2], data[:, 3])
        self.alt_density_table = np.vstack((alt_table, density_table))

        self.scale_height = np.mean(data[:, 2])  # NOTE: is this reasonable?
        self.scale_density = density_table[0]

        return


class Craft:
    """Description"""
    t_eject = 0
    beta_i = 0  # TODO: UNITS????
    beta_f = 0  # TODO: make beta a function with t as an input! simplifies behavior in other parts
    x_0 = []


class Output:
    """Description"""
    body_obj = Body()
    craft_obj = Craft()
    t = []
    r = []
    v = []
    h = []
    E = []
    a_t = []
    a_dot_t = []
    rho = 0
    P_dyn = 0
    hasImpacted = 0
    hasCaptured = 0


def prep_output_file(filename):
    """"""
    file = open(filename, "wb")  # TODO: should probably have a check if it exists or something along those lines
    file.close()


def save_output(output, filename):  # TODO: literally zero security here lol
    """"""
    file = open(filename, "ab")
    pickle.dump(output, file)
    file.close()


def load_outputs(filename):
    """"""
    file = open(filename, "rb")
    loaded_list = []
    while 1:
        try:
            loaded_list.append(pickle.load(file))
        except EOFError:
            break
    file.close()
    return np.array(loaded_list).flatten()


# Rotation Matrix Functions ############################################################################################

def M1(theta):
    return [[1.0, 0.0, 0.0], [0.0, math.cos(theta), math.sin(theta)], [0.0, -math.sin(theta), math.cos(theta)]]


def M2(theta):
    return [[math.cos(theta), 0.0, -math.sin(theta)], [0.0, 1.0, 0.0], [math.sin(theta), 0.0, math.cos(theta)]]


def M3(theta):
    return [[math.cos(theta), math.sin(theta), 0.0], [-math.sin(theta), math.cos(theta), 0.0], [0.0, 0.0, 1.0]]


# Coordinate Transformation Functions ##################################################################################

def car_2_loc(x_car, body_obj, t):
    """Description"""

    x_N = x_car[0:3]
    v_N = x_car[3:6]

    # Inertial to Planet-Fixed
    omega_p = body_obj.omega
    Omega = omega_p * t
    IN = M3(Omega)
    x_I = np.matmul(IN, x_N)
    v_I = np.matmul(IN, v_N)

    # Planet-Fixed to Positional
    theta = math.atan2(x_I[1], x_I[0])
    phi = math.asin(x_I[2] / np.linalg.norm(x_I))
    EI = np.matmul(M2(-phi), M3(theta))
    x_E = np.matmul(EI, x_I)
    v_E = np.matmul(EI, v_I)

    w_IN_E = np.matmul(EI, np.matmul(IN, [0, 0, omega_p]))

    v_E = v_E - np.cross(w_IN_E, x_E)

    # Positional to Velocity
    gamma = math.asin(v_E[0] / np.linalg.norm(v_E))
    psi = math.atan2(v_E[1], v_E[2])

    r = np.linalg.norm(x_E)
    U = np.linalg.norm(v_E)

    return [r, theta, phi, U, gamma, psi]


def loc_2_car(x_loc, body_obj, t):
    """Description"""

    r = x_loc[0]
    theta = x_loc[1]
    phi = x_loc[2]
    U = x_loc[3]
    gamma = x_loc[4]
    psi = x_loc[5]

    # Velocity to Positional
    v_S = [0, 0, U]
    ES = np.transpose(np.matmul(M2(gamma), M1(-psi)))
    v_E = np.matmul(ES, v_S)

    # Positional to Planet-Fixed
    x_E = [r, 0, 0]
    IE = np.transpose(np.matmul(M2(-phi), M3(theta)))
    x_I = np.matmul(IE, x_E)
    v_I = np.matmul(IE, v_E)

    # Planet Fixed to Inertial
    omega_p = body_obj.omega
    Omega = omega_p * t
    NI = np.transpose(M3(Omega))
    x_N = np.matmul(NI, x_I)
    v_N = np.matmul(NI, v_I)

    w_IN_N = [0, 0, omega_p]

    v_N_inertial = v_N + np.cross(w_IN_N, x_N)

    return [x_N[0], x_N[1], x_N[2], v_N_inertial[0], v_N_inertial[1], v_N_inertial[2]]


# Events ###############################################################################################################

def impact(t, x, body_obj):
    """Description"""

    # Unpack variables
    R = body_obj.R
    r = x[0]

    return R - r


# Disturbance Functions ################################################################################################
# NOTE: all disturbance calculated in E frame, may be easier to calculate in s frame for actual applications.
# TODO: investigate best frame further.

def dist_J2(t, x, body_obj, craft_obj):
    """Description"""

    # Unpack variables
    mu = body_obj.mu
    J2 = body_obj.J2
    R = body_obj.R
    r = x[0]
    phi = x[2]
    gamma = x[4]
    psi = x[5]

    # Calculate acceleration due to J2
    J2_term = 3.0 * J2 * (R**2.0) / (2.0 * (r**2.0))
    a_J2_E = np.multiply(-mu / r**2.0 * J2_term, [1.0 - 3.0 * (math.sin(phi)**2.0), 0, 2.0 * math.sin(phi) * math.cos(phi)])

    return a_J2_E


def dist_drag(t, x, body_obj, craft_obj):
    """Description"""

    # Unpack variables
    r = x[0]
    U = x[3]
    gamma = x[4]
    psi = x[5]

    if t >= craft_obj.t_eject:
        beta = craft_obj.beta_f
    else:
        beta = craft_obj.beta_i

    # Get density at specified altitude
    rho = body_obj.density_interp(r)

    # Calculate acceleration due to atmospheric drag
    a_drag_S = np.multiply(1000.0, [0.0, 0.0, -(rho * U**2.0 / 2.0 / beta)])  # NOTE: converts meters to kilometers, keeping units in km/s^2

    # Convert to S frame for consistency
    SE = np.matmul(M2(gamma), M1(-psi))  # TODO: Will need to convert this into E fram accounting for rotation of E with respect to S!
    ES = np.transpose(SE)
    a_drag_E = np.matmul(ES, a_drag_S)

    return a_drag_E


# 2 Body Problem ODE Function ##########################################################################################

def ode_two_body(t, x, body_obj, craft_obj, dist_funcs):
    """Description"""

    # Unpack variables
    r = x[0]
    theta = x[1]
    phi = x[2]
    U = x[3]
    gamma = x[4]
    psi = x[5]

    mu = body_obj.mu
    omega_p = body_obj.omega

    # Prestate DCMs
    EI = np.matmul(M2(-phi), M3(theta))
    SE = np.matmul(M2(gamma), M1(-psi))
    ES = np.transpose(SE)

    # Create vectors
    R_E = [r, 0.0, 0.0]
    R_S = np.matmul(SE, R_E)
    U_S = [0.0, 0.0, U]
    U_E = np.matmul(ES, U_S)

    # Calculate velocity based changes  # TODO: confirm this is correct (b/c it's already atmosphere relative? Still doesn't account for cross in my mind)
    r_dot = U * math.sin(gamma) #U_E[0]
    theta_dot = U * math.cos(gamma) * math.sin(psi) / r / math.cos(phi) #U_E[1] / (r * math.cos(phi))  # NOTE: could beautify this a bit by defining the r vector, then dividing component by component...
    phi_dot = U * math.cos(gamma) * math.cos(psi) / r #U_E[2] / r

    # Gather all accelerations, starting with 2 body
    U_dot_E_inertial = np.array([-mu / r**2.0, 0.0, 0.0])
    for dist_func in dist_funcs:
        U_dot_E_inertial = U_dot_E_inertial + dist_func(t, x, body_obj, craft_obj)

    # Create rotation vectors
    w_IN_I = [0.0, 0.0, omega_p]
    w_IN_E = np.matmul(EI, w_IN_I)
    w_EI_I = [phi_dot * math.sin(theta), -phi_dot * math.cos(theta), theta_dot]  # TODO: update from doc
    w_EI_E = np.matmul(EI, w_EI_I)
    w_EN_E = w_EI_E + w_IN_E
    w_EN_S = np.matmul(SE, w_EN_E)

    # Create in between vector, from derivatives of inertial velocity and velocity due to rotation rate
    temp_val_1 = r_dot * omega_p * math.cos(phi) - r * phi_dot * omega_p * math.sin(phi)
    temp_val_2 = r * omega_p * math.cos(phi)
    temp_vec_E = U_dot_E_inertial - [-temp_val_2 * w_EN_E[2], temp_val_1, temp_val_2 * w_EN_E[0]]
    temp_vec_S = np.matmul(SE, temp_vec_E)

    # Calculate acceleration based changes
    U_dot = temp_vec_S[2]
    gamma_dot = (temp_vec_S[0] / U) - w_EN_S[1]
    psi_dot = ((temp_vec_S[1] / -U) - w_EN_S[0]) / (-1 * math.cos(gamma))

    x_dot = np.zeros([6, ])
    x_dot[0] = r_dot
    x_dot[1] = theta_dot
    x_dot[2] = phi_dot
    x_dot[3] = U_dot
    x_dot[4] = gamma_dot
    x_dot[5] = psi_dot

    return x_dot


# ODE Wrapper Function #################################################################################################

def orbit_sim(x_0, body_obj, craft_obj, run_time=3600., time_res=60., method='DOP853'):
    """Description"""

    # Calculate Needed Variables
    num_time_points = int(run_time / time_res) + 1

    # Prepare ODE Variables
    times = np.linspace(0, run_time, num_time_points)   # NOTE: this is the number of points, not the step between points!
    dist_funcs = [dist_J2, dist_drag]
    ode_fun_two_body = lambda t, x: ode_two_body(t, x, body_obj, craft_obj, dist_funcs)
    lam_impact = lambda t, x: impact(t, x, body_obj)
    lam_impact.terminal = True

    # Run ODE
    solution = integrate.solve_ivp(ode_fun_two_body, (0, run_time), x_0, method=method, t_eval=times,
                                   rtol=1E-10, atol=1E-10, events=lam_impact)  # TODO: restore DOP853?

    # Analyze ODE Output
    t_arr = solution.t
    r_arr = np.zeros((3, len(t_arr)))
    v_arr = np.zeros((3, len(t_arr)))
    for idx in range(len(t_arr)):
        x_car = loc_2_car(solution.y[:, idx], body_obj, t_arr[idx])
        r_arr[:, idx] = x_car[0:3]
        v_arr[:, idx] = x_car[3:6]

    h_norm_arr = np.zeros(len(t_arr))
    E_arr = np.zeros(len(t_arr))
    rho_arr = np.zeros(len(t_arr))
    P_dyn_arr = np.zeros(len(t_arr))
    for idx in range(len(t_arr)):
        r_norm = np.linalg.norm(r_arr[:, idx])
        v_norm = np.linalg.norm(v_arr[:, idx])
        h_norm_arr[idx] = np.linalg.norm(np.cross(r_arr[:, idx], v_arr[:, idx]))
        E_arr[idx] = 0.5*(v_norm**2.) - body_obj.mu / r_norm
        rho_arr[idx] = body_obj.density_interp([r_norm])
        P_dyn_arr[idx] = 0.5*rho_arr[idx]*((v_norm*1000.)**2.)

    hasImpacted = (solution.status == 1)
    hasCaptured = (E_arr[-1] < 0)

    # Return
    ans = Output()
    ans.body_obj = body_obj
    ans.craft_obj = craft_obj
    ans.craft_obj.x_0 = x_0  # TODO: this should happen in the set up, make more streamlined
    ans.t = t_arr
    ans.r = r_arr
    ans.v = v_arr
    ans.h = h_norm_arr
    ans.E = E_arr
    ans.rho = rho_arr
    ans.P_dyn = P_dyn_arr
    ans.hasImpacted = hasImpacted
    ans.hasCaptured = hasCaptured
    # TODO: maybe should the initial orbital conditions be recorded? could be useful for plotting
    e = np.sqrt(1 + 2 * np.multiply(ans.E, np.power(h_norm_arr, 2)) / (body_obj.mu**2))
    a = np.divide(-body_obj.mu, 2*ans.E)
    ans.a_t = np.multiply(a, 1+e)

    # TODO: this loop structure is awful, clean up and speed up everything plz
    # repalce all negative ap. with infinity, replace all less than R_body with 0?
    for idx in range(len(ans.E)):
        if ans.a_t[idx] < 0:
            ans.a_t[idx] = float('inf')
        if ans.a_t[idx] < body_obj.R:
            ans.a_t[idx] = 0

    ans.a_dot_t = np.zeros(len(ans.a_t))
    for idx in range(len(ans.a_t) - 1):  # NOTE: rough linear approximation of a_dot!
        ans.a_dot_t[idx + 1] = (ans.a_t[idx + 1] - ans.a_t[idx]) / (ans.t[idx + 1] - ans.a_t[idx])

    return ans


# Post-Simulation Analysis #############################################################################################

def plot_case(ans):  # TODO: unable to format set_ylabel to get titles nice (keep overlapping w/ previous subplots)
    """Description"""

    # Unpack
    body_obj = ans.body_obj
    craft_obj = ans.craft_obj
    t_samp = ans.t
    r_samp = ans.r
    h_samp = ans.h
    E_samp = ans.E
    rho_samp = ans.rho
    P_dyn_samp = ans.P_dyn
    t_eject = craft_obj.t_eject

    # Plot the Orbit
    plt.plot(r_samp[0, :], r_samp[1, :], 'k-')
    fig1 = plt.gcf()
    ax = fig1.gca()
    ax.add_patch(plt.Circle((0, 0), body_obj.R, color='r'))
    ax.add_patch(plt.Circle((0, 0), body_obj.R + body_obj.max_atm_height, color='b', fill=False))
    ax.axis('equal')
    plt.title('Calculated Trajectory of Craft')

    # Plot Time Variables
    fig2, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 10))

    # NOTE: this method of finding altitude ignores terrain, etc.
    ax1.plot(t_samp, np.linalg.norm(r_samp, axis=0)-body_obj.R)  # TODO: confirm this is correct axis!!!
    ax1.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax1.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    ax1.minorticks_on()
    ax1.tick_params('x', labelbottom=False)
    ax1.xaxis.set_ticks(np.arange(t_samp[0], t_samp[-1]+1, max(((t_samp[-1]-t_samp[0])//(60*6)*60), 60)))
    ax1.axhline(0, 0, 1, color='black')
    ax1.set_ylabel("Altitude (km)")

    ax2.plot(t_samp, E_samp)
    ax2.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax2.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    ax2.minorticks_on()
    ax2.tick_params('x', labelbottom=False)
    ax2.xaxis.set_ticks(np.arange(t_samp[0], t_samp[-1]+1, max(((t_samp[-1]-t_samp[0])//(60*6)*60), 60)))
    ax2.axhline(0, 0, 1, color='black')
    ax2.set_ylabel("Orbital Energy (MJ/kg)")

    ax3.plot(t_samp, rho_samp)
    ax3.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax3.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    ax3.minorticks_on()
    ax3.tick_params('x', labelbottom=False)
    ax3.xaxis.set_ticks(np.arange(t_samp[0], t_samp[-1] + 1, max(((t_samp[-1] - t_samp[0]) // (60 * 6) * 60), 60)))
    ax3.set_ylabel("Atm. Density (kg/m^3)")

    ax4.plot(t_samp, P_dyn_samp)
    ax4.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax4.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    ax4.minorticks_on()
    ax4.xaxis.set_ticks(np.arange(t_samp[0], t_samp[-1]+1, max(((t_samp[-1]-t_samp[0])//(60*6)*60), 60)))
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Dynamic Pressure (Pa)")

    # Add lines to indicate eject time
    ax1.axvline(x = t_eject, color = 'r', label = 'axvline - full height')
    ax2.axvline(x = t_eject, color = 'r', label = 'axvline - full height')
    ax3.axvline(x = t_eject, color = 'r', label = 'axvline - full height')
    ax4.axvline(x = t_eject, color = 'r', label = 'axvline - full height')

    if ans.hasImpacted:
        fig2.suptitle(f"U_0 = {round(ans.craft_obj.x_0[3], 3)} km/s; gamma_0 = {round(math.degrees(ans.craft_obj.x_0[4]), 3)} deg; Status = Impacted")
    elif ans.hasCaptured:
        fig2.suptitle(f"U_0 = {round(ans.craft_obj.x_0[3], 3)} km/s; gamma_0 = {round(math.degrees(ans.craft_obj.x_0[4]), 3)} deg; Status = Captured; ap_f = {round(ans.a_t[-1] - body_obj.R, 1)} km alt")
    else:
        fig2.suptitle(f"U_0 = {round(ans.craft_obj.x_0[3], 3)} km/s; gamma_0 = {round(math.degrees(ans.craft_obj.x_0[4]), 3)} deg; Status = Escaped")

    plt.show()


# Eject Time Calculator ################################################################################################

def eject_calc_2(a_f_targ, x_0, body_obj, craft_obj):
    """Description"""

    # define
    run_time = 1024.  # NOTE: maybe these should be passed, or defined with an option to pass
    time_res = 1.

    # check bounds and confirm behavior falls in range
    craft_obj.t_eject = 0
    solution_L = orbit_sim(x_0, body_obj, craft_obj, run_time=run_time, time_res=time_res, method='RK45')
    a_f_max = solution_L.a_t[-1]

    if a_f_max < a_f_targ:
        return solution_L

    craft_obj.t_eject = run_time
    solution_U = orbit_sim(x_0, body_obj, craft_obj, run_time=run_time, time_res=time_res, method='RK45')
    a_f_min = solution_U.a_t[-1]

    if a_f_min > a_f_targ:
        return solution_U

    upper_bound = run_time
    lower_bound = 0
    err = 100
    err_tol = 1e-6
    solution = 0
    count = 0
    left_infinite_flag = math.isinf(a_f_max)
    while (abs(err) > err_tol) & (count < 10):  # TODO: make this count a parameter!
        # if left infinite, bisect!
        if left_infinite_flag:
            craft_obj.t_eject = (upper_bound + lower_bound) / 2.
            print(f"l_inf, {craft_obj.t_eject}")
        # else, linearize!
        else:
            inv_slope = -(upper_bound - lower_bound) / (solution_U.a_t[-1] - solution_L.a_t[-1])
            craft_obj.t_eject = (solution_L.a_t[-1] - a_f_targ) * inv_slope + lower_bound
            print(f"l_fin, {craft_obj.t_eject}")

        # try guess
        solution = orbit_sim(x_0, body_obj, craft_obj, run_time=run_time, time_res=time_res, method='RK45')
        a_f_norm = solution.a_t[-1] - a_f_targ

        left_infinite_flag = math.isinf(a_f_norm)

        # if positive, set to left bound
        if a_f_norm > 0:
            solution_L = solution
            lower_bound = craft_obj.t_eject
        # if negative, set to right bound
        else:
            solution_U = solution
            upper_bound = craft_obj.t_eject

        # set error and repeat
        err = a_f_targ - solution.a_t[-1]
        count = count + 1

    return solution


def eject_calc(a_f_targ, x_0, Mars, SC):
    """Description"""

    # define
    run_time = 1024  # NOTE: maybe these should be passed, or defined with an option to pass
    time_res = 1

    # check bounds and confirm behavior falls in range
    SC.t_eject = 0
    solution1 = orbit_sim(x_0, Mars, SC, run_time=run_time, time_res=time_res, method='RK45')
    a_f_max = solution1.a_t[-1]

    if a_f_max < a_f_targ:
        return solution1

    SC.t_eject = run_time
    solution2 = orbit_sim(x_0, Mars, SC, run_time=run_time, time_res=time_res, method='RK45')
    a_f_min = solution2.a_t[-1]

    if a_f_min > a_f_targ:
        return solution2

    # if behavior in range, root solve
    SC.t_eject = run_time / 2
    upper = run_time
    lower = 0
    err = 100
    err_tol = 0.01  # acceptable error (km) from target apoapsis
                    # NOTE: must be small so that time error converges consistently! (?)
    solution = 0
    count = 0
    while (abs(err) > err_tol) and (count <= 25):
        solution = orbit_sim(x_0, Mars, SC, run_time=run_time, time_res=time_res, method='RK45')
        err = a_f_targ - solution.a_t[-1]
        if err < 0:
            # higher than desired: increase time?
            lower = SC.t_eject
            SC.t_eject = (SC.t_eject + upper)/2
        else:
            # lower than desired: decrease time?
            upper = SC.t_eject
            SC.t_eject = (SC.t_eject + lower)/2
        count = count + 1

    return solution


def EFPA_calc(a_f_targ, x_0, body_obj, craft_obj):
    """Description"""
    # start EFPA at zero
    EFPA = 0
    x_0[4] = math.radians(EFPA)
    solution1 = copy.deepcopy(eject_calc(a_f_targ, x_0, body_obj, craft_obj))
    solution2 = 0

    # step down by degrees until t_j before max dynamic pressure (what about imapct?)
    while EFPA > -90:  # TODO: make a more stable escape condition
        EFPA = EFPA-1
        x_0[4] = math.radians(EFPA)
        solution2 = copy.deepcopy(eject_calc(a_f_targ, x_0, body_obj, craft_obj))
        if solution2.craft_obj.t_eject < solution2.t[np.argmax(solution2.P_dyn)]:
            break

        solution1 = copy.deepcopy(solution2)

    # bisect inward to find solution
    err = 1000
    err_max = 15  # within 15 s (parameter)  # TODO: explore weird behavior here?
    solution = 0
    count = 0
    while (abs(err) > err_max) and (count <= 6):
        x_0[4] = (solution1.craft_obj.x_0[4] + solution2.craft_obj.x_0[4]) / 2
        solution = eject_calc(a_f_targ, x_0, body_obj, craft_obj)
        err = solution.craft_obj.t_eject - solution.t[np.argmax(solution.P_dyn)]
        if err < 0:
            solution2 = copy.deepcopy(solution)
        else:
            solution1 = copy.deepcopy(solution)
        count = count + 1

    return solution


def map_t_eject_vs_U_0_and_gamma_0(a_f_targ, body_obj, craft_obj):
    """Description"""

    res = 20
    U_0_min = 5
    U_0_max = 6
    gamma_0_min = math.radians(-8)
    gamma_0_max = math.radians(-12)
    U_0_arr = np.linspace(U_0_min, U_0_max, res)
    gamma_0_arr = np.linspace(gamma_0_min, gamma_0_max, res)

    x_0 = [body_obj.R + 125.0, 0, 0, 0, 0, math.pi/2]  # NOTE: some parameters hardcoded, but U and gamma changed by sim

    t_ej_arr = np.zeros((res, res))  # TODO: could make different resolutions if needed?
    for idx_U_0, U_0 in enumerate(U_0_arr):
        for idx_gamma_0, gamma_0 in enumerate(gamma_0_arr):
            x_0[3] = U_0
            x_0[4] = gamma_0
            t_ej_arr[idx_U_0][idx_gamma_0] = eject_calc(a_f_targ, x_0, body_obj, craft_obj).craft_obj.t_eject  # TODO: figure out how to save objects and effectivley pull attributes into 2d array

    fig, ax = plt.subplots()

    # TODO: resolve this error that pops up
    c = ax.pcolormesh(U_0_arr, (gamma_0_arr) * (180 / math.pi), t_ej_arr, cmap='seismic')
    fig.colorbar(c, ax=ax)
    fig.suptitle("Eject Time Sensitivity [s] vs Entry Velocity and Entry Angle")
    ax.set_xlabel("Entry Velocity [km/s]")
    ax.set_ylabel("Entry Angle [degrees]")
    plt.show()


def plot_t_eject_vs_f(a_f_targ, body_obj, craft_obj):
    """Description"""

    f_min = 0.85
    f_max = 1.1
    res = 1000
    f_arr = np.linspace(f_min, f_max, res)

    scale_density_0 = body_obj.scale_density

    x_0 = [body_obj.R + 125.0, 0, 0, 5.25, math.radians(-10), math.pi/2]

    t_ej_arr = np.zeros(res)  # TODO: could make different resolutions if needed?
    for idx_f, f in enumerate(f_arr):
        body_obj.scale_density = scale_density_0 * f
        t_ej_arr[idx_f] = eject_calc(a_f_targ, x_0, body_obj, craft_obj).craft_obj.t_eject  # TODO: figure out how to save objects and effectivley pull attributes into 2d array

    plt.plot(f_arr, t_ej_arr)
    plt.show()
