import json
import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

import unit_conversion as uc


def get_thrust(speed, coefs):
    """

    :param speed:
    :param coefs:
    :return:
    """

    speed_vector = np.array([1, speed, speed**2, speed**3, speed**4, speed**5, speed**6, speed**7], dtype='object')
    _coefs = np.array(coefs)

    _thrust = np.dot(speed_vector, _coefs)

    return _thrust


# get config data
with open('config.json', 'r') as f:
    config_data = json.load(f)


# inputs
mass = 80000.0          # [kg]
conf = '1+F'
zp = 0.                 # [ft]
rho = 1.225             # [kg/m3]
lg = 'Up'
engine_state = 'OEI'

# constants
g = 9.81                # [kg*m/s2]
S = config_data["S"]    # [m2]
mu = config_data["mu"]

weight = mass*g

# get speeds
vmca = config_data['vmca']
v_stall = config_data['v_stall']

# v2min definition
vmuv2_lim = 1.17*v_stall
# v_grad = uc.ms2kt(78.47452364358746)

v2min = max(1.13*v_stall, 1.1*vmca, vmuv2_lim)#, v_grad)


# get vR min
vr = max(1.05*v_stall, 1.05*vmca)


# aero data
aoa0 = config_data['alpha_min']
aoa_max = config_data['alpha_max']
cl_max = config_data['cl_max']
cl0 = config_data['cl']
cd_max = config_data['cd_max']
cd0 = config_data['cd']
engine_coefs = config_data['engine_coefs']

q = config_data['q']
teta_target = config_data['teta_target']

t_tr = 0


f_cl = interp1d(np.array([aoa0, aoa_max]),
                np.array([cl0, cl_max]))

f_cd = interp1d(np.array([aoa0, aoa_max]),
                np.array([cd_max, cd_max]))

vr += .08*v_stall
# define vef
vef = vr - 2

v = 0.
t = 0.
x = 0.
gamma = 0
gamma_rad = 0
vz = 0
height = 0

dt = 0.125
v_kt = uc.ms2kt(v)
vef_inst = True
print('\n====== Ground Phase ======\n')
while v_kt < vr:

    if vef_inst:
        if v_kt >= vef:
            thrust = thrust / 2  # + 1e4 thrust_oei#
            vef_inst = False
        else:
            thrust = get_thrust(v, engine_coefs)  # trhust_aeo + 1.8e4

    drag = 0.5 * rho * S * cd0 * v ** 2
    lift = 0.5 * rho * S * cl0 * v ** 2
    f_roll = abs(mu * (weight - lift))
    sf_x = thrust - drag - f_roll

    accel = sf_x / mass
    dv = accel * dt
    dx = 0.5 * accel * dt ** 2 + v * dt

    v += dv
    x += dx
    t += dt

    v_kt = uc.ms2kt(v)
print(f'Velocity: {round(v_kt, 2)} kt')
print(f'Distance: {round(x, 2)} m\n')

x_ground = x

print("\n====== Transition to VLOF ======\n")

t_tr = 0
transition_phase_passed = True
while round(lift, 1) < round(weight, 1):

    aoa = q * dt + aoa0
    if aoa > aoa_max:
        print('Increase VR the a/c cannot Lift Off :(')
        transition_phase_passed = False
        break
    else:
        drag = 0.5 * rho * S * f_cd(aoa) * v ** 2
        lift = 0.5 * rho * S * f_cl(aoa) * v ** 2
    f_roll = abs(mu * (lift - weight))

    sf_x = thrust * np.cos(np.radians(aoa)) - drag - f_roll  # aero axis

    accel = sf_x / mass
    dv = accel * dt
    dx = 0.5 * accel * dt ** 2 + v * dt

    v += dv
    x += dx
    t += dt
    t_tr += dt

    v_kt = uc.ms2kt(v)
    aoa0 = aoa

if transition_phase_passed:
    print(f'Velocity: {round(v_kt, 2)} kt')
    print(f'Distance: {round(x, 2)} m')
    print(f'Aoa: {round(aoa, 2)} deg')

    vlof = v_kt
    teta0 = aoa
    aoa0 = aoa
    airborne_phase_passed = True
    print('\n====== Airborne Phase ======\n')
    while height < 35.:

        T_aero = thrust * np.cos(np.radians(aoa_max))

        drag = 0.5 * rho * S * f_cd(aoa_max) * v ** 2
        lift = 0.5 * rho * S * f_cl(aoa_max) * v ** 2

        SFx = T_aero - drag - weight * np.sin(gamma_rad)

        accel = SFx / mass

        dgamma = (thrust * np.sin(np.radians(aoa_max)) + lift - weight) / (mass * v) * dt

        gamma_rad += dgamma
        gamma = np.degrees(gamma_rad)
        teta = aoa_max + gamma

        v_z = v * np.sin(gamma_rad)

        dv = accel * dt

        dh = v * np.sin(gamma_rad) * dt  # a_z*dt**2 + v_z * dt

        dx = v * np.cos(gamma_rad) * dt  # dh * np.tan(gamma_rad)

        # update
        x += dx
        v += dv
        t += dt
        height += uc.m2ft(dh)
        teta0 = teta
        v_kt = uc.ms2kt(v)

        if gamma > 3:
            print('error: gamma overshoots')
            airborne_phase_passed = False
            break

        if v_kt >= v2min:
            print(f'V2min reached at {round(height, 2)}ft!')
            break

    if airborne_phase_passed:
        print(f'Height {height} ft')
        print(f'Gamma {gamma} deg')
        print(f'Teta {teta} deg')
        print(f'Velocity {v_kt} kt\n')
        print(f'Distance {x} m\n')
        while height < 35.:

            teta = teta_target
            aoa = teta - gamma
            drag = 0.5 * rho * S * f_cd(aoa) * v ** 2
            gamma = (thrust - drag) / weight

            lift = weight * np.cos(gamma)
            v_z = v * np.sin(gamma)
            dh = v_z * dt
            dx = v * np.cos(gamma) * dt

            # update
            x += dx
            t += dt
            height += uc.m2ft(dh)

        print(f'Height {round(height, 2)} ft')
        print(f'Distance {round(x, 2)} m')

        if uc.ms2kt(v) < v2min:
            print('Increase VR cannot reach V2min')
        else:
            print(f'V2min reached with a VRmin {vr} kt')


        cd_35ft = round(2 * drag / (rho * S * v**2), 3)

        x_data = [config_data['cd'], config_data['cd_max']]
        y_data = [config_data['alpha_min'], config_data['alpha_max']]
        f_aoa = interp1d(x_data, y_data)
        aoa_35ft = f_aoa(cd_35ft)

        # TODO create a loop in order fro ac speed to converge to the requirement(easy)
        # Extra requirement from JAR 25.121(b)
        gamma_jar_req = uc.perc2rad(2.4)
        v_min_grad = (2 * (thrust * np.cos(np.radians(aoa_35ft)) - weight*np.sin(gamma_jar_req)) / (rho * S * cd_35ft)) ** .5

        if v_min_grad > v:
            print('Increase VR. V2 does not reach the gradient requirement for ssg')
        else:
            print(f'V2 reached sucessfully with a VR {vr}kt')










