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


def find_alpha(_thrust, _alpha_vect, _zp, _v, _gamma):
    """

    :param _thrust:
    :param _alpha_vect:
    :param _zp:
    :param _v:
    :param _gamma:
    :return:
    """
    T_aero = _thrust * np.cos(np.radians(_alpha_vect))
    drag = 0.5 * rho * S * f_cd(_alpha_vect) * v ** 2
    SFx = T_aero - drag - weight * np.sin(np.radians(gamma))

    _accel = SFx / mass

    df = pd.DataFrame([_alpha_vect, _accel, drag])
    df = df.T
    df.columns = ['alpha', 'accel', 'drag']
    df = df[df.accel > 0]
    ind_min = df.accel.idxmin()

    _min_accel = df.loc[ind_min, 'accel']

    _aoa = df.loc[ind_min, 'alpha']

    _drag = df.loc[ind_min, 'drag']

    return _min_accel, _aoa, _drag


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
v2min = max(1.13*v_stall, 1.1*vmca, vmuv2_lim)

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

vr += 0.075*v_stall
# define vef
vef = vr - 2


v = 0.
v_kt = uc.ms2kt(v)
t = 0.
x = 0.

dt = 0.125
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

while round(lift, 1) < round(weight, 1):

    aoa = q * dt + aoa0
    if aoa > aoa_max:

        print('Increase VR the a/c cannot Lift Off :(')
        increase = True
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
print(f'Velocity: {round(v_kt, 2)} kt')
print(f'Distance: {round(x, 2)} m')

vlof = v_kt
gamma = 0
vz = 0

teta0 = aoa
aoa0 = aoa
height = 0
print('\n====== Airborne Phase ======\n')
while height < 35.:

    _alpha_vect = np.linspace(config_data['alpha_min'], aoa0, 100)

    min_accel, aoa, _drag = find_alpha(thrust, _alpha_vect, zp, v, gamma)

    dv = min_accel * dt

    if aoa < aoa_max:
        teta = q * dt + teta0
    else:

        teta = teta0

    gamma = teta - aoa

    v_z = v * np.sin(np.radians(gamma))

    a_z = min_accel * np.sin(np.radians(gamma))

    dh = a_z*dt**2 + v_z * dt

    dx = dh * np.tan(np.radians(gamma))

    height += uc.m2ft(dh)

    # update
    x += dx
    v += dv
    t += dt


    teta0 = teta
    aoa0 = min(teta, aoa_max)
    v_kt = uc.ms2kt(v)
    if gamma > 2:
        print('error: gamma overshoots')
        break

    if v_kt >= v2min:
        print(f'V2min reached at {round(height, 2)}ft!')
        break
print(f'Height {height} ft')
print(f'Gamma {gamma} deg')
print(f'Teta {teta} deg')
print(f'Velocity {v_kt} kt\n')
print(f'Distance {x} m\n')
while height < 35.:

    teta = teta_target
    gamma = (thrust - drag) / weight
    aoa = teta - uc.rad2deg(gamma)
    lift = weight * np.cos(gamma)
    v_z = v * np.sin(gamma)
    dh = v_z * dt
    dx = v * np.cos(gamma) * dt

    # update
    x += dx
    t += dt
    height += uc.m2ft(dh)
    # print(f'Distance {x} m\n')

print(f'Height {round(height, 2)} ft')
print(f'Distance {round(x, 2)} m')

if uc.ms2kt(v) < v2min:
    print('Increase VR cannot reach V2min')
else:
    print(f'V2min reached with a VRmin {vr} kt')


cd_35ft = 2 * drag / (rho * S * v**2)

x_data = [config_data['cd'], config_data['cd_max']]
y_data = [config_data['alpha_min'], config_data['alpha_max']]
f_aoa = interp1d(x_data, y_data)
aoa_35ft = f_aoa(cd_35ft)

# Extra requirement from JAR 25.121(b)
gamma_jar_req = uc.perc2rad(2.4)
v_min_grad = (2 * (thrust * np.cos(np.radians(aoa_35ft)) - weight*np.sin(gamma_jar_req)) / (rho * S * cd_35ft)) ** .5

if v_min_grad > v:
    print('Increase VR. V2 does not reach the gradient requirement for ssg')
else:
    print(f'V2 reached sucessfully with a VR {vr}kt')










