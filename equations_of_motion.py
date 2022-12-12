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


def find_alpha(_thrust, _alpha_vect, _v, _gamma):
    """

    :param _thrust:
    :param _alpha_vect:
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

vr += .09*v_stall
# define vef
vef = vr - 2

v = 0.
t = 0.
x = 0.
gamma = 0
gamma_rad = 0
v_z = 0
height = 0

v_kt_log = [0]
t_log = [0]
x_log = [0]
gamma_log = [0]
vz_log = [0]
height_log = [0]
thrust_log = [0]
lift_log = [0]
drag_log = [0]
teta_log = [0]
alpha_log = [0]


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

    v_kt_log.append(v_kt)
    t_log.append(t)
    x_log.append(x)
    gamma_log.append(gamma)
    vz_log.append(v_z)
    height_log.append(height)
    thrust_log.append(thrust)
    lift_log.append(lift)
    drag_log.append(drag)
    teta_log.append(0)
    alpha_log.append(0)


print(f'Velocity: {round(v_kt, 2)} kt')
print(f'Distance: {round(x, 2)} m\n')

x_ground = x

print("\n====== Transition to VLOF ======\n")

t_tr = 0

while round(lift, 1) < round(weight, 1):

    aoa = q * dt + aoa0
    assert aoa < aoa_max, 'Increase VR the a/c cannot Lift Off :('

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

    v_kt_log.append(v_kt)
    t_log.append(t)
    x_log.append(x)
    gamma_log.append(gamma)
    vz_log.append(v_z)
    height_log.append(height)
    thrust_log.append(thrust)
    lift_log.append(lift)
    drag_log.append(drag)
    teta_log.append(aoa)
    alpha_log.append(aoa)



print(f'Velocity: {round(v_kt, 2)} kt')
print(f'Distance: {round(x, 2)} m')
print(f'Aoa: {round(aoa, 2)} deg')

vlof = v_kt

alpha_vector = np.linspace(config_data['alpha_min'], config_data['alpha_max'], 100, endpoint=True)

print('\n====== Airborne Phase ======\n')
while height < 35.:

    accel, aoa, drag = find_alpha(thrust, alpha_vector, v, gamma)

    lift = 0.5 * rho * S * f_cl(aoa) * v ** 2

    dgamma = (thrust * np.sin(np.radians(aoa)) + lift - weight) / (mass * v) * dt

    gamma_rad += dgamma
    gamma = np.degrees(gamma_rad)
    teta = aoa + gamma

    v_z = v * np.sin(gamma_rad)

    dv = accel * dt

    dh = v * np.sin(gamma_rad) * dt  # a_z*dt**2 + v_z * dt

    dx = v * np.cos(gamma_rad) * dt  # dh * np.tan(gamma_rad)

    # update
    x += dx
    v += dv
    t += dt
    height += uc.m2ft(dh)
    v_kt = uc.ms2kt(v)

    assert gamma < 3, 'Gamma overshoots. Increase VR'

    v_kt_log.append(v_kt)
    t_log.append(t)
    x_log.append(x)
    gamma_log.append(gamma)
    vz_log.append(v_z)
    height_log.append(height)
    thrust_log.append(thrust)
    lift_log.append(lift)
    drag_log.append(drag)
    teta_log.append(teta)
    alpha_log.append(aoa)

    if v_kt >= v2min:
        print(f'V2min reached at {round(height, 2)}ft!')
        break

aoa_v2 = aoa
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

    v_kt_log.append(v_kt)
    t_log.append(t)
    x_log.append(x)
    gamma_log.append(gamma)
    vz_log.append(v_z)
    height_log.append(height)
    thrust_log.append(thrust)
    lift_log.append(lift)
    drag_log.append(drag)
    teta_log.append(teta)
    alpha_log.append(aoa)

print(f'Height {round(height, 2)} ft')
print(f'Distance {round(x, 2)} m')

df = pd.DataFrame([t_log, x_log, v_kt_log, thrust_log, lift_log, drag_log, alpha_log, teta_log, gamma_log, vz_log,
                   height_log]).T
df.columns = ['time', 'x_distance', 'cas', 'thrust', 'lift', 'drag', 'alpha', 'teta','gamma', 'vz', 'height']
df.to_csv('afm_log.csv')

if uc.ms2kt(v) < v2min:
    print('Increase VR cannot reach V2min')
else:
    print(f'V2min reached with a VRmin {vr} kt')

cd_35ft = 2 * drag / (rho * S * v**2)
assert cd_35ft <= config_data['cd_max'], 'The a/c stalls, increase VR'

x_data = [config_data['cd'], config_data['cd_max']]
y_data = [config_data['alpha_min'], config_data['alpha_max']]
f_aoa = interp1d(x_data, y_data)
aoa_35ft = f_aoa(cd_35ft)

# TODO create a loop in order for the a/c speed to converge to JAR requirement(easy)
# Extra requirement from JAR 25.121(b)
gamma_jar_req = uc.perc2rad(2.4)
v_min_grad = (2 * (thrust * np.cos(np.radians(aoa_35ft)) - weight*np.sin(gamma_jar_req)) /
              (rho * S * cd_35ft)) ** .5

if v_min_grad > v:
    print('Increase VR. V2 does not reach the gradient requirement for ssg')
else:
    print(f'V2 reached sucessfully with a VR {vr}kt')










