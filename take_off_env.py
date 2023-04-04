import json
import math

import pandas as pd
import numpy as np

from scipy import constants
from scipy.interpolate import interp1d

from ambiance import Atmosphere

import unit_conversion as uc


# from ground_phase import GroundRoll


class TakeOffPrep(object):
    def __init__(
            self,
            mass: float,
            conf: str,
            zp: float,
            lg: str,
            eng_state: str,
            timestep: float
    ):

        # performance data
        self.mass = mass
        self.conf = conf
        self.zp = zp
        self.lg = lg
        self.eng_state = eng_state
        self.dt = timestep

        # get config data
        with open('config.json', 'r') as f:
            self.config_data = json.load(f)

        # find density
        self.rho = Atmosphere(uc.ft2m(self.zp)).density
        # get weight
        self.weight = self.mass * constants.g

        self.constants_dict = {
            "rho": self.rho,
            "S": self.config_data['S'],
            "mu": self.config_data['mu'],
            "weight": self.weight,
            "mass": self.mass,
            "conf": self.conf,
            "zp": self.zp,
            "lg": self.lg,
            "eng_state": self.eng_state,
            # aero data
            "aoa0": self.config_data['alpha_min'],  # initial angle of attack
            "aoa_max": self.config_data['alpha_max'],  # max angle of attack
            "clmax": 2.4,
            "cl_max_op": self.config_data['cl_max'],  # max cl
            "cl0": self.config_data['cl'],  # max cl
            "cd_max_op": self.config_data['cd_max'],  # max cl
            "cd0": self.config_data['cd'],  # max cl
            "engine_coefs": self.config_data['engine_coefs'],  # max cl
            # operational data
            "q": self.config_data['q'],  # rotation rate
            "teta_target": self.config_data['teta_target']  # theta target after lift off (theta law)
        }

        self.vmca = self.config_data["vmca"]  # [kt]

        # get databases
        df_drag_polar = pd.read_excel("Database.xlsx", sheet_name="DragPolar")
        df_vmu = pd.read_excel("Database.xlsx", sheet_name="Vmu")
        df_thrust = pd.read_excel("Database.xlsx", sheet_name="Thrust")

        # get ssg_thrust
        self.thrust_ssg = df_thrust[self.eng_state].values

        # get interpolation models
        self.drag_polar = interp1d(df_drag_polar["Cz2"].values, df_drag_polar["Cx"].values)
        self.vmu_interp_model = interp1d(df_vmu["T/P"].values, df_vmu["kVs"].values)

        self.f_cl = interp1d(np.array([self.constants_dict["aoa0"], self.constants_dict["aoa_max"]]),
                             np.array([self.constants_dict["cl0"], self.constants_dict["cl_max_op"]]))

        self.f_cd = interp1d(np.array([self.constants_dict["aoa0"], self.constants_dict["aoa_max"]]),
                             np.array([self.constants_dict["cd0"], self.constants_dict["cd_max_op"]]))

        self.f_a_cd = interp1d(np.array([self.constants_dict["cd0"], self.constants_dict["cd_max_op"]]),
                               np.array([self.constants_dict["aoa0"], self.constants_dict["aoa_max"]]))

        # define extra variables that will be used later on
        self.v_sta: float
        self.vmu: float
        self.v2min: float
        self.vr: float
        self.vef: float
        self.v_target: float
        self.characteristic_instants = dict()

    def calculate_stall_speed(self) -> None:

        weight, rho, s, clmax = self.constants_dict["weight"], self.constants_dict["rho"], self.constants_dict["S"],\
                                self.constants_dict["clmax"]

        v_stall = math.sqrt((2 * weight / (rho * s * clmax)))  # [m/s]
        self.v_sta = uc.ms2kt(v_stall)  # [kt]

    def get_vmu(self) -> None:

        kVs = self.vmu_interp_model(self.thrust_ssg / self.weight)
        self.vmu = kVs * self.v_sta  # [kt]

    def define_characteristic_speeds(self) -> None:

        v2min = max(1.13*self.v_sta, 1.1*self.vmca, self.vmu)  # [kt]
        vr = max(1.05*self.v_sta, 1.05*self.vmca)  # [kt] # [kt]
        vef = vr - 2  # [kt]

        self.speeds = {
            "v2min": v2min,
            "vr": vr,
            "vef": vef,
            "vmu": self.vmu,
            "vmca": self.vmca,
            "v_stall": self.v_sta
        }

    def calculate_v2_jar(self) -> None:

        # constants
        weight, s, v2min = self.constants_dict["weight"], self.constants_dict["S"], self.speeds["v2min"]

        # zp at ssg level (400ft)
        zp_ssg = uc.ft2m(400)
        rho_ssg = Atmosphere(zp_ssg).density

        cl_ssg = 2 * weight / (rho_ssg * s * uc.kt2ms(v2min) ** 2)
        cd_ssg = self.drag_polar(cl_ssg ** 2)

        gamma_ssg = 100 * ((self.thrust_ssg / weight) - (cd_ssg / cl_ssg))

        if gamma_ssg < 2.4:
            gamma_jar_req = uc.perc2rad(2.4)
            v_target_ms = math.sqrt(2 * (self.thrust_ssg - weight * gamma_jar_req) / (rho_ssg * s * cd_ssg))
            v_target = uc.ms2kt(v_target_ms)
            print(f"\nJAR Assertion: The given v2min does not pass the JAR25.121(b) Regulation thus the V2_target has "
                  f"to change. New V2 target {v_target} kt\n")

        else:
            v_target = v2min
            print(f"\nJAR Assertion: SSG of V2min {round(gamma_ssg[0], 1)}% Passed\n")

        self.speeds["v_target"] = v_target

    def initialize_data(self) -> None:

        self.variables = {
            "t": 0.,
            "x": 0.,
            "v": 0.,
            "v_kt": 0.,

            "accel": 0,
            "dv": 0,
            "dx": 0,

            "height": 0,
            "v_z": 0,

            "gamma": 0,
            "gamma_rad": 0,
            "teta": 0,
            "aoa": self.config_data['alpha_min'],

            "sf_x": 0,
            "mass": self.mass,
            "dt": self.dt  # timestep
        }

        self.event_log = {
            "t_log": [0],
            "x_log": [0],
            "height_log": [0],
            "v_kt_log": [0],
            "vz_log": [0],
            "teta_log": [0],
            "alpha_log": [0],
            "gamma_log": [0],
            "thrust_log": [0],
            "lift_log": [0],
            "drag_log": [0]
        }

    def update_values(self) -> None:

        self.variables["v"] += self.variables["dv"]
        self.variables["x"] += self.variables["dx"]
        self.variables["t"] += self.variables["dt"]

        self.variables["v_kt"] = uc.ms2kt(self.variables["v"])

        self.variables["rho"] = Atmosphere(uc.ft2m(self.variables["height"])).density

        self.event_log["v_kt_log"].append(float(self.variables["v_kt"]))
        self.event_log["t_log"].append(float(self.variables["t"]))
        self.event_log["x_log"].append(float(self.variables["x"]))
        self.event_log["gamma_log"].append(float(self.variables["gamma"]))
        self.event_log["vz_log"].append(float(self.variables["v_z"]))
        self.event_log["height_log"].append(float(self.variables["height"]))
        self.event_log["thrust_log"].append(float(self.variables["thrust"]))
        self.event_log["lift_log"].append(float(self.variables["lift"]))
        self.event_log["drag_log"].append(float(self.variables["drag"]))
        self.event_log["teta_log"].append(float(self.variables["teta"]))
        self.event_log["alpha_log"].append(float(self.variables["aoa"]))

    def pilot_preparation(self) -> None:

        self.calculate_stall_speed()
        self.get_vmu()
        self.define_characteristic_speeds()
        self.calculate_v2_jar()
        self.initialize_data()
