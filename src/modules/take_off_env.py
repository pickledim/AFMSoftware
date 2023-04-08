import os
import json
import math
from pathlib import Path

import pandas as pd
import numpy as np

from scipy import constants
from scipy.interpolate import interp1d

from ambiance import Atmosphere

from src.tools import unit_conversion as uc


# from ground_phase import GroundRoll


class TakeOffPrep(object):
    """
    A class used to represent the preparation for a takeoff

    ...

    Attributes
    ----------
    mass : float
        the mass of the aircraft [kg]
    conf : str
        the configuration of the aircraft
    zp : float
        the altitude of the airport [ft]
    lg : str
        the landing gear
    eng_state : str
        the engine state
    timestep : float
        the time step

    Methods
    -------
    calculate_stall_speed()
        calculates the stall speed
    get_vmu()
        calculates the minimum unstick speed
    define_characteristic_speeds()
        defines the characteristic speeds
    calculate_v2_jar()
        calculates the takeoff safety speed

    """

    def __init__(
            self,
            input_variables: dict
    ):
        """
        Constructs necessary attributes for TakeOffPrep object.

        Parameters
        ----------
        mass : float
            the mass of the aircraft [kg]
        conf : str
            the configuration of the aircraft
        zp : float
            the altitude of the airport [ft]
        lg : str
            the landing gear
        eng_state : str
            the engine state
        timestep : float
            the time step

        """

        # performance data
        self.mass = input_variables["mass"]
        self.conf = input_variables["conf"]
        self.zp = input_variables["zp"]
        self.lg = input_variables["lg"]
        self.eng_state = input_variables["engine_state"]
        self.dt = input_variables["timestep"]

        # get config data
        self.data_dir = os.path.join(str(Path(__file__)).split("src")[0], "data")
        # print(self.data_dir)
        with open(os.path.join(self.data_dir, "config.json"), "r") as f:
            self.config_data = json.load(f)

        # find density
        self.rho = Atmosphere(uc.ft2m(self.zp)).density
        # get weight
        self.weight = self.mass * constants.g

        # get databases
        df_drag_polar = pd.read_excel(os.path.join(self.data_dir, "Aero_database.xlsx"), sheet_name="DragPolar")
        df_vmu = pd.read_excel(os.path.join(self.data_dir, "Database.xlsx"), sheet_name="Vmu")
        df_thrust = pd.read_excel(os.path.join(self.data_dir, "Database.xlsx"), sheet_name="Thrust")
        df_lift = pd.read_excel(os.path.join(self.data_dir, "Aero_database.xlsx"), sheet_name="Cz")

        # get ssg_thrust
        self.thrust_ssg = df_thrust[self.eng_state].values

        # get interpolation models
        self.drag_polar = interp1d(df_drag_polar["CZ2"].values, df_drag_polar["CX"].values)
        self.vmu_interp_model = interp1d(df_vmu["T/P"].values, df_vmu["kVs"].values)

        self.f_cl = interp1d(df_lift["ALPHA"].values, df_lift["CZ"].values)
        cl0 = self.f_cl(self.config_data['alpha_min'])
        cd0 = self.drag_polar(cl0**2)
        cd_min = self.drag_polar(df_drag_polar["CZ2"].min()**0.5)
        cd_max = self.drag_polar(df_drag_polar["CZ2"].max()**0.5)
        alpha_max = df_lift[df_lift["CZ"] <= (df_drag_polar["CZ2"].max()) ** 0.5]["ALPHA"].max()
        self.f_a_cd = interp1d(np.array([cd_min, cd_max]),
                               np.array([df_lift["ALPHA"].min(), df_lift["ALPHA"].max()]))

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
            "aoa_max": alpha_max,  # max angle of attack
            "clmax": self.config_data['cl_max'],  # approx
            "cl_max_op": df_lift["CZ"].max(),  # max cl
            "cl0": cl0,  # max cl
            "cd_max_op": cd_max,  # max cl
            "cd0": cd0,  # max cl
            "engine_coefs": self.config_data['engine_coefs'],
            # operational data
            "q": self.config_data['q'],  # rotation rate
            "teta_target": self.config_data['teta_target']  # theta target after lift off (theta law)
        }
        self.reached_airborne_phase = False
        self.vmca = self.config_data["vmca"]  # [kt]

        # define extra variables that will be used later on
        self.v_sta: float
        self.vmu: float
        self.v2min: float
        self.vr: float
        self.vef: float
        self.v_target: float
        self.characteristic_instants = dict()

    def calculate_stall_speed(self) -> None:
        """
        Calculates the stall speed of the aircraft in knots.

        Args:
            self (Aircraft): An instance of the Aircraft class.

        Returns:
            None: This function does not return anything. It updates the 'v_sta' attribute of the Aircraft instance with
             the calculated stall speed in knots.
        """
        weight, rho, s, clmax = self.constants_dict["weight"], self.constants_dict["rho"], \
                                self.constants_dict["S"], self.constants_dict["clmax"]

        v_stall = math.sqrt((2 * weight / (rho * s * clmax)))  # [m/s]
        self.v_sta = uc.ms2kt(v_stall)  # [kt]

    def v2_mu_limit(self) -> None:

        kVs = self.vmu_interp_model(self.thrust_ssg / self.weight)
        self.v2_mu_limit = kVs * self.v_sta  # [kt]

    def define_characteristic_speeds(self) -> None:
        """
        Calculates and sets the characteristic speeds of the aircraft using the previously calculated constants.

        This method calculates and sets the characteristic speeds of the aircraft, including V2min, VR, VEF, VMU, VMCA,
        and V_stall. The values of these speeds are calculated using the previously calculated constants for the weight,
        air density, wing area, and maximum lift coefficient.

        Returns:
            None
        """
        self.v2_mu_limit()
        v2min = max(1.13 * self.v_sta, 1.1 * self.vmca, self.v2_mu_limit)  # [kt], self.vmu
        vr = max(1.05 * self.v_sta, 1.05 * self.vmca)
        vef = vr - 2  # [kt]

        self.speeds = {
            "v2min": v2min,
            "vr": vr,
            "vef": vef,
            "vmca": self.vmca,
            "v_stall": self.v_sta
        }

    def calculate_v2_jar(self) -> None:
        """
        Calculates the target takeoff safety speed (V2) according to the JAR25.121(b) regulation.

        The JAR25.121(b) regulation requires that the minimum target takeoff safety speed (V2) must be at least
        1.13 times the stall speed (Vstall) or 1.1 times the minimum control speed in the air (Vmca), whichever is
        greater, and not less than the minimum unstick speed (Vmu).

        This method calculates V2 according to the JAR25.121(b) regulation by first determining the lift coefficient
        and drag coefficient at 400 ft above sea level. If the resulting climb gradient at this altitude is less than
        2.4%, then the V2 target must be adjusted to meet the regulation. Otherwise, V2 remains at the minimum target
        speed.

        Returns:
            None
        """

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
        """
        Initializes the variables and event logs for the simulation.

        Variables:
            - t: float, time [s]
            - x: float, distance traveled [m]
            - cas: float, calibrated airspeed [m/s]
            - tas: float, true airspeed [m/s]
            - cas_kt: float, calibrated airspeed [kt]
            - tas_kt: float, true airspeed [kt]
            - accel_cas: float, calibrated airspeed acceleration [m/s^2]
            - accel_tas: float, true airspeed acceleration [m/s^2]
            - dv: float, change in velocity [m/s]
            - dx: float, change in distance [m]

            - height: float, height above sea level [m]
            - v_z: float, vertical speed [m/s]

            - gamma: float, flight path angle [%]
            - gamma_rad: float, flight path angle [rad]
            - teta: float, pitch angle [deg]
            - aoa: float, angle of attack [deg]

            - sf_x: float, ground speed [m/s]
            - mass: float, aircraft mass [kg]
            - dt: float, time step size [s]

            - need_to_increase_Vr: boolean

        Event log:
            - t_log: list, log of time values [s]
            - x_log: list, log of distance values [m]
            - height_log: list, log of height values [m]
            - cas_kt_log: list, log of calibrated airspeed values [kt]
            - tas_kt_log: list, log of true airspeed values [kt]
            - vz_log: list, log of vertical speed values [m/s]
            - teta_log: list, log of pitch angle values [deg]
            - alpha_log: list, log of angle of attack values [deg]
            - gamma_log: list, log of flight path angle values [%]
            - thrust_log: list, log of thrust values [N]
            - lift_log: list, log of lift values [N]
            - drag_log: list, log of drag values [N]
        """

        self.variables = {
            "t": 0.,
            "x": 0.,
            "cas": 0.,
            "tas": 0.,
            "cas_kt": 0.,
            "tas_kt": 0.,
            "accel_cas": 0,
            "accel_tas": 0,
            "dv": 0,
            "dx": 0,

            "height": 0,
            "v_z": 0,

            "gamma": 0,
            "gamma_rad": 0,
            "teta": self.config_data['alpha_min'],
            "aoa": self.config_data['alpha_min'],

            "sf_x": 0,
            "mass": self.mass,
            "dt": self.dt,  # timestep

            "need_to_increase_Vr": False
        }

        self.event_log = {
            "t_log": [0.],
            "x_log": [0.],
            "height_log": [0.],
            "cas_kt_log": [0.],
            "tas_kt_log": [0.],
            "vz_log": [0],
            "teta_log": [self.variables["aoa"]],
            "alpha_log": [self.variables["aoa"]],
            "gamma_log": [0.],
            "thrust_log": [0.],
            "lift_log": [0.],
            "drag_log": [0.]
        }

    def update_values(self) -> None:
        """
        Update the values of the simulation parameters based on the current state of the aircraft.

        Modifies the following variables in self.variables dictionary:
        - cas
        - tas
        - x
        - t
        - cas_kt
        - tas_kt
        - rho

        Modifies the following variables in self.event_log dictionary:
        - cas_kt_log
        - tas_kt_log
        - t_log
        - x_log
        - gamma_log
        - vz_log
        - height_log
        - thrust_log
        - lift_log
        - drag_log
        - teta_log
        - alpha_log

        Returns:
        - None
        """

        self.variables["cas"] += self.variables["dv"]
        self.variables["tas"] += self.variables["dv"]
        self.variables["x"] += self.variables["dx"]
        self.variables["t"] += self.variables["dt"]

        self.variables["cas_kt"] = uc.ms2kt(self.variables["cas"])
        self.variables["tas_kt"] = uc.ms2kt(self.variables["tas"])

        self.variables["rho"] = Atmosphere(uc.ft2m(self.variables["height"])).density

        self.event_log["cas_kt_log"].append(float(self.variables["cas_kt"]))
        self.event_log["tas_kt_log"].append(float(self.variables["tas_kt"]))
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

    def check_if_v2_reached(self) -> None:
        """
        Checks if the a/c reached the v2 min.

        :return: None
        """

        if self.variables["need_to_increase_Vr"]:
            if self.reached_airborne_phase:
                print("Could not reach V2")
            print("Increase VR the a/c cannot reach the target :(\n")

    def pilot_preparation(self) -> None:
        """
        Prepares the pilot by calling several methods:
        - calculate_stall_speed: calculates the stall speed of the aircraft
        - get_vmu: gets the minimum unstick speed
        - define_characteristic_speeds: defines the characteristic speeds of the aircraft
        - calculate_v2_jar: calculates V2 for the given altitude and weight
        - initialize_data: initializes the data dictionary and event log for the simulation
        """

        self.calculate_stall_speed()
        self.define_characteristic_speeds()
        self.calculate_v2_jar()
        self.initialize_data()
