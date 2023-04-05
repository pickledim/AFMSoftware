import numpy as np
import pandas as pd

from ambiance import Atmosphere

import simple_pid

import unit_conversion as uc

# from rotation_phase import RotationPhase
from take_off_env import TakeOffPrep


class AirbornePhase(TakeOffPrep):

    def pid_control(self) -> None:

        # Use the PID controller to adjust the throttle based on the difference between the target and
        # current angle of attack
        throttle_adjustment = self.pid_acceleration(self.target_acceleration - self.variables["accel_tas"])

        self.variables["aoa"] += throttle_adjustment * 1e-2

    def normal_climb(self) -> None:

        tas = self.variables["tas"]
        thrust, weight = self.variables["thrust"], self.constants_dict["weight"]
        rho, S, mass = self.rho, self.constants_dict["S"], self.variables["mass"]

        # In order to reach constant climb -> W=L -> find the corresponding cl and cd
        cl = 2 * weight / (rho * S * tas ** 2)
        self.variables["lift"] = 0.5 * (self.rho * S * cl * tas ** 2)

        cd = self.drag_polar(cl ** 2)
        self.variables["drag"] = 0.5 * (self.rho * S * cd * tas ** 2)

        self.variables["aoa"] = self.f_a_cd(cd)

        # self.variables["gamma_rad"] = (self.variables["thrust"] * np.cos(uc.deg2rad(self.variables["aoa"])) /
        #                                self.constants_dict["weight"]) - (cd / cl)

        # aafp p299
        t_aero = thrust * np.cos(np.deg2rad(self.variables["aoa"]))
        self.variables["gamma_rad"] = (t_aero - self.variables["drag"]) / self.constants_dict["weight"]

        # aafp p306
        dgamma = (thrust * np.sin(np.radians(self.variables["aoa"])) + self.variables["lift"] - weight) / \
                 (mass * tas) * self.variables["dt"]

        self.variables["gamma_rad"] += dgamma
        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])


        # update

        # dh = self.variables["tas"] * np.sin(self.variables["gamma_rad"]) * self.variables["dt"]
        # self.variables["dv"] = 0
        # self.variables["dx"] = self.variables["tas"] * np.cos(self.variables["gamma_rad"]) * self.variables["dt"]
        # self.variables["v_z"] = self.variables["tas"] * np.sin(self.variables["gamma_rad"])
        # self.variables["height"] += uc.m2ft(dh)
        # self.rho = Atmosphere(uc.ft2m(self.variables["height"])).density
        # self.variables["tas"] = (self.rho0 / self.rho) ** 0.5 * self.variables["cas"]

            # super().update_values()

    def calculate_angles(self):
        aoa, gamma = self.variables['aoa'], self.variables['gamma']
        thrust, v = self.variables["thrust"], self.variables["cas"]
        rho, S = self.rho, self.constants_dict["S"]
        weight, mass = self.constants_dict['weight'], self.constants_dict['mass']

        tas = self.variables["tas"] = (rho / self.rho0) ** 0.5 * v

        T_aero = thrust * np.cos(np.radians(aoa))

        drag = 0.5 * rho * S * self.f_cd(aoa) * tas ** 2
        #
        # accel = (self.variables["tas"] - uc.kt2ms(self.event_log["tas_kt_log"][-1])) / self.variables["dt"]
        #
        self.variables["gamma_rad"] = (T_aero - drag - mass * (self.variables["dv"] / self.variables["dt"])) / weight
        # advanced ac perfo p306
        dgamma = (thrust * np.sin(np.radians(self.variables["aoa"])) + self.variables["lift"] - weight) / \
                 (self.variables["mass"] * self.variables["tas"]) * self.variables["dt"]

        self.variables["gamma_rad"] += dgamma

        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])

        self.variables["teta"] = self.variables["aoa"] + self.variables["gamma"]

    def calculate_equations_of_motion(self) -> None:

        aoa, gamma = self.variables['aoa'], self.variables['gamma']
        thrust, v = self.variables["thrust"], self.variables["cas"]
        rho, S = self.rho, self.constants_dict["S"]
        weight, mass = self.constants_dict['weight'], self.constants_dict['mass']

        tas = self.variables["tas"] = (rho/self.rho0)**0.5 * v

        T_aero = thrust * np.cos(np.radians(aoa))

        drag = 0.5 * rho * S * self.f_cd(aoa) * tas ** 2

        lift = 0.5 * rho * S * self.f_cl(aoa) * tas ** 2

        SFx = T_aero - drag - weight * np.sin(np.radians(gamma))

        accel = SFx / mass

        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_x": SFx,
        }

        return forces, accel

    def update_aerial_values(self) -> None:

        # self.variables["v_z"] = (self.variables["thrust"] - self.variables["drag"]) / weight * self.variables["v"] - \
        #                        (self.variables["v"]/9.81 * self.variables["accel"])
        self.variables["v_z"] = self.variables["tas"] * np.sin(self.variables["gamma_rad"])

        if self.reached_v2:
            self.variables["dv"] = 0
        else:
            self.variables["dv"] = self.variables["accel"] * self.variables["dt"]

        self.variables["dx"] = self.variables["tas"] * np.cos(self.variables["gamma_rad"]) * self.variables["dt"]
        # dh * np.tan(gamma_rad)

        dh = self.variables["tas"] * np.sin(self.variables["gamma_rad"]) * self.variables["dt"]  # a_z*dt**2 + v_z * dt

        # update
        self.variables["height"] += uc.m2ft(dh)
        self.rho = Atmosphere(uc.ft2m(self.variables["height"])).density
        self.variables["tas"] = (self.rho0 / self.rho) ** 0.5 * self.variables["cas"]

    def airborne_phase(self) -> None:

        # super().transition_phase()
        # TODO: Study the climb part p299 the Vz is low
        self.reached_v2 = False

        self.rho0 = Atmosphere(uc.ft2m(self.zp)).density

        # Define the target acceleration
        self.target_acceleration = 0.0  # [m/s^2]

        # Define the PID controller for acceleration
        self.pid_acceleration = simple_pid.PID(1.0, 0.1, 0.05, setpoint=self.target_acceleration)

        while self.variables["height"] < 35.:

            # CAS accelerating climb
            if not self.reached_v2:
                # control alpha in order to get 0 accel
                self.pid_control()

                # get the new angles
                self.calculate_angles()

                # get the forces
                forces, accel = self.calculate_equations_of_motion()

                self.variables.update(
                    {
                        "drag": forces["Drag"],
                        "lift": forces["Lift"],
                        "sf_x": forces["F_x"],
                        "accel": accel
                    }
                )
                # self.update_aerial_values()
            # CAS steady climb
            else:
                self.normal_climb()

            if self.variables["cas_kt"] >= self.speeds["v_target"]:
                if not self.reached_v2:
                    print(f'V2min reached @ {round(float(self.variables["height"]), 2)}ft!')
                    self.characteristic_instants["v2"] = {"Instant": self.variables["t"],
                                                          "Speed": self.variables["cas_kt"]}
                    self.reached_v2 = True

            self.update_aerial_values()
            super().update_values()

        self.characteristic_instants["v35ft"] = {"Instant": self.variables["t"],
                                                 "Speed": self.variables["cas_kt"]}
