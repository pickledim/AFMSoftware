from typing import Dict, Tuple
import numpy as np

from ambiance import Atmosphere

import simple_pid

import unit_conversion as uc

from take_off_env import TakeOffPrep


class AirbornePhase(TakeOffPrep):
    """
    Simulates the airborne phase of the aircraft's flight.
    """
    def pid_control(self) -> None:
        """
        Use a PID controller to adjust the stick input based on the difference between the target and
        current acceleration.
        """

        # Use the PID controller to adjust the throttle based on the difference between the target and
        # current angle of attack
        diff = self.target_acceleration - self.variables["accel"]
        if diff > 0:
            print()
        stick_adjustment = diff #self.pid_acceleration(diff)

        self.variables["aoa"] -= stick_adjustment *1e-1
        print(self.variables["aoa"])

    def normal_climb(self) -> None:
        """
        Calculate the aircraft's motion during a steady climb, using equations of motion for a constant climb
        and a drag polar. Update the aircraft's state.
        """

        tas = self.variables["tas"]
        thrust, weight = self.variables["thrust"], self.constants_dict["weight"]
        rho, S, mass = self.rho, self.constants_dict["S"], self.variables["mass"]

        # In order to reach constant climb -> W=L -> find the corresponding cl and cd
        cl = 2 * weight / (rho * S * tas ** 2)
        self.variables["lift"] = 0.5 * (rho * S * cl * tas ** 2)

        cd = self.drag_polar(cl ** 2)
        self.variables["drag"] = 0.5 * (rho * S * cd * tas ** 2)

        self.variables["aoa"] = self.f_a_cd(cd)

        # other method not optimal
        # self.variables["gamma_rad"] = (self.variables["thrust"] * np.cos(uc.deg2rad(self.variables["aoa"])) /
        #                                self.constants_dict["weight"]) - (cd / cl)

        # aafp p299
        t_aero = thrust * np.cos(np.deg2rad(self.variables["aoa"]))
        self.variables["gamma_rad"] = (t_aero - self.variables["drag"]) / self.constants_dict["weight"]

        # aafp p306
        # dgamma = (thrust * np.sin(np.radians(self.variables["aoa"])) + self.variables["lift"] - weight) / \
        #          (mass * tas) * self.variables["dt"]

        # self.variables["gamma_rad"] += dgamma
        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])

    def calculate_angles(self) -> None:
        """
        Calculate the aircraft's angle of attack, flight path angle, and other quantities related to the
        aircraft's motion.
        """

        aoa, gamma = self.variables['aoa'], self.variables['gamma']
        thrust, v = self.variables["thrust"], self.variables["cas"]
        rho, S = self.rho, self.constants_dict["S"]
        weight, mass = self.constants_dict['weight'], self.constants_dict['mass']

        tas = self.variables["tas"] = (rho / self.rho0) ** 0.5 * v

        T_aero = thrust * np.sin(np.radians(aoa))
        cos_gama = np.cos(np.radians(self.variables['gamma']))
        drag = 0.5 * rho * S * self.f_cd(aoa)  * tas ** 2
        lift = 0.5 * rho * S * self.f_cl(aoa)  * tas ** 2
        #
        # accel = (self.variables["tas"] - uc.kt2ms(self.event_log["tas_kt_log"][-1])) / self.variables["dt"]
        #
        # self.variables["gamma_rad"] = (T_aero - drag - mass * (self.variables["dv"] / self.variables["dt"])) / weight

        # advanced ac perfo p306
        # * np.sin(np.radians(self.variables["aoa"]))
        dgamma = (T_aero + lift - weight * cos_gama) / \
                 (self.variables["mass"] * self.variables["tas"]) * self.variables["dt"]

        self.variables["gamma_rad"] += dgamma

        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])

        self.variables["teta"] = self.variables["aoa"] + self.variables["gamma"]

    def calculate_equations_of_motion(self) -> Tuple[Dict[str, float], float]:
        """
        Calculate the forces and acceleration acting on the aircraft using equations of motion for an aircraft in
        accelerating climb.

        Return a dictionary containing the drag, lift, and net force acting on the aircraft,
        as well as the aircraft's acceleration.
        """
        aoa, gamma = self.variables['aoa'], self.variables['gamma']
        thrust, v = self.variables["thrust"], self.variables["cas"]
        rho, S = self.rho, self.constants_dict["S"]
        weight, mass = self.constants_dict['weight'], self.constants_dict['mass']
        tas = self.variables["tas"]  # = (rho/self.rho0)**0.5 * v


        T_aero = thrust * np.cos(np.radians(aoa))

        drag = 0.5 * rho * S * self.f_cd(aoa) * tas ** 2

        lift = 0.5 * rho * S * self.f_cl(aoa) * tas ** 2

        SFx = T_aero - drag - weight * np.sin(np.radians(gamma))

        accel = SFx / mass
        print(accel)

        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_x": SFx,
        }

        return forces, accel

    def update_aerial_values(self) -> None:
        """
        Update the aircraft's velocities, change in altitude, and change in horizontal position.
        """
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
        """
        Simulates the airborne phase of the aircraft's flight.

        During this phase, the aircraft climbs from takeoff to a height of 35 ft. The climb is divided into two parts:
        an accelerating climb until the minimum takeoff safety speed V2 is reached, and a steady climb after that.
        The function updates the aircraft's state variables at each iteration and saves the and characteristic instants.

        Returns:
            None.
        """

        # super().transition_phase()
        # TODO: Study the climb part p299 the Vz is low
        self.reached_v2 = False

        self.rho0 = Atmosphere(uc.ft2m(self.zp)).density

        # Define the target acceleration
        self.target_acceleration = 0.0  # [m/s^2]

        # Define the PID controller for acceleration
        # self.pid_acceleration = simple_pid.PID(0.3, 0.05, 0.02, setpoint=self.target_acceleration)

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

            # CAS steady climb
            # else:
            #     self.normal_climb()
            #
            # if self.variables["cas_kt"] >= self.speeds["v_target"]:
            #     if not self.reached_v2:
            #         print(f'V2min reached @ {round(float(self.variables["height"]), 2)}ft!')
            #         self.characteristic_instants["v2"] = {"Instant": self.variables["t"],
            #                                               "Speed": self.variables["cas_kt"]}
            #         self.reached_v2 = True

            self.update_aerial_values()
            super().update_values()

        self.characteristic_instants["v35ft"] = {"Instant": self.variables["t"],
                                                 "Speed": self.variables["cas_kt"]}
