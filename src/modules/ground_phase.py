from typing import Tuple

import numpy as np


from src.modules.take_off_env import TakeOffPrep


def get_thrust(speed: float, coefs: list) -> np.array:
    """
    Calculate the thrust produced by the aircraft's engine based on the speed and a list of coefficients.

    Args:
        speed (float): The speed of the aircraft in knots.
        coefs (list): A list of 8 coefficients representing the engine thrust as a polynomial function of speed.

    Returns:
        float: The thrust produced by the engine in Newtons.
    """

    speed_vector = np.array([1, speed, speed ** 2, speed ** 3, speed ** 4, speed ** 5, speed ** 6, speed ** 7],
                            dtype='object')
    _coefs = np.array(coefs)

    _thrust = np.dot(speed_vector, _coefs)

    return _thrust


class GroundRoll(TakeOffPrep):
    """
    A class representing the ground roll phase of an aircraft takeoff.

    Attributes:
    constants_dict (dict): A dictionary of constants used in the calculations.
    speeds (dict): A dictionary of speeds used in the calculations.
    variables (dict): A dictionary of variables used in the calculations.
    characteristic_instants (dict): A dictionary of characteristic instants in the takeoff process.

    """
    def get_thrust(self) -> float:
        """
        Get the thrust produced by the engine.

        Returns:
        float: The thrust produced by the engine.

        """

        engine_coefs = self.constants_dict["engine_coefs"]

        # check if vef has been reached and adjust thrust accordingly
        if self.vef_inst:
            if self.variables["cas_kt"] < self.speeds["vef"]:
                self.variables["thrust"] = get_thrust(self.variables["cas"], engine_coefs)
            else:
                self.variables["thrust"] = self.variables["thrust"]/2
                self.vef_inst = False

        return self.variables["thrust"]

    def calculate_forces(self,) -> Tuple[dict, float]:
        """
        Calculate the forces acting on the aircraft during the ground roll.

        Returns:
        Tuple[dict, float]: A tuple containing a dictionary of forces and the acceleration of the aircraft.

        """

        rho, S, cd0, cl0, mu = self.constants_dict["rho"], self.constants_dict["S"], \
                               self.constants_dict["cd0"], self.constants_dict["cl0"], \
                               self.constants_dict["mu"]

        weight, mass = self.constants_dict["weight"], self.constants_dict["mass"]

        # calculate thrust, drag, lift, and roll forces
        thrust = self.get_thrust()
        lift = 0.5 * rho * S * cl0 * self.variables["cas"] ** 2
        drag = 0.5 * rho * S * self.drag_polar(cl0**2) * self.variables["cas"] ** 2

        f_roll = abs(mu * (weight - lift))

        # calculate net force and acceleration
        sf_x = thrust - drag - f_roll
        accel = sf_x / mass

        forces = {
            "Thrust": thrust,
            "Drag": drag,
            "Lift": lift,
            "F_roll": f_roll,
            "F_x": sf_x
        }

        return forces, accel

    def up_to_rotation(self) -> None:
        """
        Perform calculations for the ground roll up to the rotation speed.
        """
        # super().pilot_preparation()

        dt, vr = self.variables["dt"], self.speeds["vr"]

        self.vef_inst = True
        # loop until rotation speed is reached
        while self.variables["cas_kt"] < vr:

            # calculate forces
            forces, accel = self.calculate_forces()

            # calculate change in velocity and displacement
            dv = accel * dt
            dx = 0.5 * accel * dt ** 2 + self.variables["cas"] * dt

            # update variables
            self.variables["drag"], self.variables["lift"], self.variables["sf_x"] = \
                forces["Drag"], forces["Lift"], forces["F_x"]

            self.variables["accel"], self.variables["dv"], self.variables["dx"] = accel, dv, dx,

            super().update_values()

        self.characteristic_instants["Rotation"] = {"Instant": self.variables["t"],
                                                    "Speed": self.variables["cas_kt"]}
