from typing import Dict, Tuple

import numpy as np


from take_off_env import TakeOffPrep


class RotationPhase(TakeOffPrep):
    """
    Performs the rotation phase of takeoff, increasing the aircraft's pitch angle until lift off.
    """
    def calculate_forces_in_rotation(self, aoa: float) -> Tuple[Dict[str, float], float]:
        """
        Calculates the forces acting on the aircraft during rotation phase for a given angle of attack.

        Args:
            aoa (float): Angle of attack in degrees.

        Returns:
            Tuple[Dict[str, float], float]: Tuple containing a dictionary of forces and the resulting acceleration.
        """
        # redefine the variables from the class for ease of reading
        f_cl = self.f_cl
        f_cd = self.f_cd

        rho = self.constants_dict["rho"]
        S = self.constants_dict["S"]
        mu = self.constants_dict["mu"]
        weight = self.constants_dict["weight"]
        mass = self.constants_dict["mass"]
        v = self.variables["cas"]
        thrust = self.variables["thrust"]

        # compute the acting forces
        drag = 0.5 * rho * S * f_cd(aoa) * v ** 2
        lift = 0.5 * rho * S * f_cl(aoa) * v ** 2
        f_roll = abs(mu * (lift - weight))

        # compute the resulting force
        f_x = thrust * np.cos(np.radians(aoa)) - drag - f_roll

        # compute the acceleration
        accel = f_x / mass

        # stack the forces in a dictionary
        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_roll": f_roll,
            "F_x": f_x,
        }

        return forces, accel

    def transition_phase(self) -> None:
        """
        Performs the rotation phase of takeoff, increasing the aircraft's pitch angle until lift off.

        Returns:
            None
        """

        # redefine the variables from the class for ease of reading
        aoa0 = self.variables["aoa"]
        dt = self.variables["dt"]
        weight = self.constants_dict["weight"]
        q = self.constants_dict["q"]
        alpha_max = self.constants_dict["aoa_max"]

        while round(float(self.variables["lift"]), 1) < round(weight, 1):
            aoa = q * dt + aoa0

            if aoa > alpha_max:
                self.variables["need_to_increase_Vr"] = True
                self.check_if_v2_reached()
                return

            forces, accel = self.calculate_forces_in_rotation(aoa)

            dv = accel * dt
            dx = 0.5 * accel * dt ** 2 + self.variables["cas"] * dt

            self.variables.update(
                {
                    "drag": forces["Drag"],
                    "lift": forces["Lift"],
                    "sf_x": forces["F_x"],
                    "accel": accel,
                    "dv": dv,
                    "dx": dx,
                    "aoa": aoa,
                    "teta": aoa
                }
            )

            super().update_values()
            aoa0 = aoa

        self.characteristic_instants["LiftOff"] = {"Instant": self.variables["t"],
                                                   "Speed": self.variables["cas_kt"]}
