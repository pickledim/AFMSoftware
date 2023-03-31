from typing import Tuple
import pandas as pd
import numpy as np

from ground_phase import GroundRoll
from take_off_env import TakeOffPrep


class RotationPhase(GroundRoll):

    def calculate_forces_in_rotation(self,  aoa: float) -> Tuple[dict, float]:

        # interpolation models
        f_cl = self.f_cl
        f_cd = self.f_cd

        rho, S, cd0, cl0, mu = self.constants_dict["rho"], self.constants_dict["S"], \
                               self.constants_dict["cd0"], self.constants_dict["cl0"], \
                               self.constants_dict["mu"]

        weight, mass, v = self.constants_dict["weight"], self.constants_dict["mass"],  self.variables["v"]

        drag = 0.5 * rho * S * f_cd(aoa) * v ** 2
        lift = 0.5 * rho * S * f_cl(aoa) * v ** 2
        f_roll = abs(mu * (lift - weight))

        f_x = self.variables["thrust"] * np.cos(np.radians(aoa)) - drag - f_roll  # aero axis

        accel = f_x / mass

        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_roll": f_roll,
            "F_x": f_x
        }

        return forces, accel

    def transition_phase(self):

        GroundRoll.up_to_rotation(self)

        self.variables["limit_aoa"] = False

        aoa0 = self.variables["aoa"]
        dt, weight = self.variables["dt"], self.constants_dict["weight"]
        q, alpha_max = self.constants_dict["q"], self.constants_dict["aoa_max"]

        # criterion: W=L
        while round(float(self.variables["lift"]), 1) < round(weight, 1):

            aoa = q * dt + aoa0

            if aoa > alpha_max:
                self.variables["limit_aoa"] = True
                print('\nIncrease VR the a/c cannot Lift Off :(\n')
                return

            forces, accel = self.calculate_forces_in_rotation(aoa)

            dv = accel * dt
            dx = 0.5 * accel * dt ** 2 + self.variables["v"] * dt

            # update variables
            self.variables["drag"], self.variables["lift"], self.variables["sf_x"] = forces["Drag"], \
                                                                                     forces["Lift"], forces["F_x"]
            self.variables["accel"], self.variables["dv"], self.variables["dx"], self.variables["aoa"] = \
                accel, dv, dx, aoa

            GroundRoll.update_values(self)
            aoa0 = aoa

        self.characteristic_instants["LiftOff"] = self.variables["t"]

if __name__ == "__main__":

    mass = 60000.0  # [kg]
    conf = '1+F'
    zp = 0.  # [ft]
    lg = 'Up'
    engine_state = 'OEI'

    a320_to = RotationPhase(mass, conf, zp, lg, engine_state)
    a320_to.transition_phase()