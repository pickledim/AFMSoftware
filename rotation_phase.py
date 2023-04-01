from typing import Dict, Tuple

import numpy as np

from ground_phase import GroundRoll


class RotationPhase(GroundRoll):
    def calculate_forces_in_rotation(self, aoa: float) -> Tuple[Dict[str, float], float]:
        f_cl = self.f_cl
        f_cd = self.f_cd

        rho = self.constants_dict["rho"]
        S = self.constants_dict["S"]
        mu = self.constants_dict["mu"]
        weight = self.constants_dict["weight"]
        mass = self.constants_dict["mass"]
        v = self.variables["v"]
        thrust = self.variables["thrust"]

        drag = 0.5 * rho * S * f_cd(aoa) * v ** 2
        lift = 0.5 * rho * S * f_cl(aoa) * v ** 2
        f_roll = abs(mu * (lift - weight))

        f_x = thrust * np.cos(np.radians(aoa)) - drag - f_roll

        accel = f_x / mass

        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_roll": f_roll,
            "F_x": f_x,
        }

        return forces, accel

    def transition_phase(self):
        # super().up_to_rotation()

        self.variables["need_to_increase_Vr"] = False

        aoa0 = self.variables["aoa"]
        dt = self.variables["dt"]
        weight = self.constants_dict["weight"]
        q = self.constants_dict["q"]
        alpha_max = self.constants_dict["aoa_max"]

        while round(float(self.variables["lift"]), 1) < round(weight, 1):
            aoa = q * dt + aoa0

            if aoa > alpha_max:
                self.variables["need_to_increase_Vr"] = True
                print("\nIncrease VR the a/c cannot Lift Off :(\n")
                return

            forces, accel = self.calculate_forces_in_rotation(aoa)

            dv = accel * dt
            dx = 0.5 * accel * dt ** 2 + self.variables["v"] * dt

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
                                                   "Speed": self.variables["v_kt"]}


if __name__ == "__main__":
    mass = 80000.0  # [kg]
    conf = "1+F"
    zp = 0.0  # [ft]
    lg = "Up"
    engine_state = "OEI"

    a320_to = RotationPhase(mass, conf, zp, lg, engine_state)
    a320_to.transition_phase()
