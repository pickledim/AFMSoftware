from typing import Tuple

import numpy as np


from take_off_env import TakeOffPrep


def get_thrust(speed: float, coefs: list) -> float:
    """

    :param speed:
    :param coefs:
    :return:
    """

    speed_vector = np.array([1, speed, speed ** 2, speed ** 3, speed ** 4, speed ** 5, speed ** 6, speed ** 7],
                            dtype='object')
    _coefs = np.array(coefs)

    _thrust = np.dot(speed_vector, _coefs)

    return _thrust


class GroundRoll(TakeOffPrep):

    def get_thrust(self) -> float:

        engine_coefs = self.constants_dict["engine_coefs"]

        # check if vef has been reached and adjust thrust accordingly
        if self.vef_inst:
            if self.variables["v_kt"] < self.speeds["vef"]:
                self.variables["thrust"] = get_thrust(self.variables["v"], engine_coefs)
            else:
                self.variables["thrust"] = self.variables["thrust"]/2
                self.vef_inst = False

        return self.variables["thrust"]

    def calculate_forces(self,) -> Tuple[dict, float]:

        rho, S, cd0, cl0, mu = self.constants_dict["rho"], self.constants_dict["S"], \
                               self.constants_dict["cd0"], self.constants_dict["cl0"], \
                               self.constants_dict["mu"]

        weight, mass = self.constants_dict["weight"], self.constants_dict["mass"]

        # calculate thrust, drag, lift, and roll forces
        thrust = self.get_thrust()
        drag = 0.5 * rho * S * cd0 * self.variables["v"] ** 2
        lift = 0.5 * rho * S * cl0 * self.variables["v"] ** 2
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

    def up_to_rotation(self):
        """
        Perform calculations for the ground roll up to the rotation speed.
        """
        TakeOffPrep.pilot_preparation(self)

        dt, vr = self.variables["dt"], self.speeds["vr"]

        self.vef_inst = True
        # loop until rotation speed is reached
        while self.variables["v_kt"] < vr:

            # calculate forces
            forces, accel = self.calculate_forces()

            # calculate change in velocity and displacement
            dv = accel * dt
            dx = 0.5 * accel * dt ** 2 + self.variables["v"] * dt

            # update variables
            self.variables["drag"], self.variables["lift"], self.variables["sf_x"] = forces["Drag"], \
                                                                                     forces["Lift"], forces["F_x"]
            self.variables["accel"], self.variables["dv"], self.variables["dx"] = accel, dv, dx,

            TakeOffPrep.update_values(self)


        self.characteristic_instants["Rotation"] = self.variables["t"]


if __name__ == "__main__":

    mass = 60000.0  # [kg]
    conf = '1+F'
    zp = 0.  # [ft]
    lg = 'Up'
    engine_state = 'OEI'

    a320_to = GroundRoll(mass, conf, zp, lg, engine_state)
    a320_to.up_to_rotation()