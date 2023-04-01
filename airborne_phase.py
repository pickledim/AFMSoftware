import numpy as np

import simple_pid

import unit_conversion as uc

from rotation_phase import RotationPhase


class AirbornePhase(RotationPhase):

    def control_accel(self):

        accel = self.variables["accel"]

        # Define the target acceleration
        target_acceleration = 0.0  # [m/s^2]

        # Define the PID controller for acceleration
        pid_acceleration = simple_pid.PID(1.0, 0.1, 0.05, setpoint=target_acceleration)

        # Main loop
        # while True:
        # Read the current acceleration from the aircraft's sensors
        current_acceleration = accel

        # Calculate the current angle of attack based on the current acceleration
        # current_angle_of_attack = aoa

        # Use the PID controller to adjust the throttle based on the difference between the target and current angle of attack
        throttle_adjustment = pid_acceleration(target_acceleration - current_acceleration)

        self.variables["aoa"] += throttle_adjustment * 1e-2



            # return aoa_adjustment



        # df = pd.DataFrame([_alpha_vect, accel, drag])
        # df = df.T
        # df.columns = ['alpha', 'accel', 'drag']
        # df = df[df.accel > 0]
        # ind_min = df.accel.idxmin()
        #
        # _min_accel = df.loc[ind_min, 'accel']
        #
        # _aoa = df.loc[ind_min, 'alpha']
        #
        # _drag = df.loc[ind_min, 'drag']
        #
        # return _min_accel, _aoa, _drag

    def normal_climb(self):

        v = self.variables["v"]
        rho, S = self.constants_dict["rho"], self.constants_dict["S"]
        # print(self.variables["lift"], self.constants_dict["weight"])
        # assert self.variables["lift"] == self.constants_dict["weight"]
        # self.variables["aoa"] = self.constants_dict["teta_target"]
        # self.variables["drag"] = 0.5 * rho * S * self.f_cd(self.variables["aoa"]) * v ** 2

        # TODO: in order to reach constant climb -> W=L
        cl = 2 * self.constants_dict["weight"] / (rho * S * v ** 2)
        cd = self.drag_polar(cl ** 2)

        self.variables["aoa"] = self.f_a_cd(cd)

        self.variables["gamma_rad"] = (self.variables["thrust"] * np.cos(uc.deg2rad(self.variables["aoa"])) /
                                       self.constants_dict["weight"]) - (cd / cl)
        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])

        self.variables["teta"] = self.variables["gamma"] + self.variables["aoa"]

        print(f'V2min reached @ {round(float(self.variables["height"]), 2)}ft!')
        #
        while self.variables["height"] < 35.:
            self.variables["dv"] = 0
            self.variables["dx"] = self.variables["v"] * np.cos(self.variables["gamma_rad"]) * self.variables["dt"]

            #                 if _bool:
            dh = self.variables["v"] * np.sin(self.variables["gamma_rad"]) * self.variables["dt"]
            # update
            self.variables["height"] += uc.m2ft(dh)
            super().update_values()

    def calculate_equations_of_motion(self):

        aoa, gamma = self.variables['aoa'], self.variables['gamma']
        thrust, v = self.variables["thrust"], self.variables["v"]
        rho, S = self.constants_dict["rho"], self.constants_dict["S"]
        weight, mass = self.constants_dict['weight'], self.constants_dict['mass']

        T_aero = thrust * np.cos(np.radians(aoa))

        drag = 0.5 * rho * S * self.f_cd(aoa) * v ** 2
        
        lift = 0.5 * rho * S * self.f_cl(aoa) * v ** 2

        SFx = T_aero - drag - weight * np.sin(np.radians(gamma))

        accel = SFx / mass

        forces = {
            "Drag": drag,
            "Lift": lift,
            "F_x": SFx,
        }

        return forces, accel

    def update_aerial_values(self):

        thrust = self.variables["thrust"]
        weight = self.constants_dict["weight"]

        # advanced ac perfo p306
        dgamma = self.variables["dt"] * (thrust * np.sin(np.radians(self.variables["aoa"])) + self.variables["lift"] - weight) \
                 / (self.variables["mass"] * self.variables["v"])

        self.variables["gamma_rad"] += dgamma

        self.variables["gamma"] = np.degrees(self.variables["gamma_rad"])

        # assert self.variables["gamma"] < 3, 'Gamma overshoots. Increase VR'

        self.variables["teta"] = self.variables["aoa"] + self.variables["gamma"]

        self.variables["v_z"] =(self.variables["thrust"] - self.variables["drag"])/ weight * self.variables["v"] - \
                               self.variables["v"]/9.81*self.variables["accel"]
            #self.variables["v"] * np.sin(self.variables["gamma_rad"])

        self.variables["dv"] = self.variables["accel"] * self.variables["dt"]

        self.variables["dx"] = self.variables["v"] * np.cos(self.variables["gamma_rad"]) * self.variables["dt"]
        # dh * np.tan(gamma_rad)

        dh = self.variables["v"] * np.sin(self.variables["gamma_rad"]) * self.variables["dt"]  # a_z*dt**2 + v_z * dt

        # update
        self.variables["height"] += uc.m2ft(dh)

    def airborne_phase(self):

        # super().transition_phase()
        # TODO: Study the climb part p299 the Vz is low
        self.reached_v2 = False
        while self.variables["height"] < 35.:

            # control alpha in order to get 0 accel
            self.control_accel()

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

            self.update_aerial_values()

            super().update_values()

            if self.variables["v_kt"] >= self.speeds["v_target"]:
                if not self.reached_v2:
                    self.characteristic_instants["v2"] = {"Instant": self.variables["t"],
                                                          "Speed": self.variables["v_kt"]}
                    self.reached_v2 = True
                self.normal_climb()

            self.characteristic_instants["v35ft"] = {"Instant": self.variables["t"],
                                                          "Speed": self.variables["v_kt"]}


if __name__ == "__main__":
    mass = 80000.0  # [kg]
    conf = "1+F"
    zp = 0.0  # [ft]
    lg = "Up"
    engine_state = "OEI"

    a320_to = AirbornePhase(mass, conf, zp, lg, engine_state)
    a320_to.airborne_phase()