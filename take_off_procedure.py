from src.modules.airborne_phase import AirbornePhase
from src.modules.rotation_phase import RotationPhase
from src.modules.ground_phase import GroundRoll


class TakeOffProc(GroundRoll, RotationPhase, AirbornePhase):

    def check_vr(self):

        if self.variables["need_to_increase_Vr"]:
            self.speeds["vr"] += 0.01 * self.speeds["v_stall"]
            self.speeds["vef"] = self.speeds["vr"] - 2
            print(f"VR increased at {self.speeds['vr']} kt")
            super().initialize_data()

    def takeoff(self):

        super().pilot_preparation()  # calculate all the characteristic speeds
        print(f"Vr = {self.speeds['vr']} kt")
        while self.variables["cas_kt"] < self.speeds["v_target"]:
            super().up_to_rotation()  # 1st segment of takeoff
            super().transition_phase()  # 2nd segment of takeoff
            self.check_vr()
            if not self.variables["need_to_increase_Vr"]:
                super().airborne_phase()  # 3rd segment of takeoff
                self.check_vr()


if __name__ == "__main__":

    input_variables = {
        "mass": 80e3,  # [kg]
        "conf": "1+F",
        "zp": 0,  # [ft]
        "lg": "Up",
        "engine_state": "OEI",
        "timestep": 1e-1,
    }
    a320_to = TakeOffProc(input_variables)
    a320_to.takeoff()
