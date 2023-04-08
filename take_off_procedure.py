from src.modules.airborne_phase import AirbornePhase
from src.modules.rotation_phase import RotationPhase
from src.modules.ground_phase import GroundRoll


class TakeOffProc(GroundRoll, RotationPhase, AirbornePhase):

    def takeoff(self):

        super().pilot_preparation()
        print(f"Vr = {self.speeds['vr']} kt")
        while self.variables["cas_kt"] < self.speeds["v_target"]:

            super().up_to_rotation()
            super().transition_phase()
            if self.variables["need_to_increase_Vr"]:
                self.speeds["vr"] += 0.01 * self.speeds["v_stall"]
                self.speeds["vef"] = self.speeds["vr"] - 2
                print(f"VR increased at {self.speeds['vr']} kt")
                super().initialize_data()
            else:
                super().airborne_phase()


if __name__ == "__main__":

    mass = 80000.0  # [kg]
    conf = "1+F"
    zp = 0.0  # [ft]
    lg = "Up"
    engine_state = "OEI"
    dt = 0.1
    a320_to = TakeOffProc(mass, conf, zp, lg, engine_state, dt)
    a320_to.takeoff()
