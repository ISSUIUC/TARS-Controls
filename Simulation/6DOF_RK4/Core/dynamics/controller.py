import numpy as np
import properties.properties as prop
import properties.data_loader as dataloader

# Load desired config file
config = dataloader.config

class Controller():
    """Flap extension controller class for the rocket
    
    Args:
        Kp (float): proportional gain
        Ki (float): integral gain
        Kd (float): derivative gain
        timestep (float): timestep in seconds
        des_apogee (float): desired apogee in meters
    """
    def __init__(self, Kp, Ki, Kd, timestep, des_apogee, stage_config):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.error_sum = 0
        self.error_prev = 0
        self.dt = timestep
        self.des_apogee = des_apogee
        self.prev_flap = 0
        self.stage_config = stage_config

    def get_flap_extension(self, control, pred_apogee):
        """Returns the flap extension in meters given the predicted apogee in meters
        
        Args:
            control (boolean): is control allowed
            pred_apogee (float): predicted apogee in meters
            
        Returns:
            float: flap extension in meters
        """
        if control:
            error_curr = pred_apogee - self.des_apogee
            self.error_sum += error_curr*self.dt
            error_dt = (error_curr - self.error_prev) / self.dt
            flap_ext = self.Kp * error_curr + self.Ki * self.error_sum + self.Kd * error_dt
            self.error_prev = error_curr

            flap_ext += np.sign(flap_ext - self.prev_flap)*min(abs((flap_ext - self.prev_flap)/self.dt), self.stage_config["flaps"]["max_ext_spd"])*self.dt
            
            if (flap_ext > self.stage_config["flaps"]["max_ext_length"]):
                flap_ext =  self.stage_config["flaps"]["max_ext_length"]
            elif (flap_ext < 0):
                flap_ext = 0
                
            self.prev_flap = flap_ext
        else:
            flap_ext = 0   
        return flap_ext