import properties.properties as prop


class Controller():

    def __init__(self, Kp, Ki, Kd, timestep, des_apogee):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.error_sum = 0
        self.error_prev = 0
        self.dt = timestep
        self.des_apogee = des_apogee
        # self.prop = prop

    def get_flap_extension(self, pred_apogee):
        error_curr = pred_apogee - self.des_apogee
        self.error_sum += error_curr*self.dt
        error_dt = (error_curr - self.error_prev) / self.dt
        flap_ext = self.Kp * error_curr + self.Ki * self.error_sum + self.Kd * error_dt
        self.error_prev = error_curr
        if (flap_ext > prop.max_ext_length):
            return prop.max_ext_length
        elif (flap_ext < 0):
            return 0
        return flap_ext
