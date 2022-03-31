import random


def alt_noise(z):
    
    # MS5611 barometer has an accuracy of +- 1.5 mbar which translates roughly to +- 12 m
    # Will use this range as simulated sensor noise
    
    # starting values were +- 1200
    
    alt_range = 1200
    random_decimal = random.randrange(-alt_range,alt_range)/100
    #random_decimal = random.randrange(-10000,10000)/100

    return z + random_decimal

def h_predicted(pos, vel):
    h_predict = (0.5*vel**2 + 9.81*pos)/9.81
    c = 0.0000095
    y_predict_update = h_predict  + c*(vel**2)*(9144 - h_predict)

    return y_predict_update
