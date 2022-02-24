import random


def alt_noise(z):
    
    # MS5611 barometer has an accuracy of +- 1.5 mbar which translates roughly to +- 12 m
    # Will use this range as simulated sensor noise
    
    # starting values were +- 1200
    
    alt_range = 12000
    random_decimal = random.randrange(-alt_range,alt_range)/100
    #random_decimal = random.randrange(-10000,10000)/100

    return z + random_decimal