import random


def altimeter(z):
    
    # MS5611 barometer has an accuracy of +- 1.5 mbar which translates roughly to +- 12 m
    # Will use this range as simulated sensor noise
    
    # starting values were +- 1200

    random_decimal = random.randrange(-1200,1200)/100
    #random_decimal = random.randrange(-10000,10000)/100

    return z + random_decimal