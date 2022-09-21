import numpy.random as random
import src.atmosphere as atmosphere


def alt_noise(alt_input, mach):
    
    # MS5611 barometer has an accuracy of ± 1.5 mbar (150Pa) which translates roughly to ± 12 m
    # 90% of test data falls within ± 1.5 mbar - calculate standard deviation associated with 90% (z = 1.645)
    # 150Pa = z*s/sqrt(n) where n is number of tests in data set
    # s = 150 / 1.645
    # Lembeck - thermal vacuum chamber for testing
    # Will use this range as simulated sensor noise
    
    # starting values were ± 1200
    pressure_range = 150    # Pa
    z = 1.645
    pressureRat = 1
    # add pressure ratio due to supersonic flight (may not be necessary based on hole locations for barometer)
    # if(mach >= 1):
    #     pressureRat = atmosphere.pressure_ratio(mach)

    pressure = atmosphere.pressure(alt_input) * pressureRat + random.normal(scale = pressure_range/z)
    # print(mach,  atmosphere.pressure(alt_input), pressure)
    alt = atmosphere.alt_from_pressure(pressure)
    # random_decimal = random.randrange(-alt_range,alt_range)/100
    # random_decimal = random.randrange(-10000,10000)/100
    return alt