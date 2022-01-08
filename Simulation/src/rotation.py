import numpy as np


def yaw(phe):
    # rotation about the z-axis, takes in angle phe in the function (rad)
    return np.array([[np.cos(phe), -np.sin(phe), 0],[np.sin(phe), np.cos(phe), 0],[0, 0, 1]])

def pitch(theta):
    # rotation about the y-axis, takes in angle theta in the function (rad)
    return np.array([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0],[-np.sin(theta), 0, np.cos(theta)]])

def roll(phi):
    # rotation about the x-axis, take sin angle phi in the function (rad)
    return np.array([[1, 0, 0],[0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])

def body_aero(velocity_body):
    #TODO: Disperse comments to properly explain the lines of code (atmosphere.py for reference)
    # Takes in np.array velocity_body 
    # rotation matrix that transforms the body frame to the aerodyanmic frame 
    # takes in beta (side-slip angle), and alpha (angles of attack)
    # needs the velocities in the body frame to calculate beta and alpha
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[0][1]
    vz_b = velocity_body[0][2]
    
    alpha = np.arctan2(vz_b,vx_b)
    beta = np.arctan2(vy_b, np.sqrt(vx_b**2 + vz_b**2))
    
    R_ba = np.array([[np.cos(beta)*np.cos(alpha), np.sin(beta), np.cos(beta)*np.sin(alpha)], 
                     [-np.sin(beta)*np.cos(alpha), np.cos(beta), -np.sin(beta)*np.sin(alpha)],
                     [-np.sin(alpha), 0, np.cos(alpha)]])
    return R_ba