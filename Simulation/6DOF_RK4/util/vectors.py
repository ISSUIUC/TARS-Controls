import numpy as np

def norm(x) -> np.ndarray:
    '''
    args
    nd vector

    returns normalized vector or 0 vector if magnitude is 0
    '''
    norm = np.linalg.norm(x)
    if (norm != 0):
        return x.copy()/norm
    return np.zeros(np.shape(x))

def body_to_world(roll, pitch, yaw, body_vector):
    roll = np.array([[1, 0, 0], [0, np.cos(roll), -np.sin(roll)], [0, np.sin(roll), np.cos(roll)]])
    pitch = np.array([[np.cos(pitch), 0, np.sin(pitch)], [0, 1, 0], [-np.sin(pitch), 0, np.cos(pitch)]])
    yaw = np.array([[np.cos(yaw), -np.sin(yaw), 0], [np.sin(yaw), np.cos(yaw), 0], [0, 0, 1]])
    return yaw @ pitch @ roll @ body_vector

def world_to_body(roll, pitch, yaw, world_vector):
    roll = np.array([[1, 0, 0], [0, np.cos(roll), -np.sin(roll)], [0, np.sin(roll), np.cos(roll)]])
    pitch = np.array([[np.cos(pitch), 0, np.sin(pitch)], [0, 1, 0], [-np.sin(pitch), 0, np.cos(pitch)]])
    yaw = np.array([[np.cos(yaw), -np.sin(yaw), 0], [np.sin(yaw), np.cos(yaw), 0], [0, 0, 1]])
    return (yaw @ pitch @ roll).T @ world_vector