from filterpy.common import Q_continuous_white_noise
import numpy as np

s_dt = 0
x_k = np.zeros([2,1])
F = np.zeros([2,2])
H = np.zeros([2,1])
B = np.zeros([2,1])
P_k = np.zeros([2,2])
Q = 0
R = np.zeros([1])
x_priori = np.zeros([2,1])
P_priori = np.zeros([2,2])

# Dictionary for Kalman Filter Estimations
kalman_dic = {
        "alt": [],
        "vel": [],
        "accel": []
}

def initialize(pos_f, vel_f, accel_f, time_step):
    global s_dt, x_k, F, H, B, P_k, Q, R

    s_dt = time_step
    x_k = np.array([[pos_f],
                    [vel_f],
                    [accel_f]])
    
    # F is the state space 'A' matrix
    F = np.array([[1.0, s_dt, (s_dt**2) / 2],
                  [0.0, 1.0, s_dt],
                  [0.0, 0.0, 1.0]])
    
    # H is converstion between measurement and state (analogous to C)
    H = np.array([[1.0,0.0,0.0],
                  [0.0,0.0,1.0]]) 
    
    # Standard State-Space 'B' Matrix
    # B = np.array([[1.0],
    #               [0.0]])

    # Covariance [P]
    # Check Steady State -> at end of flight
    P_k = np.array([[.018,0.009, 0.005],
                    [0.009,0.009, 0.0045],
                    [0.005, 0.0045, 10]])

    # White Noise [Q]
    Q = Q_continuous_white_noise(dim=3, dt=s_dt, spectral_density=.00899)

    # Measurement Noise Function [R]
    # High-G Accel: 49 - 195 m-Gs accuracy
    # Low-G Accel: .25 m-G accuracy
    # MUST BE SQUARE
    R = np.array([[12.0,0],
                  [0,1.4]])

# Set priori state (guess of next step)
def priori(u):
    global x_priori, P_priori
    
    # x_priori = (F @ x_k) + ((B @ u).T) #* For some reason doesnt work when B or u is = 0
    x_priori = F @ x_k
    P_priori = (F @ P_k @ F.T) + Q

# Update Kalman Gain, posteriori state (guess of current step with new data), Covariance update
def update(pos_f, accel_f, Sref_a, rho):
    global K, x_k, P_k, F
    
    # Update Kalman Gain
    if (len(R) == 1):
        K = P_priori @ H.T * np.reciprocal(H @ P_priori @ H.T + R)
    else:
        K = (P_priori @ H.T) @ np.linalg.inv(H @ P_priori @ H.T + R)   
             
    # Sensor Measurements
    y_k = np.array([[pos_f], [accel_f]])
    
    # Posteriori Update
    x_k = x_priori + K @ (y_k - H @ x_priori)
    P_k = (np.eye(len(K)) - K@H) @ P_priori

    # Add to Kalman Dictionary
    kalman_dic["alt"].append(x_k[0][0])
    kalman_dic["vel"].append(x_k[1][0])
    kalman_dic["accel"].append(x_k[2][0])
    
    # Extended KF
    # F[1][1] = 1 - (Sref_a*rho*0.58*x_k[0][1] * s_dt)