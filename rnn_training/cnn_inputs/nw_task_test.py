#!/usr/bin/env python3
from psychrnn.tasks.task import Task
import numpy as np

class uniform_curvature_task(Task):
    def __init__(self, dt=10, tau=20, T=1800, N_in=259, N_out=1, N_batch=600, in_noise=0.0, alpha=1):
        # ----------------------------------
        # Define network parameters
        # ----------------------------------
        super().__init__(N_in, N_out, dt, tau, T, N_batch)
        self.in_noise = in_noise
        self.act = np.random.rand(100,256)
        self.cc = np.random.rand(100,1)
                
    def generate_trial_params(self, batch, trial):
        # ----------------------------------
        # Define parameters of a trial
        # ----------------------------------
        condNum, stimNum = divmod(trial,100)

        params = dict()
        params['s'] = self.cc[stimNum]
        params['act'] = self.act[stimNum,:]

        if condNum == 0:
            params['a'] = -20
            params['b'] = 100
        elif condNum == 1:
            params['a'] = 0
            params['b'] = 100
        elif condNum == 2:
            params['a'] = 20
            params['b'] = 100
        elif condNum == 3:
            params['a'] = -20
            params['b'] = 140
        elif condNum == 4:
            params['a'] = 0
            params['b'] = 140
        elif condNum == 5:
            params['a'] = 20
            params['b'] = 140

        # if trial % 120<40:
        #     params['a'] = -20
        # elif trial % 120<80:
        #     params['a'] = 0
        # else:
        #     params['a'] = 20

        # if trial % 120<20:
        #     params['b'] = 100
        # elif trial % 120<40:
        #     params['b'] = 140
        # elif trial % 120<60:
        #     params['b'] = 100
        # elif trial % 120<80:
        #     params['b'] = 140
        # elif trial % 120<100:
        #     params['b'] = 100
        # else:
        #     params['b'] = 140

        params['o'] = round(params['s']*params['b'] - params['b']/2 + params['a'])
        
        params['fix_onset'] = 200
        params['s_onset'] = 400
        params['ab_onset'] = 1000
        params['fix_offset'] = 1300
        
        return params

    def trial_function(self, time, params):
        """ Compute the trial properties at the given time.
    
        Based on the params compute the trial stimulus (x_t), correct output (y_t), and mask (mask_t) at the given time.
    
        Args:
            time (int): The time within the trial (0 <= time < T).
            params (dict): The trial params produced generate_trial_params()
    
        Returns:
            tuple:
    
            x_t (ndarray(dtype=float, shape=(N_in,))): Trial input at time given params.
            y_t (ndarray(dtype=float, shape=(N_out,))): Correct trial output at time given params.
            mask_t (ndarray(dtype=bool or float, shape=(N_out,))): True or 1 if the network should train to match the y_t, False or 0 if the network should ignore y_t when training. Can also be an arbitrary scaling value.
        """
        
        # ----------------------------------
        # Initialize with input noise
        # ----------------------------------

        x_t = np.zeros(self.N_in)
        y_t = np.zeros(self.N_out)
        mask_t = np.ones(self.N_out) # weighing all timepoints and outputs equally
        
        # ----------------------------------
        # Compute values
        # ----------------------------------
        
        if time > params['fix_onset'] and time < params['fix_offset']:
            x_t[258] += 1
        
        if time > params['s_onset']:
            x_t[np.arange(256)] += params['act']
        
        if time > params['ab_onset']:
            x_t[256] += params['a']
            x_t[257] += params['b']
        
        if time > params['fix_offset']:
            y_t[0] = params['o']
            # y_t[1] = params['s']
        
        return x_t, y_t, mask_t
