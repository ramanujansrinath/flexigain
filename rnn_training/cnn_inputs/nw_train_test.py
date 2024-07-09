# This is more or less the same as the one with simple inputs. 
# The difference is that instead of generating a random 
# curvature values between 0 and 1 for each trial, 
# this rnn takes in a 256 dimensional vector which corresponds 
# to the activations of pool4 of VGG16 after dimensionality 
# reduction. Each trial now is an image and the RNN has to 
# figure out the curvature from CNN activations. 
# There are two sets of inputs during training. 
#     - 1000 random shapes
#     - 50 random shapes with 20 curvatures each
# Testing both of these training types is 5 shapes with 20 curvatures each.

from numpy import genfromtxt
from nw_task_init import curvature_saccade_task
from nw_task_test import uniform_curvature_task
from psychrnn.backend.models.basic import Basic
from psychrnn.backend.simulation import BasicSimulator
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import scipy as sp

def train_model(model_id):
  model_name = 'model' + str(model_id)
  N_trainbatch = 200 # number of trials per training update

  task = curvature_saccade_task(N_batch = N_trainbatch)
  task.act = genfromtxt('img_gen/img_data_2/act.csv', delimiter=',')
  task.cc = genfromtxt('img_gen/img_data_2/cc.csv', delimiter=',')

  network_params = task.get_task_params()
  network_params['name'] = model_name
  network_params['N_rec'] = 512 # number of hidden units
  network_params['rec_noise'] = 0.01 # recurrent noise

  # "biological constraints" (optional):
  network_params['autapses'] = True # whether or not hidden units can connect to themselves
  network_params['dale_ratio'] = None # dale ratio, e.g. 0.8 if 80% of units can only send excitatory projections and 20% can only send inhibitory projections

  # regularization (optional)
  network_params['L2_in'] = 0
  network_params['L2_rec'] = 0
  network_params['L2_out'] = 0
  network_params['L2_firing_rate'] = 0
  network_params['L1_in'] = 0
  network_params['L1_rec'] = 0
  network_params['L1_out'] = 0

  model = Basic(network_params)

  train_params = dict() # specify training parameters. default transfer function is ReLU (can change if desired).
  train_params['training_iters'] = 1000000
  train_params['learning_rate'] = 0.0001

  # TRAIN
  losses, initialization_time, training_time = model.train(task, train_params)
  model.save('weights/' + model_name)
  model.destruct()

  plt.figure(figsize=(5,5))
  plt.plot(losses)
  plt.title('Loss during training')
  plt.ylabel('Minibatch loss')
  plt.xlabel('Batch number')
  plt.savefig('plots/' + model_name + '_trainingLoss', dpi=200)

def test_model(model_id):
  model_name = 'model' + str(model_id)
  N_testbatch = 600 # number of trials to test
  # task = curvature_saccade_task(N_batch = N_testbatch)
  task = uniform_curvature_task(N_batch = N_testbatch)
  task.act = genfromtxt('img_gen/img_data_2/act_test.csv', delimiter=',')
  task.cc = genfromtxt('img_gen/img_data_2/cc_test.csv', delimiter=',')

  network_params = task.get_task_params()
  network_params['N_rec'] = 512 # number of hidden units
  network_params['rec_noise'] = 0.0 # recurrent noise
  test_inputs, target_outputs, mask, trial_params = task.get_trial_batch()
  simulator = BasicSimulator(weights_path='weights/' + model_name + '.npz', params=network_params)

  outputs, state_vars = simulator.run_trials(test_inputs)
  
  # save data in matlab format
  # collect arrays in dictionary
  savedict = {
      'outputs' : outputs,
      'state_vars' : state_vars,
      'test_inputs' : test_inputs,
      'target_outputs' : target_outputs,
      'trial_params' : trial_params
  }
  sp.io.savemat('data/test_' + model_name + '.mat', savedict)

if __name__ == '__main__':
  # train_model(1)
  test_model(1)
  # p = Pool(processes=cpu_count())
  # p.map(train_model, range(cpu_count()))
  # p.map(test_model, range(cpu_count()))