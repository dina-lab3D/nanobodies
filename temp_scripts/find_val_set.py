import numpy as np
from tensorflow.keras.models import load_model
import tensorflow as tf
import pickle

def reshape_y(y):
    return [y[:,0,:,:,0].reshape(-1,32,32,1), y[:,1,:,:,:], y[:,2,:,:,:], y[:,3,:,:,:]]


with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/nn_input_rev.pkl", "rb") as input_file:
    X = pickle.load(input_file)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/nn_labels_rev.pkl", "rb") as feature_file:
    Y = pickle.load(feature_file)
    Y = reshape_y(Y)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/pdb_names_rev.pkl", "rb") as names_file:
    pdb_names = pickle.load(names_file)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/test_names_NanoNet_model.pkl",  "rb") as test_file_names:
    test_names = pickle.load(test_file_names)


model = tf.keras.models.load_model("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_model")

loss = np.zeros(len(pdb_names) - len(test_names))
dis = np.zeros(len(pdb_names)- len(test_names))
omega = np.zeros(len(pdb_names)- len(test_names))
theta = np.zeros(len(pdb_names)- len(test_names))
phi = np.zeros(len(pdb_names)- len(test_names))
names = np.zeros(len(pdb_names)- len(test_names))
k = 0

for i in range(len(pdb_names)):
    if pdb_names[i] in test_names:
        continue
    s = model.evaluate(np.array([X[i]]), [Y[0][i,:,:,:].reshape(1,32,32,1), Y[1][i,:,:,:].reshape(1,32,32,2), Y[2][i,:,:,:].reshape(1,32,32,2), Y[3][i,:,:,:].reshape(1,32,32,2)])
    loss[k] = s[0]
    dis[k] = s[1]
    omega[k] = s[2]
    theta[k] = s[3]
    phi[k] = s[4]
    k += 1

sort_index = np.argsort(omega)

loss = loss[sort_index][::-1]
dis = dis[sort_index][::-1]
omega = omega[sort_index][::-1]
theta = theta[sort_index][::-1]
phi = phi[sort_index][::-1]

print(np.mean(loss[0:161]))
print(np.mean(dis[0:161]))
print(np.mean(omega[0:161]))
print(np.mean(theta[0:161]))
print(np.mean(phi[0:161]))