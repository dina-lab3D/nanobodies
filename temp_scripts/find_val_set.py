import numpy as np
from tensorflow.keras.models import load_model
import tensorflow as tf
import pickle
import sys
sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/nano_buddies')
from NanoNetUtils import *





# dir  = "/cs/labs/dina/tomer.cohen13/NN/RosettaFasta"
# pdb = "2EH7_1"
#
#
# x = generate_input("{}/{}/{}.fa".format(dir,pdb,pdb))
# y = generate_label("{}/{}/{}.fa".format(dir,pdb,pdb), "{}/{}/ref.pdb".format(dir,pdb))
# model = tf.keras.models.load_model("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_model")
# s = model.evaluate(np.array([x]), [y[0, :, :, 0].reshape(1, 32, 32, 1),
#                                       y[1, :, :, :].reshape(1, 32, 32, 2),
#                                       y[2, :, :, :].reshape(1, 32, 32, 2),
#                                       y[3, :, :, :].reshape(1, 32, 32, 2)])
# print(s)
# exit()

def reshape_y(y):
    return [y[:,0,:,:,0].reshape(-1,32,32,1), y[:,1,:,:,:], y[:,2,:,:,:], y[:,3,:,:,:]]


with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/nn_input_rev.pkl", "rb") as input_file:
    X = pickle.load(input_file)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/nn_labels_rev.pkl", "rb") as feature_file:
    Y = pickle.load(feature_file)
    Y = reshape_y(Y)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/pdb_names_rev.pkl", "rb") as names_file:
    pdb_names = pickle.load(names_file)
with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/test_names_NanoNet_model.pkl",  "rb") as test_file_names:
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
    if k == 50:
        break
    if not np.array_equal(X[i],generate_input("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/{}/{}.fa".format(pdb_names[i], pdb_names[i].split("/")[1]))):
        print(pdb_names[i])

    # s = model.evaluate(np.array([X[i]]), [Y[0][i,:,:,:].reshape(1,32,32,1), Y[1][i,:,:,:].reshape(1,32,32,2), Y[2][i,:,:,:].reshape(1,32,32,2), Y[3][i,:,:,:].reshape(1,32,32,2)])
    # loss[k] = s[0]
    # dis[k] = s[1]
    # omega[k] = s[2]
    # theta[k] = s[3]
    # phi[k] = s[4]
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