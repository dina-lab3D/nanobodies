import numpy as np
from tensorflow.keras.models import load_model
import tensorflow as tf
import pickle
import sys
import seaborn as sns
import os
import matplotlib.pyplot as plt
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


def plot_test(trained_model, test_X, test_Y):

    files_name = "covid2"
    names = ["distance", "omega", "theta", "phi"]
    pred_Y = trained_model.predict(test_X)

    if not os.path.exists("test_vs_pred_" + files_name):
        os.mkdir("test_vs_pred_" + files_name)
    for i in range(len(test_X)):

        fig, axes = plt.subplots(4, 2, figsize=(25,25))
        for j in range(len(names)):
            if names[j] == "distance":
                sns.heatmap(test_Y[i][j,:,:,0], ax=axes[j][0],xticklabels=False, yticklabels=False, cbar=False)
                sns.heatmap(pred_Y[j][i,:,:,0], ax=axes[j][1],xticklabels=False, yticklabels=False, cbar=False)
            else:
                sns.heatmap(test_Y[i][j,:,:,1], ax=axes[j][0],xticklabels=False, yticklabels=False, cbar=False)  # sin
                sns.heatmap(pred_Y[j][i,:,:,1], ax=axes[j][1],xticklabels=False, yticklabels=False, cbar=False)  # sin
                # sns.heatmap(np.arctan2(test_Y[i][j,:,:,0],test_Y[i][j,:,:,1]), ax=axes[j][0],xticklabels=False, yticklabels=False, cbar=False)
                # sns.heatmap(np.arctan2(pred_Y[j][i,:,:,0],pred_Y[j][i,:,:,1]), ax=axes[j][1],xticklabels=False, yticklabels=False, cbar=False)

            axes[j][0].set_title(names[j] + " true")
            axes[j][1].set_title(names[j] + " predicted")

        plt.savefig(os.path.join("test_vs_pred_" + files_name, "test_" + str(i)))
        plt.clf()



###################################    for getting training data  #####################################################


# with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/nn_input_3_no_clip.pkl", "rb") as input_file:
#     X = pickle.load(input_file)
# with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/nn_labels_3_no_clip.pkl", "rb") as feature_file:
#     Y = pickle.load(feature_file)
#     # Y = reshape_y(Y)
# with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/pdb_names_3_no_clip.pkl", "rb") as names_file:
#     pdb_names = pickle.load(names_file)
# with open("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/test_names_NanoNet_model.pkl",  "rb") as test_file_names:
#     old_test_names = pickle.load(test_file_names)
#
#
# # model = tf.keras.models.load_model("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_model")
#
#
# train_X = []
# train_Y = []
# train_names = []
#
# test_X = []
# test_Y = []
# test_names = []
#
# for i in range(len(pdb_names)):
#     if pdb_names[i] not in old_test_names:
#         train_X.append(X[i])
#         train_Y.append(Y[i])
#         train_names.append(pdb_names[i])
#     else:
#         test_X.append(X[i])
#         test_Y.append(Y[i])
#         test_names.append(pdb_names[i])
#
# train_X = np.stack(train_X, axis=0)
# train_Y = np.stack(train_Y, axis=0)
#
# test_X = np.stack(test_X, axis=0)
# test_Y = np.stack(test_Y, axis=0)
#
#
# name = "3_no_clip.pkl"
#
# train_label_file = "train_Y_{}".format(name)
# train_input_file = "train_X_{}".format(name)
# train_pdb_names_file = "train_names_{}".format(name)
#
# test_label_file = "test_Y_{}".format(name)
# test_input_file = "test_X_{}".format(name)
# test_pdb_names_file = "test_names_{}".format(name)
#
# pickle.dump(train_Y, open(train_label_file, "wb"))
# pickle.dump(train_X, open(train_input_file, "wb"))
# pickle.dump(np.array(train_names), open(train_pdb_names_file, "wb"))
#
# pickle.dump(test_Y, open(test_label_file, "wb"))
# pickle.dump(test_X, open(test_input_file, "wb"))
# pickle.dump(np.array(test_names), open(test_pdb_names_file, "wb"))



########################################################################################################################

with open("/cs/usr/tomer.cohen13/lab/NN/CovidTest/nn_input_3_covid.pkl", "rb") as input_file:
    X = pickle.load(input_file)
with open("/cs/usr/tomer.cohen13/lab/NN/CovidTest/nn_labels_3_covid.pkl", "rb") as feature_file:
    Y = pickle.load(feature_file)
    Y_for_plot = np.copy(Y)
    Y = reshape_y(Y)
with open("/cs/usr/tomer.cohen13/lab/NN/CovidTest/pdb_names_3_covid.pkl", "rb") as names_file:
    pdb_names = pickle.load(names_file)


model2 = tf.keras.models.load_model("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_model")
model = tf.keras.models.load_model("/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/DIM_3_best_test_1")

loss = []
loss2 = []

dis = []
omega = []
theta = []
phi = []
names = []


for i in range(len(pdb_names)):

    if "RBD" in pdb_names[i].upper():
        continue
    s = model.evaluate(np.array([X[i]]), [Y[0][i,:,:,:].reshape(1,32,32,1), Y[1][i,:,:,:].reshape(1,32,32,2), Y[2][i,:,:,:].reshape(1,32,32,2), Y[3][i,:,:,:].reshape(1,32,32,2)])
    s2 = model2.evaluate(np.array([X[i]]), [Y[0][i,:,:,:].reshape(1,32,32,1), Y[1][i,:,:,:].reshape(1,32,32,2), Y[2][i,:,:,:].reshape(1,32,32,2), Y[3][i,:,:,:].reshape(1,32,32,2)])

    loss2.append(s2[0])
    loss.append(s[0])
    dis.append(s[1])
    omega.append(s[2])
    theta.append(s[3])
    phi.append(s[4])



# plot_test(model2, X, Y_for_plot)

print("############")
print("total loss: {}".format(np.mean(loss)))
print(loss)
print(loss2)
print(np.mean(loss2))
print("#############")
print("distance loss: {}".format(np.mean(dis)))
print(dis)
print("#############")
print("omega loss: {}".format(np.mean(omega)))
print(omega)
print("#############")
print("theta loss: {}".format(np.mean(theta)))
print(theta)
print("#############")
print("phi loss: {}".format(np.mean(phi)))
print(phi)
print("#############")