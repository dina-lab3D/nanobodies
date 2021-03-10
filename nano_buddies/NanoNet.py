
import tensorflow as tf
from tensorflow.keras import layers
import pickle
import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.model_selection import train_test_split
import seaborn as sns
import os


def swish(x, beta=1):
    return (x*tf.keras.backend.sigmoid(beta*x))


tf.keras.utils.get_custom_objects().update({'swish': layers.Activation(swish)})

# idea - add resnets in the end?
# idea - calculate tanh of omega after combining with transpose (same for dist?)


CDR_DICT = {3:1, 1:0, 2:2}
DIM = 1  # normal is 3
CDR = 3  # normal is 3
KERNELS = 32  # normal is 32
RESNET_BLOCKS = 3  # normal is 3
RESNET_SIZE = (17, 17)  # normal is (17,17)
FIRST_RESNET_SIZE = (17, 17)  # normal is (17,17)
DIALETED_RESNET_BLOCKS = 5  # normal is 5
DIALETION = [1,2,4,8,16]  # normal is [1,2,4,8,16]
DIALETED_RESNET_SIZE = (5, 5)  # normal is (5,5)
EPOCHS = 175  # normal is 150
LR = 0.0005  # normal is 0.0005
TEST_SIZE = 50/2182  # 50 nano-bodies
VAL_SIZE = 0.075  # 150 nano-bodies (0.075)
BATCH = 32  # normal is 32
DROPOUT = 0.2  # normal is 0.2
END_CONV_SIZE = 4   # normal is 4
END_CONV_KER = (5, 5)  # normal is (5,5)
DIALETED_RESNET_KERNELS = 64  # normal is 64
ACTIVATION = "relu"  # normal relu
END_ACTIVATION = "elu"  # normal elu
LOSS = "mse"
BINS = False
POOL = False
files_name = "DIM_1_no_clip_1"


def reshape_y(y):
    return [y[:,0,:,:,0].reshape(-1,32,32,1), y[:,1,:,:,:], y[:,2,:,:,:], y[:,3,:,:,:]]


def d2_net_architecture():
    """

    :return:
    """
    if DIM == 1:
        input_layer = tf.keras.Input(shape=(32, 21, 1), name="InputLayer")
    else:
        input_layer = tf.keras.Input(shape=(32, 21, 3), name="InputLayer")

    loop_layer = layers.Conv2D(KERNELS, FIRST_RESNET_SIZE, padding='same')(input_layer)

    # first res block: 2d_conv -> relu -> 2d_conv -> batch_normalization -> add -> relu
    for i in range(RESNET_BLOCKS):
        conv_layer = layers.Conv2D(KERNELS, RESNET_SIZE, activation=ACTIVATION, padding='same')(loop_layer)
        conv_layer = layers.Conv2D(KERNELS, RESNET_SIZE, padding='same')(conv_layer)
        batch_layer = layers.BatchNormalization()(conv_layer)
        loop_layer = layers.Add()([batch_layer, loop_layer])
        loop_layer = layers.Activation(ACTIVATION)(loop_layer)

    premut_layer = layers.Permute((1,3,2))(loop_layer)

    # for 32*32*64 matrix
    conv_layer = layers.Conv2D(KERNELS, (5,5), activation=ACTIVATION, padding='same')(premut_layer)
    conv_layer_t = layers.Permute((2,1,3))(conv_layer)
    loop_layer = layers.Concatenate()([conv_layer, conv_layer_t])

    # second res block: dilated_2d_conv -> relu -> dilated_2d_conv -> batch_normalization -> add -> relu
    for block in range(DIALETED_RESNET_BLOCKS):
        for loop in DIALETION:  # dilated ResNet
            loop_conv_layer = layers.Conv2D(64, DIALETED_RESNET_SIZE, activation=ACTIVATION, dilation_rate=loop, padding="same")(loop_layer)
            loop_conv_layer = layers.Conv2D(64, DIALETED_RESNET_SIZE, dilation_rate=loop , padding="same")(loop_conv_layer)  # TODO: add activation?
            batch_layer = layers.BatchNormalization()(loop_conv_layer)
            loop_layer = layers.Add()([batch_layer, loop_layer])
            loop_layer = layers.Activation(ACTIVATION)(loop_layer)

    dropout_layer = layers.Dropout(DROPOUT)(loop_layer)

    # distance
    distances = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    distances = layers.Conv2D(1, END_CONV_KER, activation="relu", padding="same")(distances)
    distance_t = layers.Permute((2,1,3))(distances)
    distances = layers.Add(name="distances")([0.5 * distances, 0.5 * distance_t])  # for symmetry

    #  omega
    omegas = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    omegas = layers.Conv2D(2, END_CONV_KER, activation="tanh", padding="same")(omegas)
    omega_t = layers.Permute((2,1,3))(omegas)
    omegas = layers.Add(name="omegas")([0.5 * omegas, 0.5 * omega_t])  # for symmetry

    # theta
    thetas = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    thetas = layers.Conv2D(2, END_CONV_KER, activation="tanh", padding="same", name="thetas")(thetas)

    # phi
    phis = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    phis = layers.Conv2D(2, END_CONV_KER, activation="tanh",padding="same", name="phis")(phis)

    return tf.keras.Model(input_layer, [distances, omegas, thetas, phis], name="NanoNet")


def plot_loss(history):

    ig, axes = plt.subplots(1, 4, figsize=(15,3))
    axes[0].plot(history.history['distances_loss'], label='Training loss')
    axes[0].plot(history.history['val_distances_loss'], label='Validation loss')
    axes[0].legend()
    corr = np.correlate(history.history['distances_loss'], history.history['val_distances_loss']).item()
    axes[0].set_title("Distance loss, correlation: %.1f"%corr)

    axes[1].plot(history.history['omegas_loss'], label='Training loss')
    axes[1].plot(history.history['val_omegas_loss'], label='Validation loss')
    axes[1].legend()
    corr = np.correlate(history.history['omegas_loss'], history.history['val_omegas_loss']).item()
    axes[1].set_title("Omega loss, correlation: %.1f"%corr)

    axes[2].plot(history.history['thetas_loss'], label='Training loss')
    axes[2].plot(history.history['val_thetas_loss'], label='Validation loss')
    axes[2].legend()
    corr = np.correlate(history.history['thetas_loss'], history.history['val_thetas_loss']).item()
    axes[2].set_title("Theta loss, correlation: %.3f"%corr)

    axes[3].plot(history.history['phis_loss'], label='Training loss')
    axes[3].plot(history.history['val_phis_loss'], label='Validation loss')
    axes[3].legend()
    corr = np.correlate(history.history['phis_loss'], history.history['val_phis_loss']).item()
    _=axes[3].set_title("Phi loss, correlation: %.3f"%corr)
    plt.savefig("lose_history_" + files_name)


def plot_test(trained_model, test_X, test_Y):

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


def print_parameters():


    print("Name: " + files_name)
    print("parameters summery:")
    print("DIM:" + str(DIM))
    print("CDR:" + str(CDR))
    print("RESNET_BLOCKS:" + str(RESNET_BLOCKS))
    print("RESNET_SIZE:" + str(RESNET_SIZE))
    print("FIRST_RESNET_SIZE: " + str(FIRST_RESNET_SIZE))
    print("DIALETED_RESNET_BLOCKS:" + str(DIALETED_RESNET_BLOCKS))
    print("DIALETED_RESNET_SIZE: " + str(DIALETED_RESNET_SIZE))
    print("DIALETION: " + str(DIALETION))
    print("EPOCHS:" + str(EPOCHS))
    print("LR:" + str(LR))
    print("TEST_SIZE:" + str(TEST_SIZE))
    print("BATCH:" + str(BATCH))
    print("BINS:" + str(BINS))
    print("POOL:" + str(POOL))
    print("ACTIVATION: " + str(ACTIVATION))
    print("END_ACTIVATION: " + str(END_ACTIVATION))
    print("DROPOUT: " + str(DROPOUT))
    print("LOSS: " + LOSS)
    print("END_CONV_SIZE: " + str(END_CONV_SIZE))
    print("END_CONV_KER: " + str(END_CONV_KER))
    print("DIALETED_RESNET_KERNELS: " + str(DIALETED_RESNET_KERNELS))


def angle(angle1, angle2):
    a = angle1 - angle2
    a = (a + np.pi) % (2*np.pi) - np.pi
    print(abs(a))


def calc_angle_std(trained_model, x, y):
    y_pred = trained_model.predict(x)

    cor = 0.000001
    omega_diff = np.arctan2(y[1][:,:,:,1] + cor, y[1][:,:,:,0]+ cor) - np.arctan2(y_pred[1][:,:,:,1]+ cor, y_pred[1][:,:,:,0]+ cor)
    theta_diff = np.arctan2(y[2][:,:,:,1]+ cor, y[2][:,:,:,0]+ cor) - np.arctan2(y_pred[2][:,:,:,1]+ cor, y_pred[2][:,:,:,0]+ cor)
    phi_diff = np.arctan2(y[3][:,:,:,1]+ cor, y[3][:,:,:,0]+ cor) - np.arctan2(y_pred[3][:,:,:,1]+ cor, y_pred[3][:,:,:,0]+ cor)

    omega_diff = ((omega_diff + np.pi) % (2*np.pi) - np.pi)**2
    theta_diff = ((theta_diff + np.pi) % (2*np.pi) - np.pi)**2
    phi_diff = ((phi_diff + np.pi) % (2*np.pi) - np.pi)**2

    print("Omega atan2 STD: {}".format(np.mean(omega_diff) ** 0.5))
    print("Theta atan2 STD: {}".format(np.mean(theta_diff) ** 0.5))
    print("Phi atan2 STD: {}".format(np.mean(phi_diff) ** 0.5))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("X", help="pickle file path")
    parser.add_argument("Y", help="pickle file path")
    parser.add_argument("pdb_names", help="the names of the pdbs")
    args = parser.parse_args()

    model = d2_net_architecture()

    model.summary()
    print_parameters()

    with open(args.X, "rb") as input_file:
        X = pickle.load(input_file)
    with open(args.Y, "rb") as feature_file:
        Y = pickle.load(feature_file)
    with open(args.pdb_names, "rb") as names_file:
        pdb_names = pickle.load(names_file)

    if DIM == 1:
        X = X[:, :, :, CDR_DICT[CDR]].reshape(-1, 32, 21, 1)

    train_index, test_index, _, _ = train_test_split(np.arange(len(X)), np.arange(len(Y)), test_size=TEST_SIZE)
    X_train, X_test, Y_train, Y_test = np.array(X[train_index,:,:,:]), np.array(X[test_index,:,:,:]), np.array(Y[train_index,:,:,:,:]), np.array(Y[test_index,:,:,:,:])
    test_names = pdb_names[test_index]
    train_names = pdb_names[train_index]
    # X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=TEST_SIZE)

    pickle.dump(X_test, open("test_X_" + files_name + ".pkl", "wb"))
    pickle.dump(Y_test, open("test_Y_" + files_name + ".pkl", "wb"))
    pickle.dump(X_train, open("train_X_" + files_name + ".pkl", "wb"))
    pickle.dump(Y_train, open("train_Y_" + files_name + ".pkl", "wb"))

    pickle.dump(test_names, open("test_names_" + files_name + ".pkl", "wb"))
    pickle.dump(train_names, open("train_names_" + files_name + ".pkl", "wb"))

    save_model = tf.keras.callbacks.ModelCheckpoint(filepath=files_name,save_best_only=True, verbose=1)

    #  all features
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=LR), loss=LOSS)
    net_history = model.fit(X_train, reshape_y(Y_train), validation_split=VAL_SIZE, epochs=EPOCHS, verbose=1, batch_size=BATCH, callbacks=[save_model])

    loss = model.evaluate(X_test, reshape_y(Y_test))
    print("##################################################################")
    print("test loss: {}".format(loss))
    print("##################################################################")
    model = tf.keras.models.load_model(files_name)
    loss = model.evaluate(X_test, reshape_y(Y_test))
    print("##################################################################")
    print("test loss: {}".format(loss))
    print("##################################################################")

    print_parameters()

    plot_loss(net_history)
    plot_test(model, X_test, Y_test)

    calc_angle_std(model, X_test, reshape_y(Y_test))
