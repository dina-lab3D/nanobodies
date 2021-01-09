
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

DIM = 2
RESNET_BLOCKS = 3
RESNET_SIZE = (17, 17)
FIRST_RESNET_SIZE = (17, 17)
DIALETED_RESNET_BLOCKS = 5
DIALETION = [1,2,4,8,16]
DIALETED_RESNET_SIZE = (5, 5)
VARIANT = 1
EPOCHS = 200
LR = 0.001
TEST_SIZE = 0.05
BATCH = 32
DROPOUT = 0.2
END_CONV_SIZE = 4   # normal is 4
END_CONV_KER = (5, 5)  # normal is 5
DIALETED_RESNET_KERNELS = 64
ACTIVATION = "relu"
END_ACTIVATION = "elu"
LOSS = "mse"
BINS = False
POOL = False
files_name = "save_9"


class PolynomialDecay:
    def __init__(self, maxEpochs=EPOCHS, initAlpha=LR, power=0.9):
        # store the maximum number of epochs, base learning rate,
        # and power of the polynomial
        self.maxEpochs = maxEpochs
        self.initAlpha = initAlpha
        self.power = power

    def __call__(self, epoch):
        # compute the new learning rate based on polynomial decay
        decay = (1 - (epoch / float(self.maxEpochs))) ** self.power
        alpha = self.initAlpha * decay
        # return the new learning rate
        return float(alpha)


class StepDecay:
    def __init__(self, initAlpha=LR, factor=0.75, dropEvery=5):
        # store the base initial learning rate, drop factor, and
        # epochs to drop every
        self.initAlpha = initAlpha
        self.factor = factor
        self.dropEvery = dropEvery

    def __call__(self, epoch):
        # compute the learning rate for the current epoch
        exp = np.floor((1 + epoch) / self.dropEvery)
        alpha = self.initAlpha * (self.factor ** exp)
        # return the learning rate
        return float(alpha)


def reshape_y(y):
    if BINS:
        return [y[:,0,:,:], y[:,1,:,:], y[:,2,:,:], y[:,3,:,:]]
    return [y[:,0,:,:,0].reshape(-1,32,32,1), y[:,1,:,:,:], y[:,2,:,:,:], y[:,3,:,:,:]]


def d2_net_architecture(variant=2):
    """

    :return:
    """
    if DIM == 1:
        input_layer = tf.keras.Input(shape=(32, 21, 1), name="InputLayer")
    else:
        input_layer = tf.keras.Input(shape=(32, 21, 3), name="InputLayer")
    #  first ResNet block (or conv->relu->pool)

    #  variant 1 (ResNet blocks)
    loop_layer = layers.Conv2D(32, FIRST_RESNET_SIZE, padding='same')(input_layer)

    for i in range(RESNET_BLOCKS):
        # for loop in [1,2,4,8,16]:
            conv_layer = layers.Conv2D(32, RESNET_SIZE, activation=ACTIVATION, padding='same')(loop_layer)
            conv_layer = layers.Conv2D(32, RESNET_SIZE, padding='same')(conv_layer)
            batch_layer = layers.BatchNormalization()(conv_layer)
            loop_layer = layers.Add()([batch_layer, loop_layer])
            loop_layer = layers.Activation(ACTIVATION)(loop_layer)

    # loop_layer = layers.Conv2D(3, FIRST_RESNET_SIZE, padding='same')(loop_layer)
    # loop_layer = layers.BatchNormalization()(loop_layer)
    # loop_layer = layers.Add()([input_layer, loop_layer])
    # loop_layer = layers.Conv2D(32, FIRST_RESNET_SIZE, padding='same', activation='relu')(loop_layer)

    premut_layer = layers.Permute((1,3,2))(loop_layer)

    # for kernel 64
    conv_layer = layers.Conv2D(32, (5,5), activation=ACTIVATION, padding='same', name='Conv2D_6')(premut_layer)
    conv_layer_t = layers.Permute((2,1,3))(conv_layer)
    loop_layer = layers.Concatenate()([conv_layer, conv_layer_t])

    for block in range(DIALETED_RESNET_BLOCKS):
        for loop in DIALETION:  # dilated ResNet
            loop_conv_layer = layers.Conv2D(64, DIALETED_RESNET_SIZE, activation=ACTIVATION, dilation_rate=loop, padding="same")(loop_layer)
            loop_conv_layer = layers.Conv2D(64, DIALETED_RESNET_SIZE, dilation_rate=loop , padding="same")(loop_conv_layer)  # TODO: add activation?
            batch_layer = layers.BatchNormalization()(loop_conv_layer)
            loop_layer = layers.Add()([batch_layer, loop_layer])
            loop_layer = layers.Activation(ACTIVATION)(loop_layer)

    # loop_layer = layers.Conv2D(21, DIALETED_RESNET_SIZE, padding='same')(loop_layer)
    # loop_layer = layers.BatchNormalization()(loop_layer)
    # loop_layer = layers.Add()([premut_layer, loop_layer])

    dropout_layer = layers.Dropout(DROPOUT)(loop_layer)

    # distances = layers.Conv2D(1, END_CONV_KER, padding="same")(dropout_layer)
    # distances_res = layers.Conv2D(1, END_CONV_KER, activation=END_ACTIVATION, padding="same")(distances)
    # distances_res = layers.Conv2D(1, END_CONV_KER , padding="same")(distances_res)
    # distances_res = layers.BatchNormalization()(distances_res)
    # distances = layers.Add()([distances_res, distances])
    # distances = layers.Activation("relu")(distances)

    # distance
    distances = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    # distances = layers.PReLU()(distances)

    dims = [1]
    if BINS:
        dims = [41]
    for loop in dims:
        # distances = layers.Conv2D(4, (5, 5), activation="relu",padding="same")(distances)
        distances = layers.Conv2D(loop, END_CONV_KER, activation="relu", padding="same")(distances)
        distance_t = layers.Permute((2,1,3))(distances)
        distances = layers.Add(name="distances")([0.5 * distances, 0.5 * distance_t])  # for symmetry

        # distances = layers.Activation(name="distances", activation="relu")(distances)

    # distance_t = layers.Permute((2,1,3))(distances)
    # if BINS:
    #     distances = layers.Add()([0.5 * distances, 0.5 * distance_t])  # for symmetry
    #     distances = layers.Softmax(name="distances")(distances)
    # else:
    #     distances = layers.Add(name="distances")([0.5 * distances, 0.5 * distance_t])  # for symmetry
    #     # distances = layers.Reshape((32,32), name="distances")(distances)  # for 1D again


    # omegas = layers.Conv2D(2, END_CONV_KER, padding="same")(dropout_layer)
    # omegas_res = layers.Conv2D(2, END_CONV_KER, activation=END_ACTIVATION, padding="same")(omegas)
    # omegas_res = layers.Conv2D(2, END_CONV_KER , padding="same")(omegas_res)
    # omegas_res = layers.BatchNormalization()(omegas_res)
    # omegas = layers.Add()([omegas_res, omegas])
    # omegas = layers.Activation("tanh")(omegas)

    #  omega
    omegas = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    # omegas = layers.PReLU()(omegas)

    dims = [2]
    if BINS:
        dims = [24]
    for loop in dims:
        # omegas = layers.Conv2D(4, (5,5), activation="relu", padding="same")(omegas)
        omegas = layers.Conv2D(loop, END_CONV_KER, activation="tanh", padding="same")(omegas)
        omega_t = layers.Permute((2,1,3))(omegas)
        omegas = layers.Add(name="omegas")([0.5 * omegas, 0.5 * omega_t])  # for symmetry

        # omegas = layers.Activation(name="omegas", activation="tanh")(omegas)


    # omega_t = layers.Permute((2,1,3))(omegas)
    # if BINS:
    #     omegas = layers.Add()([0.5 * omegas, 0.5 * omega_t])  # for symmetry
    #     omegas = layers.Softmax(name="omegas")(omegas)
    # else:
    #     omegas = layers.Add(name="omegas")([0.5 * omegas, 0.5 * omega_t])  # for symmetry


    # thetas = layers.Conv2D(2, END_CONV_KER, padding="same")(dropout_layer)
    # thetas_res = layers.Conv2D(2, END_CONV_KER, activation=END_ACTIVATION, padding="same")(thetas)
    # thetas_res = layers.Conv2D(2, END_CONV_KER , padding="same")(thetas_res)
    # thetas_res = layers.BatchNormalization()(thetas_res)
    # thetas = layers.Add()([thetas_res, thetas])
    # thetas = layers.Activation("tanh")(thetas)


    # theta
    thetas = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    # thetas = layers.PReLU()(thetas)

    dims = [2]
    if BINS:
        dims = [24]
    for loop in dims:
        if loop == 2 and not BINS:
            # thetas = layers.Conv2D(4, (5, 5), activation="relu", padding="same")(thetas)
            thetas = layers.Conv2D(loop, END_CONV_KER, activation="tanh", padding="same", name="thetas")(thetas)
        else:
            thetas = layers.Conv2D(loop, END_CONV_KER, activation="tanh", padding="same")(thetas)
    # if BINS:
    #     thetas = layers.Softmax(name="thetas")(thetas)



    # phis = layers.Conv2D(2, END_CONV_KER, padding="same")(dropout_layer)
    # phis_res = layers.Conv2D(2, END_CONV_KER, activation=END_ACTIVATION, padding="same")(phis)
    # phis_res = layers.Conv2D(2, END_CONV_KER , padding="same")(phis_res)
    # phis_res = layers.BatchNormalization()(phis_res)
    # phis = layers.Add()([phis_res, phis])
    # phis = layers.Activation("tanh")(phis)


    # phi
    phis = layers.Conv2D(END_CONV_SIZE, END_CONV_KER, padding="same", activation=END_ACTIVATION)(dropout_layer)
    # phis = layers.PReLU()(phis)

    dims = [2]
    if BINS:
        dims = [12]
    for loop in dims:
        if loop == 2 and not BINS:
            # phis = layers.Conv2D(4, (5, 5), activation="relu",padding="same")(phis)
            phis = layers.Conv2D(loop, END_CONV_KER, activation="tanh",padding="same", name="phis")(phis)
        else:
            phis = layers.Conv2D(loop, END_CONV_KER, activation="tanh", padding="same")(phis)
    # if BINS:
    #     phis = layers.Softmax(name="phis")(phis)

    return tf.keras.Model(input_layer, [distances, omegas, thetas, phis], name="NanoNet2d")


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
                sns.heatmap(test_Y[i][j,:,:,1], ax=axes[j][0],xticklabels=False, yticklabels=False, cbar=False)
                sns.heatmap(pred_Y[j][i,:,:,1], ax=axes[j][1],xticklabels=False, yticklabels=False, cbar=False)
                # sns.heatmap(tf.math.atan2(x=test_Y[i][j,:,:,0],y=test_Y[i][j,:,:,1]), ax=axes[j][0],xticklabels=False, yticklabels=False, cbar=False)
                # sns.heatmap(tf.math.atan2(x=pred_Y[j][i,:,:,0],y=pred_Y[j][i,:,:,1]), ax=axes[j][1],xticklabels=False, yticklabels=False, cbar=False)

            axes[j][0].set_title(names[j] + " true")
            axes[j][1].set_title(names[j] + " predicted")

        plt.savefig(os.path.join("test_vs_pred_" + files_name, "test_" + str(i)))
        plt.clf()


def print_parameters():


    print("Name: " + files_name)
    print("parameters summery:")
    print("DIM:" + str(DIM))
    print("RESNET_BLOCKS:" + str(RESNET_BLOCKS))
    print("RESNET_SIZE:" + str(RESNET_SIZE))
    print("FIRST_RESNET_SIZE: " + str(FIRST_RESNET_SIZE))
    print("DIALETED_RESNET_BLOCKS:" + str(DIALETED_RESNET_BLOCKS))
    print("DIALETED_RESNET_SIZE: " + str(DIALETED_RESNET_SIZE))
    print("DIALETION: " + str(DIALETION))
    print("VARIANT:" + str(VARIANT))
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("X", help="pickle file path")
    parser.add_argument("Y", help="pickle file path")
    parser.add_argument("pdb_names", help="the names of the pdbs")
    args = parser.parse_args()

    model = d2_net_architecture(VARIANT)

    model.summary()
    print_parameters()

    with open(args.X, "rb") as input_file:
        X = pickle.load(input_file)
    with open(args.Y, "rb") as feature_file:
        Y = pickle.load(feature_file)
    with open(args.pdb_names, "rb") as names_file:
        pdb_names = pickle.load(names_file)

    if DIM == 1:
        X = X[:,:,:,0].reshape(-1,32,21,1)

    train_index, test_index, _, _ = train_test_split(np.arange(len(X)), np.arange(len(Y)), test_size=TEST_SIZE)
    X_train, X_test, Y_train, Y_test = np.array(X[train_index,:,:,:]), np.array(X[test_index,:,:,:]), np.array(Y[train_index,:,:,:,:]), np.array(Y[test_index,:,:,:,:])
    test_names = pdb_names[test_index]

    # X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=TEST_SIZE)

    pickle.dump(X_test, open("test_X_" + files_name + ".pkl", "wb"))
    pickle.dump(Y_test, open("test_Y_" + files_name + ".pkl", "wb"))
    pickle.dump(test_names, open("test_names_" + files_name + ".pkl", "wb"))

    poly_decay = tf.keras.callbacks.LearningRateScheduler(PolynomialDecay())
    step_decay = tf.keras.callbacks.LearningRateScheduler(StepDecay())
    platau_decay = tf.keras.callbacks.ReduceLROnPlateau()
    save_model = tf.keras.callbacks.ModelCheckpoint(filepath=files_name,save_best_only=True, verbose=1)
    #  all features
    if BINS:
        loss = "sparse_categorical_crossentropy"
    else:
        loss = LOSS
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=LR), loss=loss)  # TODO: use huber loss on angles?
    net_history = model.fit(X_train, reshape_y(Y_train), validation_data=(X_test, reshape_y(Y_test)), epochs=EPOCHS, verbose=1, batch_size=BATCH, callbacks=[save_model])
    # model.save("model" + "_" + str(VARIANT))
    # model = tf.keras.models.load_model("model1")

    loss = model.evaluate(X_test, reshape_y(Y_test))
    print("##################################################################")
    print("test loss: {}".format(loss))
    print("##################################################################")
    
    print_parameters()

    plot_loss(net_history)
    plot_test(model, X_test, Y_test)

