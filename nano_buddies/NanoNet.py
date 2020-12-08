
import tensorflow as tf
from tensorflow.keras import layers
import pickle
import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.model_selection import train_test_split
import seaborn as sns

DIM = 2
DIALETED_RESNET_BLOCKS = 10
VARIANT = 1
EPOCHS = 100
LR = 0.01
TEST_SIZE = 0.023


class PolynomialDecay:
    def __init__(self, maxEpochs=EPOCHS, initAlpha=LR, power=1):
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

    return [y[:,0,:,:,0], y[:,1,:,:,:], y[:,2,:,:,:], y[:,3,:,:,:]]


def d1_net_architecture():
    """

    :return:
    """
    input = tf.keras.Input(shape=(32,21), name="loop")
    x = layers.Conv1D(32, 5, activation="relu", padding="same")(input)
    x = layers.MaxPooling1D(2)(x)
    x = layers.Conv1D(48, 5, activation="relu", padding="same")(x)
    x = layers.MaxPooling1D(2)(x)
    x = layers.Conv1D(64, 5, activation="relu", padding="same")(x)
    x = layers.UpSampling1D(size=2)(x)
    x = layers.Conv1D(48, 5, activation="relu", padding="same")(x)
    x = layers.UpSampling1D(size=2)(x)
    x = layers.Conv1D(32, 5, activation="relu", padding="same")(x)
    x = layers.Reshape((1,32,32))(x)

    for i in range(5):
        x = layers.Concatenate(axis=1)([x,x])
    x_T = layers.Permute((2,1,3))(x)
    x_2d = layers.Concatenate()([x, x_T])


    for i in range(5):
        x_2d = layers.Conv2D(64, 5, activation="relu", dilation_rate=2**i, padding="same")(x_2d)
        x_2d = layers.Conv2D(64, 5, activation="linear", dilation_rate=2**i , padding="same")(x_2d)
        x_2d = layers.BatchNormalization()(x_2d)

    x_3 = layers.Dropout(.2)(x_2d)

    distances = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    distance_T = layers.Permute((2,1,3))(distances)
    distances = layers.Add()([distances, distance_T])
    distances = layers.Reshape((32,32), name="distances")(distances)

    omegas = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    omega_T = layers.Permute((2,1,3))(omegas)
    omegas = layers.Add()([omegas, omega_T])
    omegas = layers.Concatenate(name="omegas")([tf.math.cos(omegas), tf.math.sin(omegas)])

    thetas = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    thetas = layers.Concatenate(name="thetas")([tf.math.cos(thetas), tf.math.sin(thetas)])

    phis = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    phis = layers.Concatenate(name="phis")([tf.math.cos(phis), tf.math.sin(phis)])

    return tf.keras.Model(input, [distances, omegas, thetas, phis], name="NanoNet1d")


def d2_net_architecture(variant=2):
    """

    :return:
    """

    input_layer = tf.keras.Input(shape=(32, 21, 3), name="InputLayer")
#  first ResNet block (or conv->relu->pool)

    #  variant 1 (ResNet blocks)
    if variant == 1:
        loop_layer = layers.Conv2D(32, (15,15), padding='same', name='Conv2D_1')(input_layer)
        for i in range(4):
            conv_layer = layers.Conv2D(32, (15,15), activation='relu', padding='same')(loop_layer)
            conv_layer = layers.Conv2D(32, (15,15), padding='same')(conv_layer)
            loop_layer = layers.BatchNormalization()(conv_layer)

        relu_layer = layers.ReLU()(loop_layer)
        premut_layer = layers.Permute((1,3,2))(relu_layer)

    #  variant 2 (conv->relu->pool)
    elif variant == 2:

        conv_layer = layers.Conv2D(32, (5,5), activation='relu', padding='same', name='Conv2D_1')(input_layer)
        max_pool_layer = layers.MaxPooling2D((2,2))(conv_layer)

        conv_layer = layers.Conv2D(48, (5,5), activation='relu', padding='same', name='Conv2D_2')(max_pool_layer)
        max_pool_layer = layers.MaxPooling2D((2,2))(conv_layer)

        conv_layer = layers.Conv2D(64, (5,5), activation='relu', padding='same', name='Conv2D_3')(max_pool_layer)
        up_sampling_layer = layers.UpSampling2D((2,2))(conv_layer)

        conv_layer = layers.Conv2D(48, (5,5), activation='relu', padding='same', name='Conv2D_4')(up_sampling_layer)
        up_sampling_layer = layers.UpSampling2D((2,2))(conv_layer)

        conv_layer = layers.Conv2D(32, (5,5), activation='relu', padding='same', name='Conv2D_5')(up_sampling_layer)
        premut_layer = layers.Permute((1,3,2))(conv_layer)

    else:
        raise ValueError

    conv_layer = layers.Conv2D(32, (5,5), activation='relu', padding='same', name='Conv2D_6')(premut_layer)
    conv_layer_t = layers.Permute((2,1,3))(conv_layer)
    loop_layer = layers.Concatenate()([conv_layer, conv_layer_t])

    for block in range(DIALETED_RESNET_BLOCKS):
        for loop in range(5):  # dilated ResNet
            loop_conv_layer = layers.Conv2D(64, (5,5), activation="relu", dilation_rate=2**loop, padding="same")(loop_layer)
            loop_conv_layer = layers.Conv2D(64, (5,5), dilation_rate=2**loop , padding="same")(loop_conv_layer)  # TODO: add activation?
            loop_layer = layers.BatchNormalization()(loop_conv_layer)

    relu_layer = layers.ReLU()(loop_layer)
    dropout_layer = layers.Dropout(0.2)(relu_layer)


    #  distance
    distances = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    for loop in [4,1]:
        distances = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(distances)

    distance_t = layers.Permute((2,1,3))(distances)
    distances = layers.Add()([distances, distance_t])  # for symmetry
    distances = layers.Reshape((32,32), name="distances")(distances)  # for 1D again

    # #  omega
    # omegas = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    # for loop in [4,1]:
    #     omegas = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(omegas)
    #
    # omega_t = layers.Permute((2,1,3))(omegas)
    # omegas = layers.Add()([omegas, omega_t])  # for symmetry
    # omegas = layers.Concatenate(name="omegas")([tf.math.cos(omegas), tf.math.sin(omegas)])
    #
    # # theta
    # thetas = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    # for loop in [4,1]:
    #     thetas = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(thetas)
    # thetas = layers.Concatenate(name="thetas")([tf.math.cos(thetas), tf.math.sin(thetas)])
    #
    # # phi
    # phis = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    # for loop in [4,1]:
    #     phis = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(phis)
    # phis = layers.Concatenate(name="phis")([tf.math.cos(phis), tf.math.sin(phis)])

    return tf.keras.Model(input_layer, distances, name="NanoNet2d")


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
    plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("X_train", help="pickle file path")
    parser.add_argument("Y_train", help="pickle file path")
    args = parser.parse_args()

    if DIM == 2:
        model = d2_net_architecture(VARIANT)
    else:
        model = d1_net_architecture()

    with open(args.X_train, "rb") as input_file:
        X = pickle.load(input_file)
    with open(args.Y_train, "rb") as feature_file:
        Y = pickle.load(feature_file)

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=TEST_SIZE)

    pickle.dump(X_test, open("test_X.pkl", "wb"))
    pickle.dump(Y_test, open("test_Y.pkl", "wb"))

    Y_train, Y_test = reshape_y(Y_train), reshape_y(Y_test)

    poly_decay = tf.keras.callbacks.LearningRateScheduler(PolynomialDecay())
    step_decay = tf.keras.callbacks.LearningRateScheduler(StepDecay())
    platau_decay = tf.keras.callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.25, patience=5, verbose=1, min_delta=0.0001)

    #  only distance
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.01, decay=0.0), loss="mse")
    net_history = model.fit(X_train, Y_train[0], validation_split=0.05, epochs=50, verbose=1, batch_size=32, callbacks=[poly_decay])


    #  all features
    # model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.01, decay=0.98), loss=['mse', 'mse','mse','mse'], loss_weights=[0.05,1,1,1])  # TODO: use huber loss on angles?
    # net_history = model.fit(X_train, Y_train, validation_split=0.05, epochs=50, verbose=1, batch_size=32)
    model.save("model2")
    # plot_loss(net_history)

    loss = model.evaluate(X_test, Y_test[0])
    print("test loss: {}".format(loss))

    # model = tf.keras.models.load_model("model1")

    # r = sns.heatmap(Y_test[0][0,:,:], cmap='BuPu')
    # plt.show()
    #
    # z = model.predict(X_test[0:2,:,:,:])[0][0]
    # r2 = sns.heatmap(z, cmap='BuPu')
    # plt.show()
    #
    # r3 = sns.heatmap(tf.math.atan2(x=Y_test[3][0,:,:,0], y=Y_test[3][0,:,:,1]), cmap='BuPu')
    # plt.show()
    #
    # z2 = model.predict(X_test[0:2,:,:,:])[3][0]
    # r5 = sns.heatmap(tf.math.atan2(x=z2[:,:,0],y=z2[:,:,1]), cmap='BuPu')
    #
    # plt.show()
    #
    #
    # r6 = sns.heatmap(Y_test[0][1,:,:], cmap='BuPu')
    # plt.show()
    #
    # z = model.predict(X_test[0:2,:,:,:])[0][1]
    # r9 = sns.heatmap(z, cmap='BuPu')
    # plt.show()
    #
    # r89 = sns.heatmap(tf.math.atan2(x=Y_test[3][1,:,:,0], y=Y_test[3][1,:,:,1]), cmap='BuPu')
    # plt.show()
    #
    # z2 = model.predict(X_test[0:2,:,:,:])[3][1]
    # r99 = sns.heatmap(tf.math.atan2(x=z2[:,:,0],y=z2[:,:,1]), cmap='BuPu')
    #
    # plt.show()
