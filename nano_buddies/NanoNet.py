
import tensorflow as tf
from tensorflow.keras import layers, models



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
    distances = layers.Reshape((32,32))(distances)

    omegas = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    omega_T = layers.Permute((2,1,3))(omegas)
    omegas = layers.Add()([omegas, omega_T])
    omegas = layers.Reshape((32,32))(omegas)

    thetas = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    thetas = layers.Reshape((32,32))(thetas)

    phis = layers.Conv2D(1, 17, activation="linear", padding="same")(x_3)
    phis = layers.Reshape((32,32))(phis)
    tf.keras.Model(input, [distances, omegas, thetas, phis], name="halufit").summary()


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
    distances = layers.Reshape((32,32))(distances)  # for 1D again

    #  omega
    omegas = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    for loop in [4,1]:
        omegas = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(omegas)

    omega_t = layers.Permute((2,1,3))(omegas)
    omegas = layers.Add()([omegas, omega_t])  # for symmetry
    omegas = layers.Reshape((32,32))(omegas)  # for 1D again

    # theta
    thetas = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    for loop in [4,1]:
        thetas = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(thetas)
    thetas = layers.Reshape((32,32))(thetas)

    # phi
    phis = layers.Conv2D(16, (5,5), activation="relu", padding="same")(dropout_layer)
    for loop in [4,1]:
        phis = layers.Conv2D(loop, (5,5), activation="relu", padding="same")(phis)
    phis = layers.Reshape((32,32))(phis)

    tf.keras.Model(input_layer, [distances, omegas, thetas, phis], name="halufit").summary()





d1_net_architecture()
d2_net_architecture(1)