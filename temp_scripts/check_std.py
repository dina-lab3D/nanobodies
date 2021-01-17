import numpy as np

OMEGA_VAR = 0.036
THETA_VAR = 0.0153
PHI_VAR = 0.0082


if __name__ == '__main__':
    angles_vars = [OMEGA_VAR, THETA_VAR, PHI_VAR]
    angles = ["Omega", "Theta", "Phi"]

    for angle, angle_var in zip(angles, angles_vars):

        stds = []
        sanity_stds = []

        for i in range(1000):
            std = angle_var**0.5
            sin_true = np.random.uniform(-1,1.0000001, 161)
            cos_true = np.random.uniform(-1,1.0000001, 161)
            real_tan = np.arctan2(sin_true, cos_true)

            sin = np.zeros(161)
            cos = np.zeros(161)
            for i in range(161):
                sin[i] = np.random.normal(sin_true[i], std, 1)
                cos[i] = np.random.normal(cos_true[i], std, 1)
            tan = np.arctan2(sin, cos)
            sanity_stds.append((np.mean((sin_true - sin)**2) + np.mean((cos_true-cos)**2))/2)
            diff = ((tan - real_tan) + np.pi) % (2 * np.pi) - np.pi
            tan_std = (np.sum(diff**2) / 160) ** 0.5
            stds.append(tan_std)

        mean_std = np.mean(stds)
        sanity_std = np.mean(sanity_stds)

        print("{} atan std: {}".format(angle, mean_std))
        print("{} sanity std: {}".format(angle, sanity_std))

