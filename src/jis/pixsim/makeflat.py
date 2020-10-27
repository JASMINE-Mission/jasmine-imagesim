import numpy as np

def gaussian_flat(Nside=1024,sigma=0.01):
    """
    Summary:
        This function makes an interpixel flat
        consists of gaussian noise.

    Args:
        Nside (int)  : Number of pixels on a side.
                       (Default: 1024)
        sigma (float): Sigma of the gaussian noise.
                       (Default: 0.01)

    Returns:
        flat (ndarray): Calculated gaussian noise
                        interpixel flat (Nside x Nside).

    """

    flat = np.random.normal(1.0, sigma, Nside*Nside)
    flat = flat.reshape((Nside,Nside))

    return flat


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    flat = gaussian_flat()
    print(np.mean(flat))
    print(np.std(flat))

    a=plt.imshow(flat)
    plt.colorbar(a)
    plt.savefig("flat.png")
