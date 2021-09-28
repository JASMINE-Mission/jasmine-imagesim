import numpy as np

def gaussian_flat(Nside=1024,sigma=0.01, seed=1):
    """
    Summary:
        This function makes an interpixel flat
        consists of gaussian noise.

    Args:
        Nside (int)  : Number of pixels on a side.
                       (Default: 1024)
        sigma (float): Sigma of the gaussian noise.
                       (Default: 0.01)
        seed (int) : random seed

    Returns:
        flat (ndarray): Calculated gaussian noise
                        interpixel flat (Nside x Nside).

    """

    np.random.seed(seed)
    flat = np.random.normal(1.0, sigma, Nside*Nside)
    flat = flat.reshape((Nside,Nside))
    np.random.seed() #reset random seed
    
    return flat


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    flat = gaussian_flat()
    a=plt.imshow(flat)
    plt.colorbar(a)
    plt.savefig("flat.png")
