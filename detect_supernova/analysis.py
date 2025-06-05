from scipy.stats import sigmaclip

def estimate_background(array, low=5.0, high=5.0):
    array = [~np.isnan(array)]
    clipped, _, _ = sigmaclip(array, low=high, high=high)
    return clipped
