from scipy import stats
import matplotlib.pyplot as plt
import mplhep

def z_for_probability(probability, two_sided=True):
    """
    Calculate the number of standard deviations for a given probability.
    
    Parameters:
    -----------
    probability : float
        The probability (between 0 and 1)
    two_sided : bool
        If True, calculates for a two-sided interval (default)
        If False, calculates for one-sided (upper tail)
    
    Returns:
    --------
    float : The number of standard deviations (z-score)
    """
    if two_sided:
        # For two-sided, we want P(-z < X < z) = probability
        # So we need the quantile at (1 + probability) / 2
        alpha = (1 - probability) / 2
        z = stats.norm.ppf(1 - alpha)
    else:
        # For one-sided, we want P(X < z) = probability
        z = stats.norm.ppf(probability)
    
    return z


# Examples
if __name__ == "__main__":
    # Common confidence levels (two-sided)
    print("Two-sided intervals:")
    print(f"68.27% (1σ): {z_for_probability(0.6827):.4f}")
    print(f"90%:         {z_for_probability(0.90):.4f}")
    print(f"95%:         {z_for_probability(0.95):.4f}")
    print(f"95.45% (2σ): {z_for_probability(0.9545):.4f}")
    print(f"99%:         {z_for_probability(0.99):.4f}")
    print(f"99.73% (3σ): {z_for_probability(0.9973):.4f}")
    print(f"99.99% (4σ): {z_for_probability(0.9999):.4f}")
    print(f"5σ:          {z_for_probability(0.9999994267):.4f}")
    
    print("\nOne-sided (upper tail):")
    print(f"95%:         {z_for_probability(0.95, two_sided=False):.4f}")
    print(f"99%:         {z_for_probability(0.99, two_sided=False):.4f}")