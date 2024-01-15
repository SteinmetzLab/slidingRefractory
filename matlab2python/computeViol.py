import numpy as np
from scipy.stats import poisson

def computeViol(obsViol, firingRate, spikeCount, refDur, contaminationProp, recDur):
    """
    Compute the confidence score and expected violations.

    Parameters:
    - obsViol: Observed violations
    - firingRate: Firing rate of the neuron
    - spikeCount: Total spike count
    - refDur: Refractory period duration in seconds
    - contaminationProp: Query contamination level as a proportion (e.g., 0.1 for 10% contamination)
    - recDur: Recording duration in seconds

    Returns:
    - confidenceScore: The confidence score
    - expectedViol: The expected number of violations
    """

    # Calculate Nc and Nb based on spikeCount and contaminationProp
    Nc = spikeCount * contaminationProp
    Nb = spikeCount * (1 - contaminationProp)

    # Calculate expected violations
    expectedViol = 2 * refDur / recDur * Nc * (Nb + (Nc - 1) / 2)

    # Calculate the confidence score
    confidenceScore = 1 - poisson.cdf(obsViol, expectedViol)

    return confidenceScore, expectedViol
