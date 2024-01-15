import numpy as np

def computeMatrix(spikeTimes, params):
    """
    Core function for Sliding RP metric.

    Inputs:
    - spikeTimes: a vector of times in seconds
    - params: a dictionary which may have:
       - recDur: recording duration in seconds, if not specified will be taken
         as the max spike time - highly recommended to specify
       - cont: a vector of contamination levels at which to test the neuron.

    Outputs:
    - confMatrix: result of the test, dimensions [contamination levels, RP duration]
    - cont: scalar or vector of the contamination levels tested
    - rp: vector of the RP durations tested
    - nACG: the ACG function of the neuron at the values in rp. Units are counts
    """

    # Set cont and recDur based on params
    cont = params.get('cont', np.arange(0.5, 35.5, 0.5))
    recDur = params.get('recDur', max(spikeTimes))

    # Define rpEdges
    rpBinSize = 1 / 30000
    rpEdges = np.arange(0, 10 / 1000 + rpBinSize, rpBinSize)

    # Compute firing rate and spike count
    spikeCount = len(spikeTimes)
    firingRate = []  # Firing rate no longer required with updated calculation

    nACG, rp = histdiff(spikeTimes, spikeTimes, rpEdges)

    confMatrix = np.full((len(cont), len(rp)), np.nan)
    for rpIdx in range(len(rp)):
        # Compute observed violations
        obsViol = np.sum(nACG[0, :rpIdx + 1])

        for cidx in range(len(cont)):
            confMatrix[cidx, rpIdx] = 100 * computeViol(
                obsViol, firingRate, spikeCount, rp[rpIdx] + rpBinSize / 2, cont[cidx] / 100, recDur)

    return confMatrix, cont, rp, nACG
