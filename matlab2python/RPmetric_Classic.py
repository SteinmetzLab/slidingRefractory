import numpy as np

def RPmetric_Classic(spikeTimes, params):
    """
    Compute the classic refractory period (RP) metric.

    Inputs:
    - spikeTimes: a vector of times in seconds
    - params: a dictionary containing various parameters

    Outputs:
    - passTest: Boolean indicating if the test is passed
    - estContam: Estimated contamination
    - rp: Refractory periods
    - nACG: Auto-correlogram
    """

    # Extract parameters from params, with defaults
    metricType = params.get('metricType', 'Llobet')
    contaminationThresh = params.get('contaminationThresh', 10)
    recDur = params.get('recDur', max(spikeTimes))
    RPdur = params.get('RPdur', 0.002)

    # Check if pre-calculated ACG is provided
    if 'nACG' in params and 'rp' in params and 'spikeCount' in params:
        nACG = params['nACG']
        rp = params['rp']
        spikeCount = params['spikeCount']

        rpIdx = np.where(rp > RPdur)[0][0]
        obsViol = np.sum(nACG[:rpIdx])
    else:
        # Calculate observed violations
        rpBinSize = 1 / 30000
        rpEdges = np.arange(0, RPdur + rpBinSize, rpBinSize)
        
        spikeCount = len(spikeTimes)
        
        # Assuming histdiff is a function that computes the histogram difference
        nACG, rp = histdiff(spikeTimes, spikeTimes, rpEdges)
        
        obsViol = np.sum(nACG)

    expectedNc = spikeCount * contaminationThresh / 100
    expectedNb = spikeCount * (1 - contaminationThresh / 100)

    # Compute expected violations and estimated contamination
    if metricType == 'Llobet':
        expectedViol = 2 * RPdur / recDur * expectedNc * (expectedNb + (expectedNc - 1) / 2)
        estContam = 1 - np.sqrt(1 - obsViol * recDur / (spikeCount ** 2 * RPdur))
    elif metricType == 'Hill':
        expectedViol = 2 * RPdur / recDur * expectedNc * (expectedNb + expectedNc)
        estContam = 1 / 2 * (1 - np.sqrt(1 - 2 * obsViol * recDur / (spikeCount ** 2 / RPdur)))

    passTest = obsViol <= expectedViol

    return passTest, estContam, rp, nACG

