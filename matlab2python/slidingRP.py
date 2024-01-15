import numpy as np

def slidingRP(spikeTimes, **params):
    """
    Compute the metric for a single cluster in a recording.

    Inputs:
    - spikeTimes: a vector of times in seconds
    - params: a dictionary which may contain:
        - contaminationThresh: acceptable contamination level to pass the neuron, default 10 (%)
        - confidenceThresh: required confidence level to decide that a given contamination level is not breached, default 90 (%)
        - recDur: recording duration in seconds
        - cont: a vector of contamination levels at which to test the neuron

    Outputs:
    - passTest: binary indicating whether the neuron's contamination level is below contaminationThresh with at least confidenceThresh level of confidence
    - confidence: the confidence that you have less than the threshold level of contamination
    - contamination: the minimum contamination for which you have at least the threshold level of confidence
    - timeOfLowestCont: the time at which best score happens
    - nACGBelow2: the number of spikes in the ACG below 2 ms
    - confMatrix, cont, rp, nACG: outputs from computeMatrix
    """

    # Default parameter values
    contThresh = params.get('contaminationThresh', 10)
    confThresh = params.get('confidenceThresh', 90)

    # Assuming computeMatrix is a function that computes the required matrix
    confMatrix, cont, rp, nACG = computeMatrix(spikeTimes, params)

    # Analysis conditions
    testTimes = rp > 0.0005

    # Compute confidence
    confidence = np.max(confMatrix[cont >= contThresh, :][:, testTimes])

    # Find indices where confidence threshold is exceeded
    ii, jj = np.where(confMatrix[:, testTimes] > confThresh)
    if ii.size > 0:
        minI = np.min(ii)
        contamination = cont[minI]
        minRP = np.argmax(confMatrix[minI, testTimes]) + np.where(testTimes)[0][0]
        timeOfLowestCont = rp[minRP]
    else:
        contamination = np.nan
        timeOfLowestCont = np.nan

    # Compute number of spikes in the ACG below 2 ms
    nACGBelow2 = np.sum(nACG[:np.where(rp > 0.002)[0][0]])

    # Determine pass test
    passTest = confidence > confThresh

    return passTest, confidence, contamination, timeOfLowestCont, nACGBelow2, confMatrix, cont, rp, nACG
