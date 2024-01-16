def slidingRP_all(spikeTimes, spikeClusters, params={}):
    """
    Compute the metric for each cluster in a recording.

    Inputs:
    - spikeTimes: a vector of times in seconds
    - spikeClusters: cluster labels for each spike time
    - params: a dictionary which may contain various parameters

    Outputs:
    - rpMetrics: a list of dictionaries, each containing metrics for a cluster
    - cont: contamination levels
    - rp: refractory periods
    """

    # Default parameter values
    returnMatrix = params.get('returnMatrix', False)
    two_ms_NoSpikesCondition = params.get('2msNoSpikesCondition', False)
    FRthresh = params.get('FRthresh', 0.5) if two_ms_NoSpikesCondition else None
    verbose = params.get('verbose', True)

    # Unique cluster IDs
    cids = list(set(spikeClusters))

    rpMetrics = []

    if verbose:
        print(f'Computing metrics for {len(cids)} clusters')

    for cidx in cids:
        st = [spikeTimes[i] for i, cluster in enumerate(spikeClusters) if cluster == cidx]

        # Assuming slidingRP is a function that computes the required metrics
        metrics = slidingRP(st, params)

        cluster_metrics = {
            'cid': cidx,
            'maxConfidenceAt10Cont': metrics[0],
            'minContWith90Confidence': metrics[1],
            'timeOfLowestCont': metrics[2],
            'nSpikesBelow2': metrics[3],
        }

        # Add pass/fail value of metric
        if two_ms_NoSpikesCondition:
            if metrics[-1] < FRthresh:
                cluster_metrics['value'] = 0
            else:
                if metrics[3] == 0:
                    cluster_metrics['value'] = 1
                else:
                    cluster_metrics['value'] = 1 if metrics[1] <= 10 else 0
        else:
            cluster_metrics['value'] = 1 if metrics[1] <= 10 else 0

        if returnMatrix:
            cluster_metrics['confMatrix'] = metrics[4]

        if verbose:
            pfstring = 'PASS' if metrics[1] <= 10 else 'FAIL'
            print(f"  {cidx}: contamination = {metrics[1]:.1f}%, {pfstring}, max conf = {metrics[0]:.2f}%, time = {metrics[2]*1000:.2f} ms, n below 2 ms = {metrics[3]}")

        rpMetrics.append(cluster_metrics)

    # Assuming these are part of the returned metrics from slidingRP
    cont = metrics[5]
    rp = metrics[6]

    return rpMetrics, cont, rp
