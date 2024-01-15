import matplotlib.pyplot as plt
import numpy as np

def plotSlidingRP(spikeTimes, params):
    # Check if params is provided and has 'cidx'
    if 'cidx' not in params:
        params['cidx'] = None

    # Assuming slidingRP is a function that returns the required variables
    (maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
     nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate) = slidingRP(spikeTimes, params)

    # Set up the figure
    f, axs = plt.subplots(1, 3, figsize=(13, 3.69))
    f.subplots_adjust(wspace=0.4)
    plt.rcParams['axes.facecolor'] = 'w'

    # Subplot 1
    axs[0].bar(rp * 1000, nACG, width=1, color='k', edgecolor='none')
    axs[0].set_xlim([0, 5])
    axs[0].set_xlabel('Time from spike (ms)')
    axs[0].set_ylabel('ACG count (spks)')
    axs[0].fill_between([0, 0.5, 0.5, 0], [0, 0, 1, 1] * axs[0].get_ylim()[1], color='k', alpha=0.2)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    if params['cidx'] is not None:
        axs[0].set_title(f'Cluster #{params["cidx"]}: FR={firingRate:.2f}')
    else:
        axs[0].set_title(f'FR={firingRate:.2f}')

    # Subplot 2
    im = axs[1].imshow(confMatrix, aspect='auto', origin='lower', extent=[rp[0] * 1000, rp[-1] * 1000, cont[0], cont[-1]], vmin=0, vmax=100)
    axs[1].plot([rp[0], rp[-1]] * 1000, [10, 10], 'r')
    if not np.isnan(timeOfLowestCont):
        axs[1].plot(timeOfLowestCont * 1000 * np.ones(2), [cont[0], cont[-1]], 'r')
        # Compute and plot the conf=90 contour
        ii = np.argmax(np.vstack([np.zeros(len(rp)), confMatrix > 90]), axis=0)
        ii[ii == 0] = np.nan
        contContour = np.full_like(ii, np.nan, dtype=float)
        contContour[~np.isnan(ii)] = cont[ii[~np.isnan(ii)] - 1]
        axs[1].plot(rp * 1000, contContour, 'r', linewidth=2.0)
    axs[1].fill_between([0, 0.5, 0.5, 0], [0, 0, 1, 1] * axs[1].get_ylim()[1], color='k', alpha=0.5)
    cbar = plt.colorbar(im, ax=axs[1])
    cbar.set_label('Confidence (%)')
    axs[1].set_xlabel('Time from spike (ms)')
    axs[1].set_xlim([0, 5])
    axs[1].set_ylabel('Contamination (%)')
    axs[1].set_ylim([0, max(cont)])
    axs[1].invert_yaxis()
    title_text = f'max conf = {maxConfidenceAt10Cont:.2f}%, min cont = {minContWith90Confidence:.1f}%, time = {timeOfLowestCont*1000:.2f} ms'
    axs[1].set_title(title_text)

    # Subplot 3
    conf_10_cont = confMatrix[cont == 10, :]
    axs[2].plot(rp * 1000, conf_10_cont.T, 'k', linewidth=2.0)
    axs[2].set_xlabel('Time from spike (ms)')
    axs[2].set_ylabel('Confidence of ≤10% contamination (%)')
    axs[2].plot([0, 5], [90, 90], 'r')
    axs[2].fill_between([0, 0.5, 0.5, 0], [0, 0, 1, 1] * 100, color='k', alpha=0.2)
    axs[2].set_xlim([0, 5])
    axs[2].set_ylim([0, 100])
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)

    # Set title colors based on conditions
    if minContWith90Confidence <= 10:
        axs[0].title.set_color([34/255, 177/255, 76/255])
        axs[1].title.set_color([34/255, 177/255, 76/255])
    elif nSpikesBelow2 == 0:
        axs[0].title.set_color('b')
        axs[1].title.set_color('b')
    else:
        axs[0].title.set_color('r')
        axs[1].title.set_color('r')

    plt.show()

# Example usage
# plotSlidingRP(spikeTimes, params)
