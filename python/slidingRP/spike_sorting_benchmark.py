##
from one.api import ONE
from brainbox.io.one import SpikeSortingLoader

one = ONE(base_url='https://openalyx.internationalbrainlab.org')

benchmark_pids = ['1a276285-8b0e-4cc9-9f0a-a3a002978724',
                  '1e104bf4-7a24-4624-a5b2-c2c8289c0de7',
                  '5d570bf6-a4c6-4bf1-a14b-2c878c84ef0e',
                  '5f7766ce-8e2e-410c-9195-6bf089fea4fd',
                  '6638cfb3-3831-4fc2-9327-194b76cf22e1',
                  '749cb2b7-e57e-4453-a794-f6230e4d0226',
                  'd7ec0892-0a6c-4f4f-9d8f-72083692af5c',
                  'da8dfec1-d265-44e8-84ce-6ae9c109b8bd',
                  'dab512bd-a02d-4c1f-8dbc-9155a163efc0',
                  'dc7e9403-19f7-409f-9240-05ee57cb7aea',
                  'e8f9fba4-d151-4b00-bee7-447f0f3e752c',
                  'eebcaf65-7fa4-4118-869d-a084e84530e2',
                  'fe380793-8035-414e-b000-09bfe5ece92a']
##
# Select a PID
pid = benchmark_pids[0]

# Load spike sorting
sl = SpikeSortingLoader(pid=pid, one=one)
spikes, clusters, channels = sl.load_spike_sorting()
clusters = sl.merge_clusters(spikes, clusters, channels)

# Get AP spikeglx.Reader objects
sr_ap = sl.raw_electrophysiology(band="ap", stream=True)

# Load raw data
s_event = 100  # timepoint in recording to stream
window_secs_ap = [0, 1]
first, last = (int(window_secs_ap[0] * sr_ap.fs) + s_event, int(window_secs_ap[1] * sr_ap.fs + s_event))
raw_ap = sr_ap[first:last, :-sr_ap.nsync].T

##