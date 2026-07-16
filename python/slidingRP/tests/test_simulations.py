import numpy as np
import pytest
from slidingRP import simulations as S


def test_genST_realised_rate():
    np.random.seed(0)
    st = S.genST(10, 2000, 0.003)
    assert len(st) / 2000 == pytest.approx(10, rel=0.05)


def test_genST_respects_refractory_period():
    np.random.seed(0)
    rp = 0.004
    st = S.genST(20, 1000, rp)
    assert np.min(np.diff(st)) >= rp - 1e-12


def test_genST_zero_rp_allows_short_isis():
    np.random.seed(0)
    st = S.genST(20, 1000, 0)
    assert np.min(np.diff(st)) < 1e-3


def test_RPmetric_Classic_llobet_formula():
    # With a precomputed ACG, estContam must equal the closed-form Llobet
    # inversion and passTest must equal obsViol <= expectedViol.
    rp = np.arange(0, 10 / 1000, 1 / 30000) + (1 / 30000) / 2
    nACG = np.zeros(300)
    nACG[:15] = 1  # 15 violations below RPdur (keeps the inversion real)
    spikeCount, recDur, RPdur = 5000, 3600.0, 0.003
    params = dict(metricType='Llobet', contaminationThresh=10, recDur=recDur,
                  RPdur=RPdur, nACG=nACG, rp=rp, spikeCount=spikeCount)
    passTest, estContam = S.RPmetric_Classic(None, params)

    rpIdx = np.where(rp > RPdur)[0][0]
    obsViol = nACG[:rpIdx + 1].sum()
    Nc, Nb = spikeCount * 0.1, spikeCount * 0.9
    Ve = 2 * RPdur / recDur * Nc * (Nb + (Nc - 1) / 2)
    assert estContam == pytest.approx(1 - np.sqrt(1 - obsViol * recDur / (spikeCount**2 * RPdur)))
    assert passTest == (obsViol <= Ve)


def test_RPmetric_Classic_hill_more_conservative_than_llobet():
    # Hill Ve = Nc*Nb < Llobet Ve = Nc*(Nb+(Nc-1)/2); chosen so the observed
    # violations fall between the two expected counts -> Llobet passes, Hill
    # fails (Hill is the stricter / more conservative test).
    # contThresh=30 widens the Hill/Llobet Ve gap (Hill Ve~23.3, Llobet Ve~28.3)
    # so obsViol=26 lands between them: Llobet passes, Hill fails.
    rp = np.arange(0, 10 / 1000, 1 / 30000) + (1 / 30000) / 2
    nACG = np.zeros(300); nACG[:26] = 1  # obsViol = 26 (RPdur = 2 ms)
    base = dict(contaminationThresh=30, recDur=3600.0, RPdur=0.002,
                nACG=nACG, rp=rp, spikeCount=10000)
    llob = S.RPmetric_Classic(None, {**base, 'metricType': 'Llobet'})[0]
    hill = S.RPmetric_Classic(None, {**base, 'metricType': 'Hill'})[0]
    assert llob is True and hill is False
    assert not (hill and not llob)  # in general: Hill pass implies Llobet pass


def test_run_simulations_smoke():
    np.random.seed(0)
    sim = S.run_simulations(dict(
        nSim=40, baseRates=np.array([5.0]),
        contProps=np.array([0.0, 0.10, 0.30]),
        RPdurs=np.array([0.003]), recDurs=np.array([7200.0]),
        confThreshes=np.array([90]), contThresh=10), progress=False)
    # clean -> all pass; heavily contaminated -> none pass
    order = np.argsort(sim['cont_prop'])
    pp = sim['passPct'][order]
    assert pp[0] == 100.0
    assert pp[-1] == 0.0
    assert set(sim).issuperset({'passPct', 'passPctLlobet3', 'passPctHill3',
                                'total_rate', 'cont_prop', 'conf_level'})
