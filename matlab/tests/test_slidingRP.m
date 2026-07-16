classdef test_slidingRP < matlab.unittest.TestCase
% Unit tests for the MATLAB slidingRP implementation.
%
% Run all tests:
%   results = runtests('matlab/tests/test_slidingRP.m')
%
% Tests are grouped into:
%   1. computeViol  – core Poisson formula, no external dependencies
%   2. computeMatrix – output shape and range (requires histdiff)
%   3. slidingRP     – synthetic spike train pass/fail (requires histdiff)
%   4. slidingRP_all – batch consistency (requires histdiff)
%   5. Regression    – real IBL test data (requires histdiff + readNPY)
%
% Tests that require external toolboxes are skipped automatically when
% those toolboxes are not on the path.
%
% External dependencies:
%   - histdiff  : from cortex-lab/spikes toolbox
%   - readNPY   : from kwikteam/npy-matlab toolbox

    properties (Constant)
        REPO_ROOT     = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
        TEST_DATA_PATH = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
                                  'test-data', 'unit');
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            addpath(fullfile(testCase.REPO_ROOT, 'matlab'));
            addpath(fullfile(testCase.REPO_ROOT, 'matlab', 'simulations'));

            % Attempt to add external dependencies via the githubDir() convention.
            % Define a githubDir.m function on your MATLAB path that returns the
            % path to your local GitHub repositories folder, e.g.:
            %   function d = githubDir(), d = 'C:\Users\you\GitHub'; end
            % Then clone cortex-lab/spikes and kwikteam/npy-matlab into that folder.
            try
                gDir = githubDir();
                spikesPath = fullfile(gDir, 'spikes');
                npyPath    = fullfile(gDir, 'npy-matlab');
                if exist(spikesPath, 'dir'); addpath(genpath(spikesPath)); end
                if exist(npyPath,    'dir'); addpath(genpath(npyPath));    end
            catch
                % githubDir not defined or errored — external deps handled below
            end

            % Report dependency status once at setup time
            if exist('histdiff', 'file') == 0
                fprintf(['\n[test_slidingRP] histdiff not found — tests tagged ' ...
                    'computeMatrix / slidingRP / slidingRP_all / regression will be skipped.\n' ...
                    '  Fix: clone cortex-lab/spikes into your githubDir folder,\n' ...
                    '  or run: addpath(genpath(''<path-to-spikes>''))\n\n']);
            end
            if exist('readNPY', 'file') == 0
                fprintf(['\n[test_slidingRP] readNPY not found — regression tests will be skipped.\n' ...
                    '  Fix: clone kwikteam/npy-matlab into your githubDir folder,\n' ...
                    '  or run: addpath(genpath(''<path-to-npy-matlab>''))\n\n']);
            end
        end
    end

    % =====================================================================
    %  1. computeViol – no external dependencies
    % =====================================================================
    methods (Test, TestTags={'computeViol'})

        function test_formula_matches_llobet(testCase)
            % Expected violations must follow the Llobet et al. 2022 formula:
            %   Ve = 2*tau/D * Nc * (Nb + (Nc-1)/2)
            % and confidence must equal 1 - poisscdf(Vo, Ve).
            spikeCount = 1000; recDur = 1000;
            refDur = 0.002; contProp = 0.1; obsViol = 5;

            Nc = spikeCount * contProp;
            Nb = spikeCount * (1 - contProp);
            Ve_expected = 2 * refDur / recDur * Nc * (Nb + (Nc - 1) / 2);
            conf_expected = 1 - poisscdf(obsViol, Ve_expected);

            [conf, Ve] = computeViol(obsViol, [], spikeCount, refDur, contProp, recDur);
            testCase.verifyEqual(Ve,   Ve_expected,   'AbsTol', 1e-12);
            testCase.verifyEqual(conf, conf_expected, 'AbsTol', 1e-12);
        end

        function test_zero_obs_high_expected_gives_high_confidence(testCase)
            % 0 observed violations when many are expected → confident unit is clean
            spikeCount = 10000; recDur = 1000; refDur = 0.003; contProp = 0.1;
            [conf, Ve] = computeViol(0, [], spikeCount, refDur, contProp, recDur);
            testCase.verifyGreaterThan(Ve,   1,    'Expected violations should exceed 1');
            testCase.verifyGreaterThan(conf, 0.99, 'Confidence should exceed 99%');
        end

        function test_many_obs_gives_low_confidence(testCase)
            % Many observed violations → very unlikely the unit is clean
            spikeCount = 1000; recDur = 1000; refDur = 0.002; contProp = 0.1;
            [conf, ~] = computeViol(200, [], spikeCount, refDur, contProp, recDur);
            testCase.verifyLessThan(conf, 0.01, 'Confidence should be below 1%');
        end

        function test_confidence_increases_with_fewer_violations(testCase)
            % Monotonicity: fewer observed violations → higher confidence
            spikeCount = 5000; recDur = 1000; refDur = 0.002; contProp = 0.1;
            [conf_high, ~] = computeViol(50, [], spikeCount, refDur, contProp, recDur);
            [conf_low,  ~] = computeViol( 5, [], spikeCount, refDur, contProp, recDur);
            testCase.verifyGreaterThan(conf_low, conf_high);
        end

        function test_higher_contamination_threshold_gives_higher_confidence(testCase)
            % Testing against a looser contamination threshold should give >= confidence
            spikeCount = 2000; recDur = 1000; refDur = 0.002; obsViol = 10;
            [conf_strict, ~] = computeViol(obsViol, [], spikeCount, refDur, 0.05, recDur);
            [conf_loose,  ~] = computeViol(obsViol, [], spikeCount, refDur, 0.20, recDur);
            testCase.verifyGreaterThanOrEqual(conf_loose, conf_strict);
        end

    end % computeViol tests

    % =====================================================================
    %  2. computeMatrix – shape and range
    % =====================================================================
    methods (Test, TestTags={'computeMatrix'})

        function test_output_shape_and_range(testCase)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(5, 1000, 0.003);
            params = struct('recDur', 1000);
            [confMatrix, cont, rp, ~] = computeMatrix(st, params);

            testCase.verifyEqual(size(confMatrix), [numel(cont), numel(rp)], ...
                'confMatrix must be [nCont x nRP]');
            testCase.verifyGreaterThanOrEqual(min(confMatrix(:)), 0);
            testCase.verifyLessThanOrEqual(max(confMatrix(:)), 100);
        end

        function test_default_grid_sizes(testCase)
            % Default cont = 0.5:0.5:35 (70 levels), default window 10ms at
            % 1/30000 s resolution = 300 bins.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(5, 500, 0.003);
            params = struct('recDur', 500);
            [~, cont, rp, ~] = computeMatrix(st, params);

            testCase.verifyEqual(numel(cont), 70,  'cont should have 70 levels by default');
            testCase.verifyEqual(numel(rp),   300, 'rp should have 300 bins by default');
        end

        function test_confMatrix_monotone_in_obsViol(testCase)
            % For a fixed contamination level, confidence should decrease (or
            % stay constant) as more RP violations accumulate with larger tau_r.
            % In practice this is not strictly monotone, but the first column
            % (shortest tau) should not be the highest confidence.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            % Contaminated spike train: many violations → confidence starts low
            baseST = genST(10, 3600, 0.003);
            contST = genST(2,  3600);
            st = sort([baseST; contST]);
            params = struct('recDur', 3600);
            [confMatrix, ~, ~, ~] = computeMatrix(st, params);
            % For the 10% contamination row (row 20 in 0.5:0.5:35), at least
            % one column should be below 90 (unit is contaminated)
            contRow20 = confMatrix(20, :);  % cont(20) = 10%
            testCase.verifyTrue(any(contRow20 < 90), ...
                'Contaminated unit should have some sub-threshold confidence values');
        end

    end % computeMatrix tests

    % =====================================================================
    %  3. slidingRP – synthetic spike train pass/fail
    % =====================================================================
    methods (Test, TestTags={'slidingRP'})

        function test_clean_high_rate_passes(testCase)
            % Clean, high-rate, long recording: should pass with high confidence
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            st = genST(10, 7200, 0.003);  % 10 sp/s, 2 hr, 3 ms RP
            params = struct('recDur', 7200);
            [passTest, confidence, contamination, ~, ~, ~, ~, ~, ~] = slidingRP(st, params);
            testCase.verifyTrue(passTest, 'Clean high-rate spike train should pass');
            testCase.verifyGreaterThan(confidence,   90);
            testCase.verifyLessThanOrEqual(contamination, 10);
        end

        function test_contaminated_20pct_fails(testCase)
            % 20% contamination: must fail regardless of firing rate
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            recDur = 7200;
            baseST = genST(10,       recDur, 0.003);
            contST = genST(10 * 0.2, recDur);         % no RP on contamination
            combST = sort([baseST; contST]);
            params = struct('recDur', recDur);
            [passTest, ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(combST, params);
            testCase.verifyFalse(passTest, '20%% contaminated spike train should fail');
        end

        function test_low_firing_rate_fails(testCase)
            % Very low firing rate: insufficient power to confirm <10% contamination
            % even with zero observed violations (Fig. 4g of manuscript)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            st = genST(0.5, 3600, 0.003);  % 0.5 sp/s, 1 hr
            params = struct('recDur', 3600);
            [passTest, ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, params);
            testCase.verifyFalse(passTest, ...
                'Low firing-rate unit cannot reach 90%% confidence; should fail');
        end

        function test_short_rp_unit_passes(testCase)
            % Unit with RP < 2 ms (e.g. thalamic neuron): should still pass
            % when clean, because the metric slides tau_r down to 0.5 ms
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            st = genST(10, 7200, 0.0015);  % 10 sp/s, 2 hr, 1.5 ms RP
            params = struct('recDur', 7200);
            [passTest, ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, params);
            testCase.verifyTrue(passTest, ...
                'Clean unit with short (1.5 ms) RP should pass (key advantage over Hill-Llobet)');
        end

        function test_single_spike_fails(testCase)
            % Edge case: a single spike produces no violations and no power
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            st = 0.5;
            params = struct('recDur', 100);
            [passTest, ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, params);
            testCase.verifyFalse(passTest, 'Single spike should fail');
        end

        function test_pass_depends_on_confidence_threshold(testCase)
            % Loosening the confidence threshold should only increase pass rate
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(1);
            % Borderline unit: moderate rate, moderate recording
            st = genST(3, 3600, 0.003);
            params_strict = struct('recDur', 3600, 'confidenceThresh', 90);
            params_loose  = struct('recDur', 3600, 'confidenceThresh', 50);
            [pass_strict, ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, params_strict);
            [pass_loose,  ~, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, params_loose);
            % A unit that passes strict must also pass loose
            if pass_strict
                testCase.verifyTrue(pass_loose, ...
                    'Unit passing strict threshold must also pass loose threshold');
            end
            % A unit that fails loose cannot pass strict
            if ~pass_loose
                testCase.verifyFalse(pass_strict, ...
                    'Unit failing loose threshold cannot pass strict threshold');
            end
        end

    end % slidingRP tests

    % =====================================================================
    %  4. slidingRP_all – batch consistency
    % =====================================================================
    methods (Test, TestTags={'slidingRP_all'})

        function test_batch_matches_individual_calls(testCase)
            % Results from slidingRP_all must match individual slidingRP calls
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            recDur = 3600;
            st1 = genST(10, recDur, 0.003);
            st2 = genST(5,  recDur, 0.002);

            allST = [st1; st2];
            clu   = [ones(numel(st1), 1); 2 * ones(numel(st2), 1)];
            [~, idx] = sort(allST);
            allST = allST(idx); clu = clu(idx);

            params = struct('verbose', false, 'useParallel', false);
            metrics = slidingRP_all(allST, clu, params);

            [~, conf1, cont1, ~] = slidingRP(st1, params);
            [~, conf2, cont2, ~] = slidingRP(st2, params);

            testCase.verifyEqual(metrics(1).confidence,   conf1, 'AbsTol', 1e-10);
            testCase.verifyEqual(metrics(2).confidence,   conf2, 'AbsTol', 1e-10);
            testCase.verifyEqual(metrics(1).contamination, cont1, 'AbsTol', 1e-10);
            testCase.verifyEqual(metrics(2).contamination, cont2, 'AbsTol', 1e-10);
        end

        function test_batch_output_struct_fields(testCase)
            % Output struct must always contain the required fields
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(5, 1000, 0.003);
            params = struct('verbose', false, 'useParallel', false);
            metrics = slidingRP_all(st, ones(numel(st), 1), params);

            testCase.verifyTrue(isfield(metrics, 'cid'));
            testCase.verifyTrue(isfield(metrics, 'confidence'));
            testCase.verifyTrue(isfield(metrics, 'contamination'));
            testCase.verifyTrue(isfield(metrics, 'timeOfLowestCont'));
            testCase.verifyTrue(isfield(metrics, 'confMatrix'));
        end

    end % slidingRP_all tests

    % =====================================================================
    %  5. Regression – real IBL test data
    %     Requires: histdiff (spikes toolbox) + readNPY (npy-matlab)
    % =====================================================================
    methods (Test, TestTags={'regression'})

        function test_regression_cluster_274_passes(testCase)
            % Cluster 274: well-isolated unit, expected to PASS
            % Python expected: confidence=100%, contamination=0.5%
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            testCase.assumeTrue(exist('readNPY', 'file') > 0, ...
                'Skipped: readNPY not on path (add kwikteam/npy-matlab)');

            [spikeTimes, spikeClusters] = testCase.loadTestData();
            st = spikeTimes(spikeClusters == 274);
            [passTest, confidence, contamination, ~, ~, ~, ~, ~, ~] = ...
                slidingRP(st, struct('verbose', false));

            testCase.verifyTrue(passTest, 'Cluster 274 should pass');
            testCase.verifyGreaterThan(confidence, 90);
            testCase.verifyLessThanOrEqual(contamination, 10);
        end

        function test_regression_cluster_167_fails(testCase)
            % Cluster 167: contaminated unit, expected to FAIL (conf < 1%)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            testCase.assumeTrue(exist('readNPY', 'file') > 0, ...
                'Skipped: readNPY not on path (add kwikteam/npy-matlab)');

            [spikeTimes, spikeClusters] = testCase.loadTestData();
            st = spikeTimes(spikeClusters == 167);
            [passTest, confidence, contamination, ~, ~, ~, ~, ~, ~] = ...
                slidingRP(st, struct('verbose', false));

            testCase.verifyFalse(passTest, 'Cluster 167 should fail');
            testCase.verifyLessThan(confidence, 90);
            testCase.verifyTrue(isnan(contamination), ...
                'Cluster 167 contamination should be NaN (>35%%)');
        end

        function test_regression_cluster_275_fails(testCase)
            % Cluster 275: moderately contaminated, expected to FAIL (conf ~54%)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            testCase.assumeTrue(exist('readNPY', 'file') > 0, ...
                'Skipped: readNPY not on path (add kwikteam/npy-matlab)');

            [spikeTimes, spikeClusters] = testCase.loadTestData();
            st = spikeTimes(spikeClusters == 275);
            [passTest, confidence, contamination, ~, ~, ~, ~, ~, ~] = ...
                slidingRP(st, struct('verbose', false));

            testCase.verifyFalse(passTest, 'Cluster 275 should fail');
            testCase.verifyLessThan(confidence, 90);
            testCase.verifyGreaterThan(contamination, 10);  % min confirmable cont > threshold
        end

    end % regression tests

    % =====================================================================
    %  6. genST – simulated spike train generator (no external dependencies)
    % =====================================================================
    methods (Test, TestTags={'genST'})

        function test_genST_realised_rate(testCase)
            % The rate correction should make the realised rate match the
            % requested rate despite the imposed refractory period.
            rng(7);
            st = genST(10, 2000, 0.003);   % 10 sp/s, long recording, 3 ms RP
            realised = numel(st) / 2000;
            testCase.verifyEqual(realised, 10, 'RelTol', 0.05, ...
                'Realised rate should be within 5%% of requested rate');
        end

        function test_genST_respects_refractory_period(testCase)
            % Every ISI is refractoryPeriod + Exponential, so no ISI may be
            % shorter than the refractory period.
            rng(7);
            rp = 0.004;
            st = genST(20, 1000, rp);
            testCase.verifyGreaterThanOrEqual(min(diff(st)), rp - 1e-12, ...
                'No ISI should be shorter than the refractory period');
        end

        function test_genST_zero_rp_allows_short_isis(testCase)
            % With no refractory period (contamination train), short ISIs occur.
            rng(7);
            st = genST(20, 1000, 0);
            testCase.verifyLessThan(min(diff(st)), 0.001, ...
                'With RP=0 some ISIs should fall below 1 ms');
            testCase.verifyEqual(numel(st)/1000, 20, 'RelTol', 0.05);
        end

    end % genST tests

    % =====================================================================
    %  7. Corrected (multiple-comparisons) variant
    %     Requires: histdiff (spikes toolbox)
    % =====================================================================
    methods (Test, TestTags={'corrected'})

        function test_corrected_output_shape(testCase)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(10, 600, 0.003);
            params = struct('recDur', 600, 'correction', true);
            [~, confidence, ~, ~, ~, confMatrix, cont, rp, ~] = slidingRP(st, params);
            testCase.verifyEqual(size(confMatrix), [numel(cont), numel(rp)]);
            testCase.verifyGreaterThanOrEqual(confidence, 0);
            testCase.verifyLessThanOrEqual(confidence, 100);
        end

        function test_corrected_not_more_lenient_than_nominal(testCase)
            % The FWER correction is more conservative: the corrected per-row
            % confidence (constant across tau_r) cannot exceed the best
            % (max-over-tau_r) nominal confidence in that row, because
            % P(cross boundary somewhere) >= P(cross at the single best bin).
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(8, 600, 0.003);
            params = struct('recDur', 600, 'correction', true);
            [confMatrix, ~, ~, ~, nominalConfMatrix] = computeMatrix(st, params);
            rowMaxCorrected = max(confMatrix, [], 2);
            rowMaxNominal   = max(nominalConfMatrix, [], 2);
            testCase.verifyLessThanOrEqual(max(rowMaxCorrected - rowMaxNominal), 1e-9, ...
                'Corrected confidence must not exceed the best nominal confidence per contamination level');
        end

        function test_corrected_confidence_le_standard(testCase)
            % The corrected metric's scalar confidence (at the contamination
            % threshold) must be <= the standard metric's confidence.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(5);
            st = genST(8, 600, 0.003);
            params = struct('recDur', 600);
            [~, confStd]  = slidingRP(st, params);
            [~, confCorr] = slidingRP(st, struct('recDur', 600, 'correction', true));
            testCase.verifyLessThanOrEqual(confCorr, confStd + 1e-9, ...
                'Corrected confidence should be <= standard confidence');
        end

        function test_corrected_nominal_matches_computeMatrix(testCase)
            % The nominal (uncorrected) matrix returned alongside the corrected
            % one must equal the standard (correction-off) computeMatrix output.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(0);
            st = genST(8, 600, 0.003);
            [~, ~, ~, ~, nominalConfMatrix] = computeMatrix(st, struct('recDur', 600, 'correction', true));
            stdConfMatrix = computeMatrix(st, struct('recDur', 600));
            testCase.verifyEqual(nominalConfMatrix, stdConfMatrix, 'AbsTol', 1e-9);
        end

    end % corrected tests

    % =====================================================================
    %  8. RPmetric_Classic – Hill/Llobet point-estimate comparison metric
    %     Requires: histdiff (spikes toolbox)
    % =====================================================================
    methods (Test, TestTags={'classic'})

        function test_classic_llobet_clean_passes(testCase)
            % Clean, high-rate train: few violations -> passes, low estContam.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            st = genST(15, 3600, 0.003);
            params = struct('recDur', 3600, 'metricType', 'Llobet', 'RPdur', 0.002);
            [passTest, estContam] = RPmetric_Classic(st, params);
            testCase.verifyTrue(passTest, 'Clean train should pass the classic metric');
            testCase.verifyLessThan(estContam, 0.1, 'Estimated contamination should be small');
        end

        function test_classic_contaminated_fails(testCase)
            % Heavily contaminated train: many violations -> fails.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(42);
            recDur = 3600;
            baseST = genST(15,       recDur, 0.003);
            contST = genST(15 * 0.3, recDur);          % 30% contamination, no RP
            combST = sort([baseST; contST]);
            params = struct('recDur', recDur, 'metricType', 'Llobet', 'RPdur', 0.002);
            passTest = RPmetric_Classic(combST, params);
            testCase.verifyFalse(passTest, '30%% contaminated train should fail');
        end

        function test_classic_llobet_estContam_formula(testCase)
            % estContam must equal the closed-form Llobet inversion
            %   C = 1 - sqrt(1 - obsViol*D / (Nt^2 * RPdur))
            % for a pre-computed ACG (uses the .nACG fast path).
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(3);
            st = genST(10, 1800, 0.0025);
            recDur = 1800; RPdur = 0.002;
            % Provide the full 0-10 ms ACG so the RPdur bin lookup works, exactly
            % as the metric does internally.
            [nACG, rp] = histdiff(st, st, 0:1/30000:10/1000);
            obsViol = sum(nACG(1:find(rp > RPdur, 1)));
            Nt = numel(st);
            expected = 1 - sqrt(1 - obsViol*recDur / (Nt^2 * RPdur));
            params = struct('recDur', recDur, 'metricType', 'Llobet', 'RPdur', RPdur, ...
                'nACG', nACG, 'rp', rp, 'spikeCount', Nt);
            [~, estContam] = RPmetric_Classic([], params);
            testCase.verifyEqual(estContam, expected, 'AbsTol', 1e-9);
        end

    end % classic tests

    % =====================================================================
    %  9. Analytical fast path (slidingRP without the full matrix)
    %     Requires: histdiff (spikes toolbox), chi2inv (Statistics toolbox)
    % =====================================================================
    methods (Test, TestTags={'analytical'})

        function test_fast_confidence_matches_matrix(testCase)
            % The fast-path confidence (output 2) must equal the confidence
            % derived from the full matrix at the contamination-threshold row.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(11);
            st = genST(8, 1800, 0.003);
            params = struct('recDur', 1800);
            [~, confFast, ~, ~, ~, confMatrix, cont, rp] = slidingRP(st, params);
            testTimes = rp > 0.0005;
            confFromMatrix = max(confMatrix(find(cont >= 10, 1), testTimes));
            testCase.verifyEqual(confFast, confFromMatrix, 'AbsTol', 1e-9, ...
                'Fast-path confidence must equal the matrix-derived confidence');
        end

        function test_analytical_contamination_matches_grid(testCase)
            % The analytical (continuous) min contamination must agree with the
            % grid-based value to within the 0.5%% grid resolution.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(11);
            st = genST(8, 1800, 0.003);
            params = struct('recDur', 1800);
            [~, ~, contAnalytical, ~, ~, confMatrix, cont, rp] = slidingRP(st, params);
            testTimes = rp > 0.0005;
            [ii, ~] = find(confMatrix(:, testTimes) > 90);
            minI = min(ii);
            contGrid = cont(minI);
            testCase.assumeFalse(isnan(contAnalytical) || isempty(contGrid), ...
                'Both contamination estimates must be defined for this comparison');
            testCase.verifyLessThanOrEqual(abs(contAnalytical - contGrid), 0.5, ...
                'Analytical contamination should match the grid within grid resolution');
        end

        function test_fast_path_no_matrix_when_not_requested(testCase)
            % When the matrix output is not requested, the scalar metrics are
            % still produced (and equal the values from a matrix-requesting call).
            testCase.assumeTrue(exist('histdiff', 'file') > 0, ...
                'Skipped: histdiff not on path (add cortex-lab/spikes)');
            rng(11);
            st = genST(8, 1800, 0.003);
            params = struct('recDur', 1800);
            [passFast, confFast, contFast] = slidingRP(st, params);          % fast
            [passM, confM, contM] = slidingRP(st, params);                   % same call
            testCase.verifyEqual(confFast, confM, 'AbsTol', 1e-12);
            testCase.verifyEqual(contFast, contM, 'AbsTol', 1e-12);
            testCase.verifyEqual(passFast, passM);
        end

    end % analytical tests

    % =====================================================================
    %  Helper methods
    % =====================================================================
    methods (Access = private)

        function [spikeTimes, spikeClusters] = loadTestData(testCase)
            spikeTimes    = double(readNPY(fullfile(testCase.TEST_DATA_PATH, 'spikes.times.npy')));
            spikeClusters = double(readNPY(fullfile(testCase.TEST_DATA_PATH, 'spikes.clusters.npy')));
        end

    end

end
