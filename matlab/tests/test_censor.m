classdef test_censor < matlab.unittest.TestCase
% TEST_CENSOR  Tests for params.censor, the sorter censor (duplicate-removal)
% window. Mirrors python/slidingRP/tests/test_censor.py. Requires histdiff
% (cortex-lab/spikes) on the path for the tests that build an ACG.
%
%   results = runtests('matlab/tests/test_censor.m')

    properties (Constant)
        REPO_ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        DUR = 3600;
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            addpath(fullfile(testCase.REPO_ROOT, 'matlab'));
            addpath(fullfile(testCase.REPO_ROOT, 'matlab', 'simulations'));
        end
    end

    methods (Static)
        function st = train(rate, cont, rp, seed, censor)
            rng(seed);
            base = genST(rate * (1 - cont), 3600, rp, []);
            if cont > 0; other = genST(rate * cont, 3600, 0, []); else; other = []; end
            st = sort([base(:); other(:)]);
            if censor > 0
                keep = true(size(st)); last = -Inf;
                for i = 1:numel(st)
                    if st(i) - last < censor; keep(i) = false; else; last = st(i); end
                end
                st = st(keep);
            end
        end
    end

    methods (Test, TestTags={'power'})
        function test_power_functions_shift_by_the_censor(testCase)
            n = 18000; w = 0.00025;
            testCase.verifyEqual(tauPass0(n, 3600, 10, 90, w), tauPass0(n, 3600) + w, 'RelTol', 1e-12);
            testCase.verifyEqual(minPassingFR(3600, 0.003, 10, 90, w), ...
                minPassingFR(3600, 0.003 - w), 'RelTol', 1e-12);
            testCase.verifyEqual(minPassingFR(3600, 0.0002, 10, 90, w), Inf);
        end
    end

    methods (Test, TestTags={'slidingRP'})
        function test_default_is_censor_zero(testCase)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, 'histdiff not on path');
            st = test_censor.train(5, 0.10, 0.003, 1, 0);
            a = cell(1, 10); b = cell(1, 10);
            [a{:}] = slidingRP(st, struct('recDur', testCase.DUR));
            [b{:}] = slidingRP(st, struct('recDur', testCase.DUR, 'censor', 0));
            for k = 1:10
                testCase.verifyEqual(b{k}, a{k});
            end
        end

        function test_censor_uses_the_observable_window(testCase)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, 'histdiff not on path');
            w = 0.0005;
            st = test_censor.train(8, 0.05, 0.002, 2, w);
            [cm, cont, rp, nACG] = computeMatrix(st, struct('recDur', testCase.DUR, 'censor', w));
            refDur = rp + (rp(2) - rp(1)) / 2;
            i = find(abs(cont - 10) < 1e-9, 1);
            expect = 100 * computeViol(cumsum(nACG), [], numel(st), max(refDur - w, 0), 0.10, testCase.DUR);
            testCase.verifyEqual(cm(i, :), expect, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(all(cm(:, refDur <= w) == 0)));
        end

        function test_ignoring_a_censor_is_anti_conservative(testCase)
            % Units at 12% contamination with a 0.5 ms censor: accepted almost
            % always by the published method, rarely once it is accounted for.
            testCase.assumeTrue(exist('histdiff', 'file') > 0, 'histdiff not on path');
            w = 0.0005; nSim = 60; ign = 0; cor = 0;
            for k = 1:nSim
                st = test_censor.train(5, 0.12, 0.003, 100 + k, w);
                ign = ign + slidingRP(st, struct('recDur', testCase.DUR));
                cor = cor + slidingRP(st, struct('recDur', testCase.DUR, 'censor', w));
            end
            testCase.verifyGreaterThan(ign / nSim, 0.9);
            testCase.verifyLessThan(cor / nSim, 0.3);
        end

        function test_hill_llobet_censor_shortens_the_window(testCase)
            testCase.assumeTrue(exist('histdiff', 'file') > 0, 'histdiff not on path');
            st = test_censor.train(5, 0.10, 0.003, 5, 0.0005);
            p = struct('recDur', testCase.DUR, 'RPdur', 0.003);
            [~, e0] = RPmetric_Classic(st, p);
            p.censor = 0.0005;
            [~, ec] = RPmetric_Classic(st, p);
            testCase.verifyGreaterThan(real(ec), real(e0));
        end
    end
end
