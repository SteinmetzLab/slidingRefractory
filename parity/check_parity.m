% CHECK_PARITY  Compare the MATLAB implementation with the Python one.
%
%   Recomputes every quantity stored by make_cases.py (Sliding RP scalars,
%   the full confidence matrix, the FWER-corrected confidence, Hill-Llobet,
%   tauPass0 and minPassingFR) on the same spike trains and option sets, and
%   reports the worst disagreement for each. Errors at the end if any
%   comparison fails its tolerance.
%
%   Requires histdiff (cortex-lab/spikes) on the path, or set spikesPath.
%   Run from anywhere:  run('<repo>/parity/check_parity.m')

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo, 'matlab'));
spikesPath = 'D:\Dropbox\code\spikes\analysis\helpers';
if isempty(which('histdiff')) && exist(spikesPath, 'dir'); addpath(spikesPath); end
assert(~isempty(which('histdiff')), 'histdiff (cortex-lab/spikes) must be on the path');

S = load(fullfile(here, 'parity_cases.mat'));
nT = numel(S.trains); nO = numel(S.opt_cont);
names = cellstr(S.names);
fails = {};
note = @(msg) fprintf('  %s\n', msg);

fprintf('\nPython/MATLAB parity: %d trains x %d option sets\n\n', nT, nO);
worst = struct('conf', 0, 'cmin', 0, 'rpmin', 0, 'tp0', 0, 'mat', 0);
nMis = struct('pass', 0, 'nshort', 0, 'nan', 0);
for i = 1:nT
    st = double(S.trains{i}(:));
    for j = 1:nO
        p = struct('recDur', S.recDur(i), 'contaminationThresh', S.opt_cont(j), ...
                   'confidenceThresh', S.opt_conf(j), 'rpReject', S.opt_rpReject(j), ...
                   'censor', S.opt_censor(j));
        [passTest, conf, cmin, rpmin, nShort, confMatrix, ~, ~, ~, tp0] = slidingRP(st, p);
        worst.conf = max(worst.conf, abs(conf - S.py_max_conf(i, j)));
        pc = S.py_min_cont(i, j); pr = S.py_rp_min_val(i, j);
        if xor(isnan(cmin), isnan(pc)) || xor(isnan(rpmin), isnan(pr))
            nMis.nan = nMis.nan + 1;
            note(sprintf('C_min NaN mismatch: %s / %s  (MATLAB %g, Python %g)', ...
                names{i}, S.opt_names{j}, cmin, pc));
        elseif ~isnan(cmin)
            worst.cmin = max(worst.cmin, abs(cmin - pc));
            worst.rpmin = max(worst.rpmin, abs(rpmin - pr));
        end
        if passTest ~= logical(S.py_passes(i, j))
            nMis.pass = nMis.pass + 1;
            note(sprintf('pass mismatch: %s / %s (conf %.10f)', names{i}, S.opt_names{j}, conf));
        end
        if nShort ~= S.py_n_below2(i, j); nMis.nshort = nMis.nshort + 1; end
        worst.tp0 = max(worst.tp0, abs(tp0 - S.py_tau_pass0(i, j)) / max(tp0, eps));
        worst.mat = max(worst.mat, max(abs(confMatrix(:) - S.py_matrix{i, j}(:))));
    end
end

fprintf('%-44s %12s  %s\n', 'quantity', 'worst diff', 'tolerance');
chk = @(name, v, tol) fprintf('%-44s %12.3g  %-8.0e %s\n', name, v, tol, ternary(v <= tol));
chk('max confidence (%)', worst.conf, 1e-8);
chk('C_min (%)', worst.cmin, 1e-6);
chk('tau_r at C_min (s)', worst.rpmin, 1e-12);
chk('tau_pass0 (relative)', worst.tp0, 1e-12);
chk('confidence matrix (%)', worst.mat, 1e-8);
fprintf('%-44s %12d\n', 'pass/fail mismatches', nMis.pass);
fprintf('%-44s %12d\n', 'nViolShort mismatches', nMis.nshort);
fprintf('%-44s %12d\n', 'C_min NaN mismatches', nMis.nan);
if worst.conf > 1e-8;  fails{end+1} = 'max confidence'; end
if worst.cmin > 1e-6;  fails{end+1} = 'C_min'; end
if worst.rpmin > 1e-12; fails{end+1} = 'tau_r at C_min'; end
if worst.tp0 > 1e-12;  fails{end+1} = 'tau_pass0'; end
if worst.mat > 1e-8;   fails{end+1} = 'confidence matrix'; end
if nMis.pass > 0;  fails{end+1} = 'pass/fail'; end
if nMis.nshort > 0; fails{end+1} = 'nViolShort'; end
if nMis.nan > 0;   fails{end+1} = 'C_min NaN pattern'; end

% FWER-corrected confidence
wc = 0;
for a = 1:numel(S.corr_idx)
    i = S.corr_idx(a); st = double(S.trains{i}(:));
    for b = 1:2
        w = [0, 0.00025];
        p = struct('recDur', S.recDur(i), 'censor', w(b), 'correction', true);
        [~, conf] = slidingRP(st, p);
        wc = max(wc, abs(conf - S.py_corrected(a, b)));
    end
end
chk('FWER-corrected confidence (%)', wc, 1e-6);
if wc > 1e-6; fails{end+1} = 'FWER-corrected confidence'; end

% Hill-Llobet, spike-time path
wp = 0; we = 0;
metrics = cellstr(S.hl_metric);
for i = 1:nT
    st = double(S.trains{i}(:));
    for k = 1:numel(S.hl_rp)
        p = struct('metricType', metrics{k}, 'RPdur', S.hl_rp(k), 'recDur', S.recDur(i), ...
                   'censor', S.hl_censor(k), 'contaminationThresh', 10);
        [ps, est] = RPmetric_Classic(st, p);
        wp = wp + (ps ~= logical(S.py_hl_pass(i, k)));
        pe = S.py_hl_est(i, k);
        if ~(isnan(est) && isnan(pe)); we = max(we, abs(real(est) - pe)); end
    end
end
fprintf('%-44s %12d\n', 'Hill-Llobet pass/fail mismatches', wp);
chk('Hill-Llobet estimate', we, 1e-12);
if wp > 0; fails{end+1} = 'Hill-Llobet pass/fail'; end
if we > 1e-12; fails{end+1} = 'Hill-Llobet estimate'; end

% power functions
g = S.pow_grid;
tpm = arrayfun(@(r) tauPass0(g(r,1), g(r,2), g(r,3), g(r,4), g(r,5)), 1:size(g,1))';
frm = arrayfun(@(r) minPassingFR(g(r,2), 0.002, g(r,3), g(r,4), g(r,5)), 1:size(g,1))';
wtp = max(abs(tpm - S.py_pow_tau_pass0(:)) ./ S.py_pow_tau_pass0(:));
wfr = max(abs(frm - S.py_min_fr(:)) ./ S.py_min_fr(:));
chk('tauPass0 function (relative)', wtp, 1e-12);
chk('minPassingFR function (relative)', wfr, 1e-12);
if wtp > 1e-12; fails{end+1} = 'tauPass0 function'; end
if wfr > 1e-12; fails{end+1} = 'minPassingFR function'; end

if isempty(fails)
    fprintf('\nALL PARITY CHECKS PASSED\n');
else
    error('Parity failures: %s', strjoin(fails, ', '));
end

function s = ternary(ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
end
