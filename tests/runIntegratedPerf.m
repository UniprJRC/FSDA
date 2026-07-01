% runIntegratedPerf.m — Run FSDAPerformanceSuite via runperf (baseline vs optimized)
%
% Uses matlab.perftest framework for proper warmup and statistical sampling.

fsda_orig = 'C:/Users/anoushn/AppData/Roaming/MathWorks/MATLAB Add-Ons/Toolboxes/FSDA';
integ_opt = 'C:/Users/anoushn/claude-workspace/community-toolboxes/fsda/integrated/optimized';
mex_dir = fullfile(integ_opt, 'mex');
test_dir = 'C:/Users/anoushn/claude-workspace/community-toolboxes/fsda/integrated/tests';

%% BASELINE
fprintf('=== BASELINE ===\n');
restoredefaultpath;
addpath(genpath(fsda_orig));
addpath(test_dir);

baselineResults = runperf('FSDAPerformanceSuite');

fprintf('\n=== OPTIMIZED ===\n');
restoredefaultpath;
addpath(genpath(fsda_orig));
addpath(fullfile(integ_opt, 'multivariate'));
addpath(fullfile(integ_opt, 'regression'));
addpath(fullfile(integ_opt, 'clustering'));
addpath(fullfile(integ_opt, 'utilities'));
addpath(fullfile(integ_opt, 'utilities_stat'));
addpath(mex_dir);
addpath(test_dir);

optimizedResults = runperf('FSDAPerformanceSuite');

%% Compare
fprintf('\n\n==================== RESULTS ====================\n');
fprintf('%-60s %10s %10s %8s\n', 'Test', 'Baseline', 'Optimized', 'Speedup');
fprintf('%s\n', repmat('-', 1, 92));

for i = 1:numel(baselineResults)
    bMean = mean(baselineResults(i).Samples.MeasuredTime);
    oMean = mean(optimizedResults(i).Samples.MeasuredTime);
    sp = bMean / oMean;
    fprintf('%-60s %10.4f %10.4f %7.2fx\n', baselineResults(i).Name, bMean, oMean, sp);
end

save(fullfile(test_dir, 'IntegratedResults.mat'), 'baselineResults', 'optimizedResults');
fprintf('\nResults saved.\n');
