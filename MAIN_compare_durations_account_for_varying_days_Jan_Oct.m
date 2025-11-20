%% MRI Study Comparison Between 2023 and 2025 (Jan 1–Oct 12, Hourly Rate + Daily Productivity)
% Jan–Oct version
% - Period: 1.1.–12.10. for both 2023 and 2025
% - 2023 files: *_01_01_2023-12_31_2023.xlsx  (Optima2 for B0_Vida)
% - 2025 files: *_01_01_2025-10_12_2025.xlsx or *_01_01_2025-10_12_2023.xlsx
% - Writes to a separate results folder for this Jan–Oct analysis
%
% Adds:
%  - mean ± SD and median [95% CI] for hourly rates (Jan–Oct)
%  - % change (2025 vs 2023) for mean and median per scanner
%  - option to ignore days with idle gaps > threshold
%  - Optima1 & Optima2 hourly productivity (Jan–Oct 2023)
%  - Fleet-wide daily throughput (patients/hour = Exams / AvailHours)
%  - Fleet EXAM & BETWEEN-PATIENT fractions and % changes
%  - Tables:
%       * scannerTable: all scanners + Optimas
%       * fleetSummaryTable: fleet throughput + EXAM/GAP fractions
%       * scannerModalityTotalTable: scanner-wise + B0 + G0 + Total summary
%       * presentationTable: formatted human-readable summary table

clear; clc;
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%% ===================== User configuration =====================
dataFolder    = '/Users/mikaeljuntunen/Documents_Mikael/Work/Publications/Post-doc/DRB/Practical impact of DLR/Data';
resultsFolder = '/Users/mikaeljuntunen/Documents_Mikael/Work/Publications/Post-doc/DRB/Practical impact of DLR/Results';

% --- Options for daily productivity filtering ---
ignoreDaysWithLargeIdleGaps = true;   % set false to disable filtering
gapThresholdMinutes         = inf;    % threshold for "large" idle gap (minutes)

% --- Jan–Oct date window (1.1.–12.10.) ---
janOctStart2023 = datetime(2023,1,1);
janOctEnd2023   = datetime(2023,10,12);
janOctStart2025 = datetime(2025,1,1);
janOctEnd2025   = datetime(2025,10,12);

% Output folder for this Jan–Oct analysis
baseSubfolderName = 'Results_2025_01_01-10_12';
if ignoreDaysWithLargeIdleGaps
    outputFolder = fullfile(resultsFolder, ...
        sprintf('%s_IgnoreGap%imin', baseSubfolderName, gapThresholdMinutes));
else
    outputFolder = fullfile(resultsFolder, baseSubfolderName);
end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Scanners and data columns
scanners = {'B0_Vida','B0_Sola','G0_Solafit','G0_Vida','G0_Vidafit'};
dtFormat = 'yyyy-MM-dd HH:mm:ss';
variablesToRead = {'Date_Time','Description','x_OfSeries','Duration_min_','MRScanDuration_min_'};

% Codes for boxplots
allowedCodes  = ["AA1BG","AA1CG","KE1CG","NA1BG","NA3BG","NG1BG","NH1BG"];
excludedCodes = ["AA5FG"];  % research etc.

% Default fonts
set(groot, 'DefaultAxesFontName', 'Calibri');
set(groot, 'DefaultTextFontName', 'Calibri');
set(groot, 'DefaultAxesFontSize', 22);

%% ===================== Initialize accumulators (Jan–Oct) =====================
allGrandSummary = table();   % now: Total Duration + MR Scan Duration
capacitySummary = table();
dailyCapacity   = table();  % per-day productivity across scanners

%% ===================== Main loop over scanners (Jan–Oct 2023 vs 2025) =====================
for s = 1:numel(scanners)
    scanner = scanners{s};
    disp("Processing scanner (Jan–Oct): " + scanner);

    %% ---- File paths (2023) ----
    if strcmp(scanner,'B0_Vida')
        % B0 Vida 2023 is Optima2
        file2023 = fullfile(dataFolder,'Optima2_01_01_2023-12_31_2023.xlsx');
    else
        file2023 = fullfile(dataFolder,sprintf('%s_01_01_2023-12_31_2023.xlsx',scanner));
    end

    %% ---- File paths (2025: extended Jan–Oct) ----
    % Try both naming conventions for every scanner:
    %   <scanner>_01_01_2025-10_12_2025.xlsx
    %   <scanner>_01_01_2025-10_12_2023.xlsx
    cand1 = fullfile(dataFolder, sprintf('%s_01_01_2025-10_12_2025.xlsx', scanner));
    cand2 = fullfile(dataFolder, sprintf('%s_01_01_2025-10_12_2023.xlsx', scanner));

    if exist(cand1, 'file') == 2
        file2025 = cand1;
    elseif exist(cand2, 'file') == 2
        file2025 = cand2;
    else
        warning('File not found for 2025 extended period (tried):\n  %s\n  %s', cand1, cand2);
        file2025 = '';  % will lead to empty 2025 data
    end

    %% ---- Read & process 2023 (Jan–Oct) ----
    T2023 = read_study_table(file2023,variablesToRead,dtFormat);
    [desc2023,dur2023,mrdur2023,date2023,code2023] = process_table(T2023);

    % Exclude research / non-routine codes
    mask = ~ismember(code2023, excludedCodes);
    desc2023 = desc2023(mask); dur2023 = dur2023(mask);
    mrdur2023 = mrdur2023(mask); date2023 = date2023(mask); code2023 = code2023(mask);

    % Restrict to 1.1.–12.10.2023
    mask2023 = date2023 >= janOctStart2023 & date2023 <= janOctEnd2023;
    desc2023 = desc2023(mask2023); dur2023 = dur2023(mask2023);
    mrdur2023 = mrdur2023(mask2023); date2023 = date2023(mask2023);
    code2023 = code2023(mask2023);

    %% ---- Read & process 2025 (Jan–Oct) ----
    if isempty(file2025) || exist(file2025,'file') ~= 2
        warning('No valid 2025 file for scanner %s (Jan–Oct).', scanner);
        desc2025 = strings(0,1); dur2025 = []; mrdur2025 = []; date2025 = datetime.empty(0,1); code2025 = strings(0,1);
    else
        T2025 = read_study_table(file2025,variablesToRead,dtFormat);
        [desc2025,dur2025,mrdur2025,date2025,code2025] = process_table(T2025);

        mask = ~ismember(code2025, excludedCodes);
        desc2025 = desc2025(mask); dur2025 = dur2025(mask);
        mrdur2025 = mrdur2025(mask); date2025 = date2025(mask); code2025 = code2025(mask);

        % Restrict to 1.1.–12.10.2025
        mask2025 = date2025 >= janOctStart2025 & date2025 <= janOctEnd2025;
        desc2025 = desc2025(mask2025); dur2025 = dur2025(mask2025);
        mrdur2025 = mrdur2025(mask2025); date2025 = date2025(mask2025);
        code2025 = code2025(mask2025);
    end

    %% ---- Grand totals (per scanner, Jan–Oct) ----
    % === A) Total examination duration (existing) ===
    grand_N2023_tot = numel(dur2023);
    grand_N2025_tot = numel(dur2025);

    grand_median2023_tot = median(dur2023);
    grand_median2025_tot = median(dur2025);

    grand_absChange_tot  = grand_median2025_tot - grand_median2023_tot;
    grand_pctChange_tot  = (grand_median2025_tot - grand_median2023_tot) / grand_median2023_tot * 100;

    if grand_N2023_tot >= 3 && grand_N2025_tot >= 3
        grand_pval_tot = ranksum(dur2023, dur2025);
    else
        grand_pval_tot = NaN;
    end

    grandTable_tot = table(string(scanner), "Total Duration", ...
        grand_N2023_tot, grand_median2023_tot, ...
        grand_N2025_tot, grand_median2025_tot, ...
        grand_absChange_tot, grand_pctChange_tot, grand_pval_tot, ...
        'VariableNames', {'Scanner','DurationType', ...
        'N_2023','Median_2023','N_2025','Median_2025', ...
        'AbsChange','PctChange','Pval'});

    % === B) MR scan duration (NEW, analogous to total) ===
    grand_N2023_mr = numel(mrdur2023);
    grand_N2025_mr = numel(mrdur2025);

    grand_median2023_mr = median(mrdur2023);
    grand_median2025_mr = median(mrdur2025);

    grand_absChange_mr  = grand_median2025_mr - grand_median2023_mr;
    grand_pctChange_mr  = (grand_median2025_mr - grand_median2023_mr) / grand_median2023_mr * 100;

    if grand_N2023_mr >= 3 && grand_N2025_mr >= 3
        grand_pval_mr = ranksum(mrdur2023, mrdur2025);
    else
        grand_pval_mr = NaN;
    end

    grandTable_mr = table(string(scanner), "MR Scan Duration", ...
        grand_N2023_mr, grand_median2023_mr, ...
        grand_N2025_mr, grand_median2025_mr, ...
        grand_absChange_mr, grand_pctChange_mr, grand_pval_mr, ...
        'VariableNames', {'Scanner','DurationType', ...
        'N_2023','Median_2023','N_2025','Median_2025', ...
        'AbsChange','PctChange','Pval'});

    % Append both rows (Total + MR) to the grand summary
    allGrandSummary = [allGrandSummary; grandTable_tot; grandTable_mr];

    %% ---- Availability-normalized mean hourly capacity (Jan–Oct) ----
    [availH_2023, days_2023] = computeAvailabilityHours(date2023, dur2023);
    [availH_2025, days_2025] = computeAvailabilityHours(date2025, dur2025);

    exams_2023 = grand_N2023_tot;
    exams_2025 = grand_N2025_tot;

    examsPerHour_2023 = NaN; examsPerHour_2025 = NaN;
    if availH_2023 > 0, examsPerHour_2023 = exams_2023 / availH_2023; end
    if availH_2025 > 0, examsPerHour_2025 = exams_2025 / availH_2025; end

    pctChange_examsPerHour = (examsPerHour_2025 - examsPerHour_2023) ./ examsPerHour_2023 * 100;

    capRow = table(string(scanner), exams_2023, exams_2025, ...
        days_2023, days_2025, ...
        availH_2023, availH_2025, ...
        examsPerHour_2023, examsPerHour_2025, ...
        pctChange_examsPerHour, ...
        'VariableNames', {'Scanner','Exams_2023','Exams_2025', ...
        'ActiveDays_2023','ActiveDays_2025', ...
        'AvailHours_2023','AvailHours_2025', ...
        'ExamsPerHour_2023','ExamsPerHour_2025', ...
        'PctChange_ExamsPerHour'});

    capacitySummary = [capacitySummary; capRow];

    %% ---- Daily productivity (Jan–Oct) + EXAM/GAP contributions ----
    optsDP = struct('ignoreLargeGaps', ignoreDaysWithLargeIdleGaps, ...
                    'gapThresholdMinutes', gapThresholdMinutes);

    [daily2023, excluded2023] = computeDailyProductivity(date2023, dur2023, optsDP);
    [daily2025, excluded2025] = computeDailyProductivity(date2025, dur2025, optsDP);

    if ~isempty(daily2023)
        daily2023.Scanner = repmat(string(scanner), height(daily2023), 1);
        daily2023.Period  = repmat("JanOct_2023",   height(daily2023), 1);

        dailyCapacity = [dailyCapacity; ...
            daily2023(:, {'Scanner','Period','Day','Exams','AvailHours','ExamsPerHour', ...
                          'ExamHours','GapHours','ExamFrac','GapFrac'})];
    end

    if ~isempty(daily2025)
        daily2025.Scanner = repmat(string(scanner), height(daily2025), 1);
        daily2025.Period  = repmat("JanOct_2025",   height(daily2025), 1);

        dailyCapacity = [dailyCapacity; ...
            daily2025(:, {'Scanner','Period','Day','Exams','AvailHours','ExamsPerHour', ...
                          'ExamHours','GapHours','ExamFrac','GapFrac'})];
    end

    fprintf('[%s] Excluded days (Jan–Oct 2023): %d | (Jan–Oct 2025): %d  (idle gap > %d min, ignore=%s)\n', ...
        scanner, excluded2023, excluded2025, gapThresholdMinutes, string(ignoreDaysWithLargeIdleGaps));

    %% ---- Boxplots for durations (Jan–Oct) ----
    plotBox(scanner, allowedCodes, dur2023, dur2025, code2023, code2025, ...
        "Total Study Duration (min)", ...
        fullfile(outputFolder, sprintf('Boxplot_%s_TotalDurations_JanOct.png', scanner)));

    plotBox(scanner, allowedCodes, mrdur2023, mrdur2025, code2023, code2025, ...
        "MR Scan Duration (min)", ...
        fullfile(outputFolder, sprintf('Boxplot_%s_MRDurations_JanOct.png', scanner)));
end

%% === Daily productivity for Optima1 (Jan–Oct 2023) ==================
nExamsOpt1_2023      = NaN;
availH_Optima1_2023  = NaN;

fileOptima1 = fullfile(dataFolder,'Optima1_01_01_2023-12_31_2023.xlsx');
if exist(fileOptima1,'file') == 2
    TOpt1 = read_study_table(fileOptima1,variablesToRead,dtFormat);
    [~, durOpt1, ~, dateOpt1, codeOpt1] = process_table(TOpt1);

    mask = ~ismember(codeOpt1, excludedCodes);
    durOpt1  = durOpt1(mask);
    dateOpt1 = dateOpt1(mask);

    maskJanOct = dateOpt1 >= janOctStart2023 & dateOpt1 <= janOctEnd2023;
    durOpt1  = durOpt1(maskJanOct);
    dateOpt1 = dateOpt1(maskJanOct);

    nExamsOpt1_2023 = numel(durOpt1);
    [availH_Optima1_2023, ~] = computeAvailabilityHours(dateOpt1, durOpt1);

    optsDP = struct('ignoreLargeGaps', ignoreDaysWithLargeIdleGaps, ...
                    'gapThresholdMinutes', gapThresholdMinutes);
    [dailyOptima1, excludedOpt1] = computeDailyProductivity(dateOpt1, durOpt1, optsDP);

    fprintf('[Optima1] Excluded days (Jan–Oct 2023): %d  (idle gap > %d min, ignore=%s)\n', ...
        excludedOpt1, gapThresholdMinutes, string(ignoreDaysWithLargeIdleGaps));

    if ~isempty(dailyOptima1)
        dailyOptima1.Scanner = repmat("Optima1", height(dailyOptima1), 1);
        dailyOptima1.Period  = repmat("JanOct_2023", height(dailyOptima1), 1);
        dailyCapacity = [dailyCapacity; ...
            dailyOptima1(:, {'Scanner','Period','Day','Exams','AvailHours','ExamsPerHour', ...
                             'ExamHours','GapHours','ExamFrac','GapFrac'})];
    end
else
    warning('Optima1 file not found: %s', fileOptima1);
    dailyOptima1 = table(); %#ok<NASGU>
end

%% ---- Save per-scanner summaries (Jan–Oct) ----
writetable(allGrandSummary, fullfile(outputFolder, 'Grand_Totals_All_Scanners_JanOct2025.xlsx'));
writetable(capacitySummary, fullfile(outputFolder, 'Capacity_Summary_All_Scanners_JanOct2025.xlsx'));
writetable(dailyCapacity,   fullfile(outputFolder, 'Daily_Productivity_All_Scanners_JanOct2025.xlsx'));
disp('[Jan–Oct] All summaries saved.');

%% ---- Stats + Command Window printout (+ % changes by scanner, Jan–Oct) ----
statsTable = summarizeHourlyStats_JanOct(dailyCapacity, scanners);
writetable(statsTable, fullfile(outputFolder, 'Hourly_Rate_Stats_ByScanner_JanOct2025.xlsx'));
disp('[Jan–Oct] Hourly rate stats saved (mean±SD, median with 95% CI).');

fprintf('\n========== [Jan–Oct] Hourly Rate Summary (with %% changes) ==========\n');
for s = 1:numel(scanners)
    sc = string(scanners{s});
    t  = statsTable(statsTable.Scanner == sc, :);

    t23 = t(t.Period=="JanOct_2023", :);
    t25 = t(t.Period=="JanOct_2025", :);

    n23 = getVal(t23,'N');     n25 = getVal(t25,'N');
    mu23 = getVal(t23,'Mean'); mu25 = getVal(t25,'Mean');
    sd23 = getVal(t23,'SD');   sd25 = getVal(t25,'SD');
    md23 = getVal(t23,'Median'); md25 = getVal(t25,'Median');
    lo23 = getVal(t23,'Median_CI_Low'); hi23 = getVal(t23,'Median_CI_High');
    lo25 = getVal(t25,'Median_CI_Low'); hi25 = getVal(t25,'Median_CI_High');

    pctMean   = NaN;
    pctMedian = NaN;
    if ~isnan(mu23) && mu23~=0 && ~isnan(mu25)
        pctMean = 100*(mu25 - mu23)/mu23;
    end
    if ~isnan(md23) && md23~=0 && ~isnan(md25)
        pctMedian = 100*(md25 - md23)/md23;
    end

    fprintf('--- %s (Jan–Oct) ---\n', sc);
    fprintf('2023 Jan–Oct | N=%3d | mean = %.3f ± %.3f | median = %.3f [%.3f, %.3f]\n', n23, mu23, sd23, md23, lo23, hi23);
    fprintf('2025 Jan–Oct | N=%3d | mean = %.3f ± %.3f | median = %.3f [%.3f, %.3f]\n', n25, mu25, sd25, md25, lo25, hi25);
    fprintf('Δ 2025 vs 2023 | mean: %+6.2f%% | median: %+6.2f%%\n\n', pctMean, pctMedian);
end
fprintf('=====================================================================\n\n');

%% === Hourly productivity for Optima1 & Optima2 (Jan–Oct 2023) =======

% Optima1: from dailyOptima1 (if available)
if exist('dailyOptima1','var') && ~isempty(dailyOptima1)
    rO1 = dailyOptima1.ExamsPerHour;
    [nO1, muO1, sdO1, medO1, loO1, hiO1] = summarizeVector(rO1);

    fprintf('--- Optima1 (Jan–Oct 2023) ---\n');
    fprintf('N=%3d | mean = %.3f ± %.3f | median = %.3f [%.3f, %.3f]\n\n', ...
        nO1, muO1, sdO1, medO1, loO1, hiO1);
else
    fprintf('--- Optima1 (Jan–Oct 2023) ---\nNo valid daily productivity data.\n\n');
    nO1 = NaN; muO1 = NaN; sdO1 = NaN; medO1 = NaN; loO1 = NaN; hiO1 = NaN;
end

% Optima2: 2023 data are used as B0_Vida JanOct_2023 in this script
tOpt2 = statsTable(statsTable.Scanner=="B0_Vida" & statsTable.Period=="JanOct_2023", :);
if ~isempty(tOpt2)
    nO2  = tOpt2.N;
    muO2 = tOpt2.Mean;
    sdO2 = tOpt2.SD;
    mdO2 = tOpt2.Median;
    loO2 = tOpt2.Median_CI_Low;
    hiO2 = tOpt2.Median_CI_High;

    fprintf('--- Optima2 (Jan–Oct 2023) ---\n');
    fprintf('N=%3d | mean = %.3f ± %.3f | median = %.3f [%.3f, %.3f]\n\n', ...
        nO2, muO2, sdO2, mdO2, loO2, hiO2);
else
    fprintf('--- Optima2 (Jan–Oct 2023) ---\nNo valid stats (check B0_Vida 2023 data).\n\n');
    nO2 = NaN; muO2 = NaN; sdO2 = NaN; mdO2 = NaN; loO2 = NaN; hiO2 = NaN;
end

%% === Combined scanner + Optima hourly productivity table (Jan–Oct) ===
rowsAB = {};

for s = 1:numel(scanners)
    sc = string(scanners{s});
    t_sc = statsTable(statsTable.Scanner == sc, :);
    if isempty(t_sc), continue; end

    t23 = t_sc(t_sc.Period=="JanOct_2023", :);
    t25 = t_sc(t_sc.Period=="JanOct_2025", :);

    if ~isempty(t23)
        rowsAB = [rowsAB; {sc, "JanOct_2023", ...
            t23.N, t23.Mean, t23.SD, t23.Median, ...
            t23.Median_CI_Low, t23.Median_CI_High, ...
            NaN, NaN}]; %#ok<AGROW>
    end

    if ~isempty(t25)
        pctMean = NaN;
        pctMed  = NaN;
        if ~isempty(t23) && ~isnan(t23.Mean) && t23.Mean ~= 0 && ~isnan(t25.Mean)
            pctMean = 100*(t25.Mean - t23.Mean)/t23.Mean;
        end
        if ~isempty(t23) && ~isnan(t23.Median) && t23.Median ~= 0 && ~isnan(t25.Median)
            pctMed = 100*(t25.Median - t23.Median)/t23.Median;
        end

        rowsAB = [rowsAB; {sc, "JanOct_2025", ...
            t25.N, t25.Mean, t25.SD, t25.Median, ...
            t25.Median_CI_Low, t25.Median_CI_High, ...
            pctMean, pctMed}]; %#ok<AGROW>
    end
end

if ~isnan(nO1)
    rowsAB = [rowsAB; {"Optima1","JanOct_2023", ...
        nO1, muO1, sdO1, medO1, loO1, hiO1, ...
        NaN, NaN}]; %#ok<AGROW>
end

if ~isnan(nO2)
    rowsAB = [rowsAB; {"Optima2","JanOct_2023", ...
        nO2, muO2, sdO2, mdO2, loO2, hiO2, ...
        NaN, NaN}]; %#ok<AGROW>
end

scannerTable = cell2table(rowsAB, 'VariableNames', ...
    {'Scanner','Period','N','Mean','SD','Median', ...
     'Median_CI_Low','Median_CI_High', ...
     'PctChange_Mean','PctChange_Median'});

outFile_scanners = fullfile(outputFolder, 'Scanner_Optima_Hourly_Productivity_JanOct2025.xlsx');
writetable(scannerTable, outFile_scanners);
fprintf('Scanner + Optima hourly productivity table (Jan–Oct) saved to:\n  %s\n\n', outFile_scanners);

%% === Fleet-level hourly productivity & EXAM/GAP fractions (Jan–Oct) ===

fleet2023 = dailyCapacity(dailyCapacity.Period=="JanOct_2023", :);
fleet2025 = dailyCapacity(dailyCapacity.Period=="JanOct_2025", :);

G23 = table(); G25 = table();
if ~isempty(fleet2023)
    G23 = groupsummary(fleet2023, "Day", "sum", ...
        ["Exams","ExamHours","GapHours","AvailHours"]);
end
if ~isempty(fleet2025)
    G25 = groupsummary(fleet2025, "Day", "sum", ...
        ["Exams","ExamHours","GapHours","AvailHours"]);
end

if ~isempty(G23)
    valid23 = G23.sum_AvailHours > 0;
    dailyFleetRate23  = G23.sum_Exams(valid23) ./ G23.sum_AvailHours(valid23); % patients/hour
    dailyExamHours23  = G23.sum_ExamHours(valid23);
    dailyGapHours23   = G23.sum_GapHours(valid23);
    dailyAvailHours23 = G23.sum_AvailHours(valid23);
    dailyExamFrac23   = dailyExamHours23 ./ dailyAvailHours23;
    dailyGapFrac23    = dailyGapHours23  ./ dailyAvailHours23;
else
    dailyFleetRate23  = [];
    dailyExamFrac23   = [];
    dailyGapFrac23    = [];
end

if ~isempty(G25)
    valid25 = G25.sum_AvailHours > 0;
    dailyFleetRate25  = G25.sum_Exams(valid25) ./ G25.sum_AvailHours(valid25);
    dailyExamHours25  = G25.sum_ExamHours(valid25);
    dailyGapHours25   = G25.sum_GapHours(valid25);
    dailyAvailHours25 = G25.sum_AvailHours(valid25);
    dailyExamFrac25   = dailyExamHours25 ./ dailyAvailHours25;
    dailyGapFrac25    = dailyGapHours25  ./ dailyAvailHours25;
else
    dailyFleetRate25  = [];
    dailyExamFrac25   = [];
    dailyGapFrac25    = [];
end

[nFR23, muFR23, sdFR23, medFR23, loFR23, hiFR23] = summarizeVector(dailyFleetRate23);
[nFR25, muFR25, sdFR25, medFR25, loFR25, hiFR25] = summarizeVector(dailyFleetRate25);

[~, muEF23, sdEF23, medEF23, loEF23, hiEF23] = summarizeVector(dailyExamFrac23);
[~, muEF25, sdEF25, medEF25, loEF25, hiEF25] = summarizeVector(dailyExamFrac25);

[~, muGF23, sdGF23, medGF23, loGF23, hiGF23] = summarizeVector(dailyGapFrac23);
[~, muGF25, sdGF25, medGF25, loGF25, hiGF25] = summarizeVector(dailyGapFrac25);

deltaFR_mean   = NaN;
deltaFR_median = NaN;
if ~isnan(muFR23) && muFR23 ~= 0 && ~isnan(muFR25)
    deltaFR_mean = 100*(muFR25 - muFR23)/muFR23;
end
if ~isnan(medFR23) && medFR23 ~= 0 && ~isnan(medFR25)
    deltaFR_median = 100*(medFR25 - medFR23)/medFR23;
end

deltaEF_mean   = NaN;
deltaEF_median = NaN;
if ~isnan(muEF23) && muEF23 ~= 0 && ~isnan(muEF25)
    deltaEF_mean = 100*(muEF25 - muEF23)/muEF23;
end
if ~isnan(medEF23) && medEF23 ~= 0 && ~isnan(medEF25)
    deltaEF_median = 100*(medEF25 - medEF23)/medEF23;
end

deltaGF_mean   = NaN;
deltaGF_median = NaN;
if ~isnan(muGF23) && muGF23 ~= 0 && ~isnan(muGF25)
    deltaGF_mean = 100*(muGF25 - muGF23)/muGF23;
end
if ~isnan(medGF23) && medGF23 ~= 0 && ~isnan(medGF25)
    deltaGF_median = 100*(medGF25 - medGF23)/medGF23;
end

Period = ["JanOct_2023"; "JanOct_2025"];
Ndays  = [nFR23; nFR25];

FleetThroughput_Mean   = [muFR23;   muFR25];
FleetThroughput_SD     = [sdFR23;   sdFR25];
FleetThroughput_Median = [medFR23;  medFR25];
FleetThroughput_CI_Low = [loFR23;   loFR25];
FleetThroughput_CI_High= [hiFR23;   hiFR25];

ExamFrac_Mean   = [muEF23;   muEF25];
ExamFrac_SD     = [sdEF23;   sdEF25];
ExamFrac_Median = [medEF23;  medEF25];
ExamFrac_CI_Low = [loEF23;   loEF25];
ExamFrac_CI_High= [hiEF23;   hiEF25];

GapFrac_Mean    = [muGF23;   muGF25];
GapFrac_SD      = [sdGF23;   sdGF25];
GapFrac_Median  = [medGF23;  medGF25];
GapFrac_CI_Low  = [loGF23;   loGF25];
GapFrac_CI_High = [hiGF23;   hiGF25];

Delta_FleetThroughput_Mean   = [NaN; deltaFR_mean];
Delta_FleetThroughput_Median = [NaN; deltaFR_median];
Delta_ExamFrac_Mean          = [NaN; deltaEF_mean];
Delta_ExamFrac_Median        = [NaN; deltaEF_median];
Delta_GapFrac_Mean           = [NaN; deltaGF_mean];
Delta_GapFrac_Median         = [NaN; deltaGF_median];

fleetSummaryTable = table(Period, Ndays, ...
    FleetThroughput_Mean, FleetThroughput_SD, FleetThroughput_Median, ...
    FleetThroughput_CI_Low, FleetThroughput_CI_High, ...
    ExamFrac_Mean, ExamFrac_SD, ExamFrac_Median, ...
    ExamFrac_CI_Low, ExamFrac_CI_High, ...
    GapFrac_Mean, GapFrac_SD, GapFrac_Median, ...
    GapFrac_CI_Low, GapFrac_CI_High, ...
    Delta_FleetThroughput_Mean, Delta_FleetThroughput_Median, ...
    Delta_ExamFrac_Mean, Delta_ExamFrac_Median, ...
    Delta_GapFrac_Mean, Delta_GapFrac_Median);

outFile_fleet = fullfile(outputFolder, 'Fleet_Productivity_Exam_Gap_Fractions_JanOct2025.xlsx');
writetable(fleetSummaryTable, outFile_fleet);
fprintf('Fleet throughput + EXAM/GAP fraction table (Jan–Oct) saved to:\n  %s\n\n', outFile_fleet);

%% === Scanner / B0 / G0 / Total summary table (Jan–Oct) ===

rowsSummary = {};

% ---------- 1) Per-scanner rows (core scanners) ----------
for s = 1:numel(scanners)
    sc = string(scanners{s});

    % Volumes from capacitySummary
    capRow = capacitySummary(capacitySummary.Scanner == sc, :);
    if isempty(capRow)
        Exams23 = NaN; Exams25 = NaN;
    else
        Exams23 = capRow.Exams_2023;
        Exams25 = capRow.Exams_2025;
    end

    % Hourly rate stats from statsTable (per scanner)
    t23 = statsTable(statsTable.Scanner == sc & statsTable.Period == "JanOct_2023", :);
    t25 = statsTable(statsTable.Scanner == sc & statsTable.Period == "JanOct_2025", :);

    % Default NaNs
    mu23 = NaN; sd23 = NaN; med23 = NaN; lo23 = NaN; hi23 = NaN;
    mu25 = NaN; sd25 = NaN; med25 = NaN; lo25 = NaN; hi25 = NaN;

    if ~isempty(t23)
        mu23  = t23.Mean;
        sd23  = t23.SD;
        med23 = t23.Median;
        lo23  = t23.Median_CI_Low;
        hi23  = t23.Median_CI_High;
    end
    if ~isempty(t25)
        mu25  = t25.Mean;
        sd25  = t25.SD;
        med25 = t25.Median;
        lo25  = t25.Median_CI_Low;
        hi25  = t25.Median_CI_High;
    end

    % Special rule: B0_Vida has no 2023 results in this summary (Optima2 instead)
    if sc == "B0_Vida"
        Exams23 = NaN;
        mu23 = NaN; sd23 = NaN; med23 = NaN; lo23 = NaN; hi23 = NaN;
    end

    % Percentage changes (2025 vs 2023)
    deltaMean = NaN;
    deltaMed  = NaN;
    if ~isnan(mu23) && mu23 ~= 0 && ~isnan(mu25)
        deltaMean = 100 * (mu25 - mu23) / mu23;
    end
    if ~isnan(med23) && med23 ~= 0 && ~isnan(med25)
        deltaMed = 100 * (med25 - med23) / med23;
    end

    rowsSummary = [rowsSummary; { ...
        sc, ...
        Exams23, Exams25, ...
        mu23, sd23, med23, lo23, hi23, ...
        mu25, sd25, med25, lo25, hi25, ...
        deltaMean, deltaMed ...
    }]; %#ok<AGROW>
end

% ---------- 2) Optima1 and Optima2 rows (B department legacy) ----------
opt2CapRow = capacitySummary(capacitySummary.Scanner=="B0_Vida", :);
if ~isempty(opt2CapRow)
    ExamsOpt2_2023 = opt2CapRow.Exams_2023;
    AvailOpt2_2023 = opt2CapRow.AvailHours_2023;
else
    ExamsOpt2_2023 = NaN;
    AvailOpt2_2023 = NaN;
end

if exist('nO1','var') && (~isnan(nO1) || ~isnan(nExamsOpt1_2023))
    rowsSummary = [rowsSummary; { ...
        "Optima1", ...
        nExamsOpt1_2023, NaN, ...
        muO1, sdO1, medO1, loO1, hiO1, ...
        NaN, NaN, NaN, NaN, NaN, ...
        NaN, NaN ...
    }]; %#ok<AGROW>
end

if exist('nO2','var') && (~isnan(nO2) || ~isnan(ExamsOpt2_2023))
    rowsSummary = [rowsSummary; { ...
        "Optima2", ...
        ExamsOpt2_2023, NaN, ...
        muO2, sdO2, mdO2, loO2, hiO2, ...
        NaN, NaN, NaN, NaN, NaN, ...
        NaN, NaN ...
    }]; %#ok<AGROW>
end

% ---------- 3) Aggregated groups: B0 imaging, G0 imaging, Total ----------

% --- B0 imaging (includes Optima1 & Optima2 in 2023, Sola + Vida in 2025) ---
B0Cap     = capacitySummary(startsWith(capacitySummary.Scanner, "B0_"), :);
B0VidaRow = B0Cap(B0Cap.Scanner=="B0_Vida", :);
B0SolaRow = B0Cap(B0Cap.Scanner=="B0_Sola", :);

Exams_B0_2023 = 0;
Exams_B0_2025 = 0;
Avail_B0_2023 = 0;
Avail_B0_2025 = 0;

if ~isempty(B0SolaRow)
    Exams_B0_2023 = Exams_B0_2023 + B0SolaRow.Exams_2023;
    Exams_B0_2025 = Exams_B0_2025 + B0SolaRow.Exams_2025;
    Avail_B0_2023 = Avail_B0_2023 + B0SolaRow.AvailHours_2023;
    Avail_B0_2025 = Avail_B0_2025 + B0SolaRow.AvailHours_2025;
end

% Optima1 (2023 only)
if ~isnan(nExamsOpt1_2023) && ~isnan(availH_Optima1_2023)
    Exams_B0_2023 = Exams_B0_2023 + nExamsOpt1_2023;
    Avail_B0_2023 = Avail_B0_2023 + availH_Optima1_2023;
end

% Optima2 (2023 only)
if exist('ExamsOpt2_2023','var') && ~isnan(ExamsOpt2_2023) && ...
   exist('AvailOpt2_2023','var') && ~isnan(AvailOpt2_2023)
    Exams_B0_2023 = Exams_B0_2023 + ExamsOpt2_2023;
    Avail_B0_2023 = Avail_B0_2023 + AvailOpt2_2023;
end

% B0 Vida (2025 only)
if ~isempty(B0VidaRow)
    Exams_B0_2025 = Exams_B0_2025 + B0VidaRow.Exams_2025;
    Avail_B0_2025 = Avail_B0_2025 + B0VidaRow.AvailHours_2025;
end

% Department capacity as sum of scanner-level exams/hour
B0_rate23_parts = [];
B0_rate25_parts = [];

if ~isempty(B0SolaRow)
    if ~isnan(B0SolaRow.ExamsPerHour_2023)
        B0_rate23_parts(end+1) = B0SolaRow.ExamsPerHour_2023; %#ok<AGROW>
    end
    if ~isnan(B0SolaRow.ExamsPerHour_2025)
        B0_rate25_parts(end+1) = B0SolaRow.ExamsPerHour_2025; %#ok<AGROW>
    end
end

if ~isempty(B0VidaRow)
    if ~isnan(B0VidaRow.ExamsPerHour_2025)
        B0_rate25_parts(end+1) = B0VidaRow.ExamsPerHour_2025; %#ok<AGROW>
    end
end

% Optima1 (2023 only)
if ~isnan(nExamsOpt1_2023) && ~isnan(availH_Optima1_2023) && availH_Optima1_2023>0
    B0_rate23_parts(end+1) = nExamsOpt1_2023 / availH_Optima1_2023; %#ok<AGROW>
end

% Optima2 (2023 only)
if exist('ExamsOpt2_2023','var') && ~isnan(ExamsOpt2_2023) && ...
   exist('AvailOpt2_2023','var') && ~isnan(AvailOpt2_2023) && AvailOpt2_2023>0
    B0_rate23_parts(end+1) = ExamsOpt2_2023 / AvailOpt2_2023; %#ok<AGROW>
end

if ~isempty(B0_rate23_parts)
    B0_rate23 = sum(B0_rate23_parts, 'omitnan');
else
    B0_rate23 = NaN;
end

if ~isempty(B0_rate25_parts)
    B0_rate25 = sum(B0_rate25_parts, 'omitnan');
else
    B0_rate25 = NaN;
end

B0_deltaMean = NaN;
if ~isnan(B0_rate23) && B0_rate23 ~= 0 && ~isnan(B0_rate25)
    B0_deltaMean = 100*(B0_rate25 - B0_rate23)/B0_rate23;
end

rowsSummary = [rowsSummary; { ...
    "B0 imaging", ...
    Exams_B0_2023, Exams_B0_2025, ...
    B0_rate23, NaN, NaN, NaN, NaN, ...
    B0_rate25, NaN, NaN, NaN, NaN, ...
    B0_deltaMean, NaN ...
}]; %#ok<AGROW>

% --- G0 imaging ---
G0Cap         = capacitySummary(startsWith(capacitySummary.Scanner, "G0_"), :);
Exams_G0_2023 = sum(G0Cap.Exams_2023,      'omitnan');
Exams_G0_2025 = sum(G0Cap.Exams_2025,      'omitnan');
Avail_G0_2023 = sum(G0Cap.AvailHours_2023, 'omitnan');
Avail_G0_2025 = sum(G0Cap.AvailHours_2025, 'omitnan');

G0_rate23 = NaN;
G0_rate25 = NaN;
if ~isempty(G0Cap)
    if any(~isnan(G0Cap.ExamsPerHour_2023))
        G0_rate23 = sum(G0Cap.ExamsPerHour_2023, 'omitnan');
    end
    if any(~isnan(G0Cap.ExamsPerHour_2025))
        G0_rate25 = sum(G0Cap.ExamsPerHour_2025, 'omitnan');
    end
end

G0_deltaMean = NaN;
if ~isnan(G0_rate23) && G0_rate23 ~= 0 && ~isnan(G0_rate25)
    G0_deltaMean = 100*(G0_rate25 - G0_rate23)/G0_rate23;
end

rowsSummary = [rowsSummary; { ...
    "G0 imaging", ...
    Exams_G0_2023, Exams_G0_2025, ...
    G0_rate23, NaN, NaN, NaN, NaN, ...
    G0_rate25, NaN, NaN, NaN, NaN, ...
    G0_deltaMean, NaN ...
}]; %#ok<AGROW>

% --- Total fleet: B0 + G0 ---
Exams_Total_2023 = Exams_B0_2023 + Exams_G0_2023;
Exams_Total_2025 = Exams_B0_2025 + Exams_G0_2025;
Avail_Total_2023 = Avail_B0_2023 + Avail_G0_2023;
Avail_Total_2025 = Avail_B0_2025 + Avail_G0_2025;

Tot_rate23 = NaN;
Tot_rate25 = NaN;
if ~isnan(B0_rate23) || ~isnan(G0_rate23)
    Tot_rate23 = sum([B0_rate23, G0_rate23], 'omitnan');
end
if ~isnan(B0_rate25) || ~isnan(G0_rate25)
    Tot_rate25 = sum([B0_rate25, G0_rate25], 'omitnan');
end

Tot_deltaMean = NaN;
if ~isnan(Tot_rate23) && Tot_rate23 ~= 0 && ~isnan(Tot_rate25)
    Tot_deltaMean = 100*(Tot_rate25 - Tot_rate23)/Tot_rate23;
end

rowsSummary = [rowsSummary; { ...
    "Total", ...
    Exams_Total_2023, Exams_Total_2025, ...
    Tot_rate23, NaN, NaN, NaN, NaN, ...
    Tot_rate25, NaN, NaN, NaN, NaN, ...
    Tot_deltaMean, NaN ...
}]; %#ok<AGROW>

% ---------- 4) Build table and save ----------
scannerModalityTotalTable = cell2table(rowsSummary, 'VariableNames', { ...
    'Group', ...
    'Exams_2023', 'Exams_2025', ...
    'MeanRate_2023', 'SDRate_2023', 'MedianRate_2023', ...
    'MedianCI_Low_2023', 'MedianCI_High_2023', ...
    'MeanRate_2025', 'SDRate_2025', 'MedianRate_2025', ...
    'MedianCI_Low_2025', 'MedianCI_High_2025', ...
    'Delta_MeanRate_pct', 'Delta_MedianRate_pct' });

outFile_group = fullfile(outputFolder, 'Scanner_B0_G0_Total_HourlySummary_JanOct2025.xlsx');
writetable(scannerModalityTotalTable, outFile_group);

fprintf('Scanner/B0/G0/Total summary table (incl. Optimas, B0_Vida 2025-only, Jan–Oct) saved to:\n  %s\n\n', outFile_group);

% ---------- 5) Nicely formatted presentation table (indented years) ----------
rowsPres = {};

for i = 1:height(scannerModalityTotalTable)
    r = scannerModalityTotalTable(i,:);
    scannerName = string(r.Group);

    % 2023 row
    have23 = ~(isnan(r.Exams_2023) & isnan(r.MeanRate_2023) & isnan(r.MedianRate_2023));
    if have23
        exams23   = r.Exams_2023;
        meanSD23  = formatMeanSD(r.MeanRate_2023,  r.SDRate_2023);
        medCI23   = formatMedianCI(r.MedianRate_2023, ...
                                   r.MedianCI_Low_2023, r.MedianCI_High_2023);

        rowsPres(end+1, :) = {scannerName, "Jan–Oct 2023", exams23, ... %#ok<AGROW>
                              meanSD23, medCI23, "", ""};
        scannerName = "";
    end

    % 2025 row
    have25 = ~(isnan(r.Exams_2025) & isnan(r.MeanRate_2025) & isnan(r.MedianRate_2025));
    if have25
        exams25   = r.Exams_2025;
        meanSD25  = formatMeanSD(r.MeanRate_2025,  r.SDRate_2025);
        medCI25   = formatMedianCI(r.MedianRate_2025, ...
                                   r.MedianCI_Low_2025, r.MedianCI_High_2025);

        dMeanStr  = "";
        dMedStr   = "";
        if ~isnan(r.Delta_MeanRate_pct)
            dMeanStr = sprintf('%+.2f', r.Delta_MeanRate_pct);
        end
        if ~isnan(r.Delta_MedianRate_pct)
            dMedStr = sprintf('%+.2f', r.Delta_MedianRate_pct);
        end

        rowsPres(end+1, :) = {scannerName, "Jan–Oct 2025", exams25, ... %#ok<AGROW>
                              meanSD25, medCI25, string(dMeanStr), string(dMedStr)};
    end
end

presentationTable = cell2table(rowsPres, 'VariableNames', ...
    {'Scanner','Period','Exams', ...
     'Mean_SD','Median_CI','Delta_Mean_pct','Delta_Median_pct'});

outFile_pres = fullfile(outputFolder, 'Scanner_B0_G0_Total_HourlySummary_Presentation_JanOct2025.xlsx');
writetable(presentationTable, outFile_pres);

fprintf('Presentation-style hourly summary (Jan–Oct) saved to:\n  %s\n\n', outFile_pres);

%% ---- Plot daily hourly productivity per scanner (Jan–Oct) ----
plotDailyProductivityAllScanners_JanOct(dailyCapacity, scanners, ...
    fullfile(outputFolder, 'Daily_Productivity_ByScanner_JanOct2025.png'));
disp('Daily productivity figure (Jan–Oct) saved.');

%% ===================== Helper functions =====================
function s = formatMeanSD(mu, sd)
    if isnan(mu)
        s = "";
    elseif isnan(sd) || sd == 0
        s = string(sprintf('%.2f', mu));
    else
        s = string(sprintf('%.2f ± %.2f', mu, sd));
    end
end

function s = formatMedianCI(med, lo, hi)
    if isnan(med)
        s = "";
    elseif isnan(lo) || isnan(hi)
        s = string(sprintf('%.2f', med));
    else
        s = string(sprintf('%.2f [%.2f, %.2f]', med, lo, hi));
    end
end

function T = read_study_table(file, vars, fmt)
    if exist(file,'file') ~= 2
        error('File not found: %s', file);
    end
    opts = detectImportOptions(file);
    opts.SelectedVariableNames = vars;
    opts = setvartype(opts, vars, 'char');
    T = readtable(file, opts);
    T.DateTimeParsed = datetime(T.Date_Time, 'InputFormat', fmt);
end

function [desc, dur, mrdur, date, code] = process_table(T)
    studyMarkers = ~cellfun(@isempty, T.x_OfSeries);
    studyGroup = cumsum(studyMarkers);
    nRows = height(T);

    durs   = nan(nRows,1);
    mrdurs = nan(nRows,1);

    for i = 1:nRows
        s = T.Duration_min_{i};
        m = T.MRScanDuration_min_{i};
        if ~isempty(s) && contains(s, ':')
            p = split(s, ':');
            durs(i) = str2double(p{1}) + str2double(p{2})/60;
        end
        if ~isempty(m) && contains(m, ':')
            p = split(m, ':');
            mrdurs(i) = str2double(p{1}) + str2double(p{2})/60;
        end
    end

    uStudies = unique(studyGroup);
    nStudies = numel(uStudies);

    desc = strings(nStudies,1);
    code = strings(nStudies,1);
    dur  = nan(nStudies,1);
    mrdur= nan(nStudies,1);
    date = NaT(nStudies,1);

    for i = 1:nStudies
        mask = (studyGroup == uStudies(i));
        idxMarker = find(mask & studyMarkers, 1, 'first');
        txt = strtrim(T.Description{idxMarker});

        if strlength(txt) == 0
            desc(i) = "Unknown";
            code(i) = "Unknown";
        else
            desc(i) = txt;
            if strlength(txt) >= 5
                code(i) = extractBetween(txt, 1, 5);
            else
                code(i) = "Unknown";
            end
        end

        dur(i)   = sum(durs(mask),   'omitnan');
        mrdur(i) = sum(mrdurs(mask), 'omitnan');
        date(i)  = T.DateTimeParsed(find(mask,1,'first'));
    end
end

function plotBox(scanner, codes, val23, val25, code23, code25, label, fname)
    f = figure('Name',char(scanner + " - " + label),'Position',[100,100,1200,500], 'Color','w','Theme','light');

    allDur = [];
    allGrp = strings(0,1);
    allYr  = strings(0,1);

    allDur = [allDur; val23(:); val25(:)];
    allGrp = [allGrp; repmat("Total", numel(val23) + numel(val25),1)];
    allYr  = [allYr;  repmat("2023 Jan–Oct", numel(val23),1); repmat("2025 Jan–Oct", numel(val25),1)];

    for i = 1:numel(codes)
        c = codes(i);
        d23 = val23(code23 == c);
        d25 = val25(code25 == c);
        if ~isempty(d23) || ~isempty(d25)
            allDur = [allDur; d23(:); d25(:)];
            allGrp = [allGrp; repmat(c, numel(d23)+numel(d25),1)];
            allYr  = [allYr;  repmat("2023 Jan–Oct",numel(d23),1); repmat("2025 Jan–Oct",numel(d25),1)];
        end
    end

    allGrp(strlength(allGrp)==0) = "Unknown";
    cats = ["Total"; codes(:)];
    allGrp = categorical(allGrp, cats);

    boxchart(allGrp, allDur, 'GroupByColor', allYr);
    xlabel('Code');
    ylabel(label);
    legend('2023 Jan–Oct','2025 Jan–Oct','Location','northeast','FontSize',16);
    grid on;
    ylim([0,120]);

    exportgraphics(f, fname, 'Resolution',300);
    close(f);
end

function [totalHours, nDays] = computeAvailabilityHours(startTimes, durationsMin)
    if isempty(startTimes)
        totalHours = 0; nDays = 0; return;
    end

    startTimes   = startTimes(:);
    durationsMin = durationsMin(:);
    valid = ~isnat(startTimes) & ~isnan(durationsMin);
    startTimes   = startTimes(valid);
    durationsMin = durationsMin(valid);

    if isempty(startTimes)
        totalHours = 0; nDays = 0; return;
    end

    dayOnly    = dateshift(startTimes,'start','day');
    uniqueDays = unique(dayOnly);
    nDays = numel(uniqueDays);

    totalHours = 0;
    for i = 1:nDays
        mask = (dayOnly == uniqueDays(i));
        firstStart = min(startTimes(mask));
        lastEnd    = max(startTimes(mask) + minutes(durationsMin(mask)));
        hrs = hours(lastEnd - firstStart);
        if hrs > 0.25
            totalHours = totalHours + min(hrs, 24);
        end
    end
end

function [dailyTable, nExcludedDays] = computeDailyProductivity(startTimes, durationsMin, opts)
    if nargin < 3 || isempty(opts), opts = struct; end
    if ~isfield(opts, 'ignoreLargeGaps'), opts.ignoreLargeGaps = false; end
    if ~isfield(opts, 'gapThresholdMinutes'), opts.gapThresholdMinutes = 30; end

    dailyTable = table();
    nExcludedDays = 0;

    if isempty(startTimes)
        return;
    end

    startTimes   = startTimes(:);
    durationsMin = durationsMin(:);
    valid = ~isnat(startTimes) & ~isnan(durationsMin);
    startTimes   = startTimes(valid);
    durationsMin = durationsMin(valid);

    if isempty(startTimes)
        return;
    end

    dayOnly    = dateshift(startTimes,'start','day');
    uniqueDays = unique(dayOnly);

    nD = numel(uniqueDays);
    Day          = NaT(nD,1);
    Exams        = nan(nD,1);
    AvailHours   = nan(nD,1);
    ExamsPerHour = nan(nD,1);

    ExamHours    = nan(nD,1);
    GapHours     = nan(nD,1);
    ExamFrac     = nan(nD,1);
    GapFrac      = nan(nD,1);

    row = 0;
    for i = 1:nD
        d = uniqueDays(i);
        dmask = (dayOnly == d);
        if ~any(dmask), continue; end

        starts    = startTimes(dmask);
        durations = durationsMin(dmask);
        ends      = starts + minutes(durations);

        [starts, order] = sort(starts);
        ends            = ends(order);
        durations       = durations(order);

        gaps = starts(2:end) - ends(1:end-1);
        maxGapMin = 0;
        gapMinSum = 0;
        if ~isempty(gaps)
            gapMin   = max(minutes(gaps),0);
            maxGapMin = max(gapMin);
            gapMinSum = sum(gapMin);
        end

        if opts.ignoreLargeGaps && (maxGapMin > opts.gapThresholdMinutes)
            nExcludedDays = nExcludedDays + 1;
            continue;
        end

        firstStart = min(starts);
        lastEnd    = max(ends);
        hrs        = hours(lastEnd - firstStart);
        nEx        = sum(dmask);

        if hrs > 0.25
            row = row + 1;
            Day(row)          = d;
            Exams(row)        = nEx;
            AvailHours(row)   = min(hrs, 24);
            ExamsPerHour(row) = nEx / AvailHours(row);

            examMin           = sum(durations, 'omitnan');
            ExamHours(row)    = examMin / 60;
            GapHours(row)     = gapMinSum / 60;

            if AvailHours(row) > 0
                ExamFrac(row) = ExamHours(row) / AvailHours(row);
                GapFrac(row)  = GapHours(row) / AvailHours(row);
            end
        end
    end

    validRows = ~isnat(Day);
    dailyTable = table(Day(validRows), Exams(validRows), AvailHours(validRows), ...
        ExamsPerHour(validRows), ExamHours(validRows), GapHours(validRows), ...
        ExamFrac(validRows), GapFrac(validRows), ...
        'VariableNames', {'Day','Exams','AvailHours','ExamsPerHour', ...
                          'ExamHours','GapHours','ExamFrac','GapFrac'});
end

function plotDailyProductivityAllScanners_JanOct(dailyCapacity, scanners, fname)
    if isempty(dailyCapacity)
        warning('No daily capacity data to plot.');
        return;
    end

    nScanners = numel(scanners);
    rowH = 200;
    f = figure('Name','Daily hourly productivity by scanner (Jan–Oct)', ...
               'Position',[100, 100, 1000, max(320, nScanners*rowH)], ...
               'Color','w','Theme','light');

    tiledlayout(nScanners, 1, 'Padding','compact', 'TileSpacing','compact');

    col23 = [0 0.4470 0.7410];      % blue
    col25 = [0.8500 0.3250 0.0980]; % orange
    legendsAdded = false;

    for i = 1:nScanners
        nexttile; hold on;

        sc = string(scanners{i});
        data_sc = dailyCapacity(dailyCapacity.Scanner == sc, :);

        r2023 = data_sc.ExamsPerHour(data_sc.Period == "JanOct_2023");
        r2025 = data_sc.ExamsPerHour(data_sc.Period == "JanOct_2025");

        h = gobjects(0);

        if ~isempty(r2023)
            b1 = boxchart(ones(size(r2023)), r2023, ...
                'BoxFaceColor', col23, 'MarkerColor', col23, 'WhiskerLineColor', col23);
            b1.JitterOutliers = 'on';
            b1.DisplayName = '2023 Jan–Oct';
            h(end+1) = b1;
        end

        if ~isempty(r2025)
            b2 = boxchart(2*ones(size(r2025)), r2025, ...
                'BoxFaceColor', col25, 'MarkerColor', col25, 'WhiskerLineColor', col25);
            b2.JitterOutliers = 'on';
            b2.DisplayName = '2025 Jan–Oct';
            h(end+1) = b2;
        end

        xlim([0.5 2.5]);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'2023 Jan–Oct','2025 Jan–Oct'});
        ylabel('Exams/hour');
        title(strrep(scanners{i}, '_', ' '));
        grid on;

        yl = ylim; yText = yl(2) - 0.05*(yl(2)-yl(1));
        if ~isempty(r2023)
            med23 = median(r2023,'omitnan');
            text(1, yText, sprintf('%.2f', med23), 'Color', col23, ...
                'HorizontalAlignment','center','VerticalAlignment','top', ...
                'FontSize', 8, 'FontWeight','bold');
        end
        if ~isempty(r2025)
            med25 = median(r2025,'omitnan');
            text(2, yText, sprintf('%.2f', med25), 'Color', col25, ...
                'HorizontalAlignment','center','VerticalAlignment','top', ...
                'FontSize', 8, 'FontWeight','bold');
        end

        hold off;

        if i==1 && ~legendsAdded && ~isempty(h)
            lg = legend(h, 'Location','northoutside', 'Orientation','horizontal', ...
                        'Box','off', 'FontSize',10);
            lg.Layout.Tile = 'north';
            legendsAdded = true;
        end
    end

    exportgraphics(f, fname, 'Resolution', 300);
    close(f);
end

function statsTable = summarizeHourlyStats_JanOct(dailyCapacity, scanners)
    rows = {};
    periods = ["JanOct_2023","JanOct_2025"];
    for s = 1:numel(scanners)
        sc = string(scanners{s});
        for per = periods
            r = dailyCapacity.ExamsPerHour(dailyCapacity.Scanner==sc & dailyCapacity.Period==per);
            r = r(~isnan(r));
            if isempty(r)
                row = {sc, per, 0, NaN, NaN, NaN, NaN, NaN};
            else
                [n, mu, sd, med, lo, hi] = summarizeVector(r);
                row = {sc, per, n, mu, sd, med, lo, hi};
            end
            rows = [rows; row]; %#ok<AGROW>
        end
    end
    statsTable = cell2table(rows, 'VariableNames', ...
        {'Scanner','Period','N','Mean','SD','Median','Median_CI_Low','Median_CI_High'});
end

function [med, ciLow, ciHigh] = medianWithCI(x, alpha, nBoot)
    if nargin<2 || isempty(alpha), alpha = 0.05; end
    if nargin<3 || isempty(nBoot), nBoot = 2000; end

    x = x(~isnan(x));
    if isempty(x)
        med = NaN; ciLow = NaN; ciHigh = NaN; return;
    end
    med = median(x);
    if numel(x) < 3
        ciLow = NaN; ciHigh = NaN; return;
    end
    if exist('bootstrp','file') == 2
        bootMeds = bootstrp(nBoot, @median, x);
    else
        bootMeds = zeros(nBoot,1);
        n = numel(x);
        for b = 1:nBoot
            idx = randi(n, n, 1);
            bootMeds(b) = median(x(idx));
        end
    end
    q = quantile(bootMeds, [alpha/2, 1-alpha/2]);
    ciLow  = q(1);
    ciHigh = q(2);
end

function v = getVal(T, name)
    if isempty(T) || height(T)==0 || all(ismissing(T.(name)))
        v = NaN;
    else
        v = T.(name)(1);
    end
end

function [n, mu, sd, med, lo, hi] = summarizeVector(x)
    x = x(~isnan(x));
    n = numel(x);
    if n == 0
        mu = NaN; sd = NaN; med = NaN; lo = NaN; hi = NaN;
        return;
    end
    mu = mean(x,'omitnan');
    sd = std(x,'omitnan');
    [med, lo, hi] = medianWithCI(x, 0.05, 2000);
end
