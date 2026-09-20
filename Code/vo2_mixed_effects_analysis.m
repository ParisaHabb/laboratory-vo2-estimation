%% VO2 and Heart-Rate Analysis for Incremental Exercise Tests
% Multi-test analysis with:
% - robust column detection
% - VO2 smoothing
% - incremental-stage segmentation
% - per-stage HR and VO2 features
% - global linear regression
% - linear mixed-effects model
% - leave-one-athlete-out validation
% - residual and Bland-Altman diagnostics
% - Excel export and diagnostic figures
%
% Save this file as:
%   vo2_mixed_effects_analysis.m
%
% Required toolboxes:
%   Statistics and Machine Learning Toolbox for fitlm and fitlme

clear; clc; close all;

%% Configuration
warmupSec = 5 * 60;
recoveryKmh = 5;
tolSpeed = 0.20;
vo2HalfWinSec = 15;

%% Select input files
[fileName, pathName] = uigetfile( ...
    {'*.xlsx;*.xls;*.csv', 'Data files (*.xlsx, *.xls, *.csv)'}, ...
    'Select exercise data files', ...
    'MultiSelect', 'on');

if isequal(fileName, 0)
    error('No input file was selected.');
end

if ischar(fileName) || isstring(fileName)
    fileName = {char(fileName)};
end

nTests = numel(fileName);
Results = struct([]);

%% Process all tests
for kTest = 1:nTests
    currentFile = fileName{kTest};
    filePath = fullfile(pathName, currentFile);
    raw = readcell(filePath);

    if size(raw, 1) < 3
        error('File "%s" is too short. At least two header rows and data are expected.', currentFile);
    end

    headers = lower(strtrim(string(raw(1, :))));

    idxTime = findHeaderIndex(headers, ["t", "time", "timestamp"], 10);
    idxVO2kg = findHeaderIndex(headers, ["vo2/kg", "vo2kg", "vo2 / kg", "vo2"], 22);
    idxHR = findHeaderIndex(headers, ["hr", "heart rate", "heartrate"], 24);
    idxGrade = findHeaderIndex(headers, ["grade", "pendenza", "incline"], 36);
    idxSpeed = findHeaderIndex(headers, ["speed", "kmh", "velocita", "velocity"], 37);

    fprintf('\n=============================================\n');
    fprintf('Test %d/%d: %s\n', kTest, nTests, currentFile);
    fprintf('Time: %d | VO2/kg: %d | HR: %d | Grade: %d | Speed: %d\n', ...
        idxTime, idxVO2kg, idxHR, idxGrade, idxSpeed);

    ageVal = findMetaNumeric(raw, ["età", "eta", "age"]);
    heightVal = findMetaNumeric(raw, ["altezza (cm)", "altezza", "height"]);
    weightVal = findMetaNumeric(raw, ["peso (kg)", "peso", "weight"]);
    hrRestVal = findMetaNumeric(raw, ...
        ["hr rest", "hrrest", "hr_rest", "hr a riposo", "frequenza a riposo", "resting hr"]);

    baseName = stripExtension(currentFile);
    subjectID = findMetaString(raw, ["id1", "subject", "athlete", "id"], baseName);

    dataRows = 3:size(raw, 1);
    tSec = nan(numel(dataRows), 1);
    for i = 1:numel(dataRows)
        tSec(i) = toSeconds(raw{dataRows(i), idxTime});
    end

    vo2kg = colToNumeric(raw(dataRows, idxVO2kg));
    hr = colToNumeric(raw(dataRows, idxHR));
    grade = colToNumeric(raw(dataRows, idxGrade));
    speed = colToNumeric(raw(dataRows, idxSpeed));

    maskKeep = isfinite(tSec) & isfinite(hr) & isfinite(speed);
    tSec = tSec(maskKeep);
    vo2kg = vo2kg(maskKeep);
    hr = hr(maskKeep);
    grade = grade(maskKeep);
    speed = speed(maskKeep);

    [tSec, order] = sort(tSec);
    vo2kg = vo2kg(order);
    hr = hr(order);
    grade = grade(order);
    speed = speed(order);

    if numel(tSec) < 4
        warning('Skipping %s: too few valid samples.', currentFile);
        continue;
    end

    tMin = tSec / 60;
    vo2kgSmooth = movingMeanByTime(tSec, vo2kg, vo2HalfWinSec);
    vo2MaxSmooth = max(vo2kgSmooth, [], 'omitnan');
    hrMax = max(hr, [], 'omitnan');

    signalTable = table( ...
        tSec, tMin, vo2kg, vo2kgSmooth, hr, grade, speed, ...
        'VariableNames', {'time_s', 'time_min', 'VO2kg', 'VO2kg_smooth', ...
        'HR', 'Grade_pct', 'Speed_kmh'});

    [regIntercept, regSpeed, regHR, regSE, regR2, regN] = ...
        fitVO2LinearRegression(signalTable.VO2kg_smooth, signalTable.Speed_kmh, signalTable.HR);

    %% Isolate incremental phase
    idxAfterWarm = find(tSec >= warmupSec, 1, 'first');
    if isempty(idxAfterWarm)
        warning('Skipping %s: no data after the warm-up period.', currentFile);
        continue;
    end

    iStartInc = find(speed(idxAfterWarm:end) > recoveryKmh + tolSpeed, 1, 'first');
    if isempty(iStartInc)
        iStartInc = idxAfterWarm;
    else
        iStartInc = idxAfterWarm + iStartInc - 1;
    end

    peakSpeed = max(speed(iStartInc:end), [], 'omitnan');
    lastPeakOffset = find(abs(speed(iStartInc:end) - peakSpeed) < tolSpeed, 1, 'last');
    if isempty(lastPeakOffset)
        lastPeakOffset = numel(speed) - iStartInc + 1;
    end
    iEndInc = iStartInc + lastPeakOffset - 1;

    maskInc = false(size(tSec));
    maskInc(iStartInc:iEndInc) = true;
    tInc = tSec(maskInc);
    hrInc = hr(maskInc);
    vo2Inc = vo2kgSmooth(maskInc);
    gradeInc = grade(maskInc);
    speedInc = speed(maskInc);

    %% Segment stages
    speedStage = round(speedInc * 10) / 10;
    changeIdx = [1; find(abs(diff(speedStage)) > tolSpeed) + 1; numel(speedStage) + 1];

    Step = [];
    Start_min = [];
    End_min = [];
    Duration_s = [];
    Speed_kmh = [];
    Grade_pct = [];
    VO2kg_mean = [];
    HR_mean_bpm = [];
    HR_start_bpm = [];
    HR_end_bpm = [];
    HR_intra_bpm = [];
    HR_std_bpm = [];
    HR_slope_bpm_s = [];
    HR_slope_bpm_min = [];
    HR_R2 = [];
    HR_cost_bpm_per_kmh = [];
    HR_skewness = [];
    HR_skewness_res = [];
    HR_RMS_residual_bpm = [];
    Nsamples = [];

    stepCount = 0;
    for s = 1:numel(changeIdx) - 1
        i1 = changeIdx(s);
        i2 = changeIdx(s + 1) - 1;

        tStep = tInc(i1:i2);
        hrStep = hrInc(i1:i2);
        vo2Step = vo2Inc(i1:i2);
        gradeStep = gradeInc(i1:i2);
        speedStep = speedInc(i1:i2);

        localMask = isfinite(tStep) & isfinite(hrStep) & isfinite(speedStep);
        tStep = tStep(localMask);
        hrStep = hrStep(localMask);
        vo2Step = vo2Step(localMask);
        gradeStep = gradeStep(localMask);
        speedStep = speedStep(localMask);

        if numel(hrStep) < 3
            continue;
        end

        thisSpeed = median(speedStep, 'omitnan');
        if thisSpeed <= recoveryKmh + tolSpeed
            continue;
        end

        stepCount = stepCount + 1;
        tLocal = tStep - tStep(1);
        polynomial = polyfit(tLocal, hrStep, 1);
        hrFit = polyval(polynomial, tLocal);

        hrMean = mean(hrStep, 'omitnan');
        hrStart = hrStep(1);
        hrEnd = hrStep(end);
        hrIntra = hrEnd - hrStart;
        hrStd = std(hrStep, 0, 'omitnan');

        ssResidual = sum((hrStep - hrFit).^2);
        ssTotal = sum((hrStep - hrMean).^2);
        if ssTotal > 0
            rSquared = 1 - ssResidual / ssTotal;
        else
            rSquared = NaN;
        end

        residual = hrStep - hrFit;
        stepCountValues = {
            stepCount, tStep(1) / 60, tStep(end) / 60, ...
            tStep(end) - tStep(1), thisSpeed, median(gradeStep, 'omitnan'), ...
            mean(vo2Step, 'omitnan'), hrMean, hrStart, hrEnd, hrIntra, hrStd, ...
            polynomial(1), polynomial(1) * 60, rSquared, hrMean / thisSpeed, ...
            localSkewness(hrStep), localSkewness(residual), ...
            sqrt(mean(residual.^2, 'omitnan')), numel(hrStep)};

        Step(end + 1, 1) = stepCountValues{1};
        Start_min(end + 1, 1) = stepCountValues{2};
        End_min(end + 1, 1) = stepCountValues{3};
        Duration_s(end + 1, 1) = stepCountValues{4};
        Speed_kmh(end + 1, 1) = stepCountValues{5};
        Grade_pct(end + 1, 1) = stepCountValues{6};
        VO2kg_mean(end + 1, 1) = stepCountValues{7};
        HR_mean_bpm(end + 1, 1) = stepCountValues{8};
        HR_start_bpm(end + 1, 1) = stepCountValues{9};
        HR_end_bpm(end + 1, 1) = stepCountValues{10};
        HR_intra_bpm(end + 1, 1) = stepCountValues{11};
        HR_std_bpm(end + 1, 1) = stepCountValues{12};
        HR_slope_bpm_s(end + 1, 1) = stepCountValues{13};
        HR_slope_bpm_min(end + 1, 1) = stepCountValues{14};
        HR_R2(end + 1, 1) = stepCountValues{15};
        HR_cost_bpm_per_kmh(end + 1, 1) = stepCountValues{16};
        HR_skewness(end + 1, 1) = stepCountValues{17};
        HR_skewness_res(end + 1, 1) = stepCountValues{18};
        HR_RMS_residual_bpm(end + 1, 1) = stepCountValues{19};
        Nsamples(end + 1, 1) = stepCountValues{20};
    end

    %% Normalize stage values
    speedMaxStep = max(Speed_kmh, [], 'omitnan');
    Speed_rel = nan(size(Speed_kmh));
    HR_rel = nan(size(HR_mean_bpm));
    HR_cost_rel = nan(size(HR_mean_bpm));

    for j = 1:numel(Speed_kmh)
        if isfinite(speedMaxStep) && speedMaxStep > 0
            Speed_rel(j) = Speed_kmh(j) / speedMaxStep;
        end

        if isfinite(hrRestVal) && isfinite(hrMax) && hrMax > hrRestVal
            HR_rel(j) = (HR_mean_bpm(j) - hrRestVal) / (hrMax - hrRestVal);
        elseif isfinite(hrMax) && hrMax > 0
            HR_rel(j) = HR_mean_bpm(j) / hrMax;
        end

        if isfinite(HR_rel(j)) && isfinite(Speed_rel(j)) && Speed_rel(j) > 0
            HR_cost_rel(j) = HR_rel(j) / Speed_rel(j);
        end
    end

    StepTable = table( ...
        Step, Start_min, End_min, Duration_s, Speed_kmh, Speed_rel, Grade_pct, VO2kg_mean, ...
        HR_mean_bpm, HR_start_bpm, HR_end_bpm, HR_intra_bpm, HR_std_bpm, ...
        HR_slope_bpm_s, HR_slope_bpm_min, HR_R2, HR_rel, HR_cost_bpm_per_kmh, HR_cost_rel, ...
        HR_skewness, HR_skewness_res, HR_RMS_residual_bpm, Nsamples);

    Results(kTest).FileName = currentFile;
    Results(kTest).SubjectLabel = subjectID;
    Results(kTest).Age = ageVal;
    Results(kTest).Height_cm = heightVal;
    Results(kTest).Weight_kg = weightVal;
    Results(kTest).HRrest = hrRestVal;
    Results(kTest).VO2max_smooth = vo2MaxSmooth;
    Results(kTest).HRmax = hrMax;
    Results(kTest).Reg_intercept = regIntercept;
    Results(kTest).Reg_coeff_speed = regSpeed;
    Results(kTest).Reg_coeff_hr = regHR;
    Results(kTest).Reg_standard_error = regSE;
    Results(kTest).Reg_R2 = regR2;
    Results(kTest).Reg_N = regN;
    Results(kTest).DatiSegnale = signalTable;
    Results(kTest).TabellaStep = StepTable;
end

if isempty(Results)
    error('No test produced valid results.');
end

%% Save summary workbook
outFile = fullfile(pathName, 'vo2_analysis_results.xlsx');
if exist(outFile, 'file')
    delete(outFile);
end

Summary = table( ...
    (1:numel(Results))', string({Results.SubjectLabel})', string({Results.FileName})', ...
    [Results.Age]', [Results.Height_cm]', [Results.Weight_kg]', [Results.HRrest]', ...
    [Results.VO2max_smooth]', [Results.HRmax]', [Results.Reg_intercept]', ...
    [Results.Reg_coeff_speed]', [Results.Reg_coeff_hr]', [Results.Reg_standard_error]', ...
    [Results.Reg_R2]', [Results.Reg_N]', ...
    'VariableNames', {'TestN', 'Subject', 'FileName', 'Age', 'Height_cm', 'Weight_kg', ...
    'HRrest', 'VO2max_smooth', 'HRmax', 'VO2reg_intercept', 'VO2reg_coeff_speed', ...
    'VO2reg_coeff_HR', 'VO2reg_standard_error', 'VO2reg_R2', 'VO2reg_N'});

writetable(Summary, outFile, 'Sheet', 'Summary');
for kTest = 1:numel(Results)
    writetable(Results(kTest).TabellaStep, outFile, ...
        'Sheet', makeValidSheetName(Results(kTest).SubjectLabel, kTest));
end

%% Build normalized sample table
AllSamplesTable = buildAllSamplesTable(Results);
writetable(AllSamplesTable, outFile, 'Sheet', 'AllSamples_Normalized');

%% Global linear model
GlobalLinearSummary = table();
GlobalLinearCoefficients = table();

if height(AllSamplesTable) >= 10
    formula = 'VO2_rel ~ Speed_rel + HR_rel';
    linearModel = fitlm(AllSamplesTable, formula);
    yTrue = AllSamplesTable.VO2_rel;
    yPred = predict(linearModel, AllSamplesTable);
    [rmse, mae, r2] = regressionMetrics(yTrue, yPred);

    GlobalLinearSummary = table(string(formula), rmse, mae, r2, ...
        linearModel.Rsquared.Ordinary, linearModel.Rsquared.Adjusted, ...
        linearModel.ModelCriterion.AIC, linearModel.ModelCriterion.BIC, ...
        height(AllSamplesTable), ...
        'VariableNames', {'Formula', 'RMSE', 'MAE', 'R2_pred', 'R2_model', ...
        'R2_adjusted', 'AIC', 'BIC', 'N_samples'});

    GlobalLinearCoefficients = table( ...
        string(linearModel.CoefficientNames(:)), ...
        linearModel.Coefficients.Estimate, linearModel.Coefficients.SE, ...
        linearModel.Coefficients.tStat, linearModel.Coefficients.pValue, ...
        'VariableNames', {'Coefficient', 'Estimate', 'SE', 'tStat', 'pValue'});

    writetable(GlobalLinearSummary, outFile, 'Sheet', 'Global_LM_Summary');
    writetable(GlobalLinearCoefficients, outFile, 'Sheet', 'Global_LM_Coefficients');
end

%% Linear mixed-effects model
MixedModelSummary = table();
MixedFixedEffects = table();
MixedRandomEffects = table();

if height(AllSamplesTable) >= 10 && numel(unique(AllSamplesTable.Athlete)) >= 2
    formula = 'VO2_rel ~ Speed_rel + HR_rel + (1 | Athlete)';
    mixedModel = fitlme(AllSamplesTable, formula);

    yTrue = AllSamplesTable.VO2_rel;
    yPredConditional = predict(mixedModel, AllSamplesTable, 'Conditional', true);
    yPredMarginal = predictFixedEffectsVO2(mixedModel, AllSamplesTable);
    [rmseCond, maeCond, r2Cond] = regressionMetrics(yTrue, yPredConditional);
    [rmseMarg, maeMarg, r2Marg] = regressionMetrics(yTrue, yPredMarginal);

    MixedModelSummary = table(string(formula), rmseCond, maeCond, r2Cond, ...
        rmseMarg, maeMarg, r2Marg, mixedModel.ModelCriterion.AIC, ...
        mixedModel.ModelCriterion.BIC, height(AllSamplesTable), ...
        'VariableNames', {'Formula', 'RMSE_conditional', 'MAE_conditional', ...
        'R2_conditional', 'RMSE_marginal', 'MAE_marginal', 'R2_marginal', ...
        'AIC', 'BIC', 'N_samples'});

    MixedFixedEffects = table( ...
        string(mixedModel.CoefficientNames(:)), fixedEffects(mixedModel), ...
        'VariableNames', {'FixedEffect', 'Estimate'});

    try
        MixedRandomEffects = randomEffects(mixedModel);
    catch
        MixedRandomEffects = table();
    end

    writetable(MixedModelSummary, outFile, 'Sheet', 'MixedModel_Summary');
    writetable(MixedFixedEffects, outFile, 'Sheet', 'MixedModel_FixedEffects');
    if istable(MixedRandomEffects)
        writetable(MixedRandomEffects, outFile, 'Sheet', 'MixedModel_RandomEffects');
    end

    %% Leave-one-athlete-out validation
    LOAO_Validation = runLOAOValidation(AllSamplesTable);
    writetable(LOAO_Validation, outFile, 'Sheet', 'LOAO_Validation');

    validRows = isfinite(LOAO_Validation.RMSE);
    if any(validRows)
        LOAO_GlobalSummary = table( ...
            mean(LOAO_Validation.RMSE(validRows), 'omitnan'), ...
            mean(LOAO_Validation.MAE(validRows), 'omitnan'), ...
            mean(LOAO_Validation.R2_pred(validRows), 'omitnan'), ...
            mean(LOAO_Validation.MeanAbsPercentError(validRows), 'omitnan'), ...
            sum(validRows), ...
            'VariableNames', {'Mean_RMSE', 'Mean_MAE', 'Mean_R2_pred', ...
            'Mean_AbsPercentError', 'N_valid_athletes'});
        writetable(LOAO_GlobalSummary, outFile, 'Sheet', 'LOAO_GlobalSummary');
    end
end

%% Residual diagnostics
ResidualDiagnostics = table();
BlandAltmanTable = table();
if exist('mixedModel', 'var')
    yTrue = AllSamplesTable.VO2_rel;
    yPred = predictFixedEffectsVO2(mixedModel, AllSamplesTable);
    modelLabel = "Mixed_fixed";
elseif exist('linearModel', 'var')
    yTrue = AllSamplesTable.VO2_rel;
    yPred = predict(linearModel, AllSamplesTable);
    modelLabel = "GlobalLM";
else
    yTrue = [];
    yPred = [];
    modelLabel = "none";
end

if ~isempty(yTrue)
    residual = yTrue - yPred;
    [rmse, mae, r2] = regressionMetrics(yTrue, yPred);
    ResidualDiagnostics = table(string(modelLabel), mean(residual, 'omitnan'), ...
        std(residual, 'omitnan'), sqrt(mean(residual.^2, 'omitnan')), ...
        localSkewness(residual), min(residual, [], 'omitnan'), ...
        max(residual, [], 'omitnan'), rmse, mae, r2, height(AllSamplesTable), ...
        'VariableNames', {'Model', 'ResidualMean', 'ResidualSTD', 'ResidualRMS', ...
        'ResidualSkewness', 'ResidualMin', 'ResidualMax', 'RMSE', 'MAE', ...
        'R2_pred', 'N_samples'});

    blandMean = (yTrue + yPred) / 2;
    blandDifference = yPred - yTrue;
    bias = mean(blandDifference, 'omitnan');
    sdDifference = std(blandDifference, 'omitnan');
    upperLoA = bias + 1.96 * sdDifference;
    lowerLoA = bias - 1.96 * sdDifference;

    BlandAltmanTable = table(blandMean, yTrue, yPred, blandDifference, ...
        repmat(bias, numel(blandDifference), 1), ...
        repmat(upperLoA, numel(blandDifference), 1), ...
        repmat(lowerLoA, numel(blandDifference), 1), ...
        'VariableNames', {'Mean_true_pred', 'VO2_rel_true', 'VO2_rel_pred', ...
        'Diff_pred_minus_true', 'Bias', 'UpperLoA', 'LowerLoA'});

    writetable(ResidualDiagnostics, outFile, 'Sheet', 'Residual_Diagnostics');
    writetable(BlandAltmanTable, outFile, 'Sheet', 'BlandAltman_VO2');
end

fprintf('\nResults saved to:\n%s\n', outFile);

%% Diagnostic figures
[nRows, nCols] = bestTileGrid(numel(Results));

figure('Name', 'Smoothed VO2 and HR', 'Color', 'w');
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for kTest = 1:numel(Results)
    nexttile;
    D = Results(kTest).DatiSegnale;
    yyaxis left;
    plot(D.time_min, D.VO2kg_smooth, 'LineWidth', 1.5);
    ylabel('VO2/kg [ml/min/kg]');
    yyaxis right;
    plot(D.time_min, D.HR, 'LineWidth', 1.5);
    ylabel('HR [bpm]');
    xlabel('Time [min]');
    title(Results(kTest).SubjectLabel, 'Interpreter', 'none');
    grid on;
end
sgtitle('Smoothed VO2/kg and heart rate');

figure('Name', 'Grade and speed', 'Color', 'w');
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for kTest = 1:numel(Results)
    nexttile;
    D = Results(kTest).DatiSegnale;
    yyaxis left;
    plot(D.time_min, D.Grade_pct, 'LineWidth', 1.5);
    ylabel('Grade [%]');
    yyaxis right;
    plot(D.time_min, D.Speed_kmh, 'LineWidth', 1.5);
    ylabel('Speed [km/h]');
    xlabel('Time [min]');
    title(Results(kTest).SubjectLabel, 'Interpreter', 'none');
    grid on;
end
sgtitle('Grade and speed over time');

%% Workspace output
AllResults = Results;

%% Local functions
function idx = findHeaderIndex(headers, candidateNames, fallbackIdx)
    idx = [];
    headers = lower(strtrim(string(headers)));
    for k = 1:numel(candidateNames)
        key = lower(strtrim(string(candidateNames(k))));
        hit = find(headers == key, 1, 'first');
        if isempty(hit)
            hit = find(contains(headers, key), 1, 'first');
        end
        if ~isempty(hit)
            idx = hit;
            return;
        end
    end
    idx = fallbackIdx;
end

function vec = colToNumeric(col)
    vec = nan(numel(col), 1);
    for i = 1:numel(col)
        value = col{i};
        if isnumeric(value) && isscalar(value)
            vec(i) = value;
        elseif isstring(value) || ischar(value)
            text = strtrim(string(value));
            if text ~= "" && text ~= "--"
                vec(i) = str2double(replace(text, ",", "."));
            end
        end
    end
end

function sec = toSeconds(value)
    sec = NaN;
    if isduration(value)
        sec = seconds(value);
        return;
    end
    if isnumeric(value) && isscalar(value) && isfinite(value)
        if value >= 0 && value < 1
            sec = value * 24 * 3600;
        else
            sec = value;
        end
        return;
    end
    if isstring(value) || ischar(value)
        text = strtrim(string(value));
        if text == "" || text == "--"
            return;
        end
        parts = split(text, ":");
        numbers = str2double(parts);
        if any(isnan(numbers))
            sec = str2double(replace(text, ",", "."));
        elseif numel(numbers) == 2
            sec = numbers(1) * 60 + numbers(2);
        elseif numel(numbers) == 3
            sec = numbers(1) * 3600 + numbers(2) * 60 + numbers(3);
        end
    end
end

function value = cellToNum(cellValue)
    if isnumeric(cellValue) && isscalar(cellValue)
        value = cellValue;
    elseif isstring(cellValue) || ischar(cellValue)
        value = str2double(replace(strtrim(string(cellValue)), ",", "."));
    else
        value = NaN;
    end
end

function value = findMetaNumeric(raw, labels)
    value = NaN;
    if size(raw, 2) < 2
        return;
    end
    firstColumn = lower(strtrim(string(raw(:, 1))));
    for k = 1:numel(labels)
        key = lower(strtrim(string(labels(k))));
        index = find(firstColumn == key, 1, 'first');
        if isempty(index)
            index = find(contains(firstColumn, key), 1, 'first');
        end
        if ~isempty(index)
            value = cellToNum(raw{index, 2});
            return;
        end
    end
end

function value = findMetaString(raw, labels, fallback)
    value = fallback;
    if size(raw, 2) < 2
        return;
    end
    firstColumn = lower(strtrim(string(raw(:, 1))));
    for k = 1:numel(labels)
        key = lower(strtrim(string(labels(k))));
        index = find(firstColumn == key, 1, 'first');
        if isempty(index)
            index = find(contains(firstColumn, key), 1, 'first');
        end
        if ~isempty(index)
            rawValue = raw{index, 2};
            if ischar(rawValue) || isstring(rawValue)
                value = char(strtrim(string(rawValue)));
            elseif isnumeric(rawValue) && isscalar(rawValue)
                value = num2str(rawValue);
            end
            return;
        end
    end
end

function smoothed = movingMeanByTime(timeSeconds, values, halfWindowSeconds)
    smoothed = nan(size(values));
    for i = 1:numel(values)
        if ~isfinite(timeSeconds(i))
            continue;
        end
        mask = isfinite(timeSeconds) & isfinite(values) & ...
            timeSeconds >= timeSeconds(i) - halfWindowSeconds & ...
            timeSeconds <= timeSeconds(i) + halfWindowSeconds;
        if any(mask)
            smoothed(i) = mean(values(mask), 'omitnan');
        end
    end
end

function skew = localSkewness(values)
    values = values(isfinite(values));
    if numel(values) < 3
        skew = NaN;
        return;
    end
    standardDeviation = std(values, 0);
    if standardDeviation == 0
        skew = 0;
    else
        skew = mean(((values - mean(values)) / standardDeviation).^3);
    end
end

function name = stripExtension(filename)
    name = regexprep(filename, '\.(xlsx|xls|csv)$', '', 'ignorecase');
end

function [nRows, nCols] = bestTileGrid(numberOfPlots)
    nCols = max(1, ceil(sqrt(numberOfPlots)));
    nRows = max(1, ceil(numberOfPlots / nCols));
end

function sheetName = makeValidSheetName(baseName, index)
    if isempty(baseName)
        baseName = sprintf('Test_%d', index);
    end
    baseName = regexprep(char(baseName), '[:\\/\?\*\[\]]', '_');
    baseName = strtrim(baseName);
    if numel(baseName) > 25
        baseName = baseName(1:25);
    end
    sheetName = sprintf('%02d_%s', index, baseName);
end

function [intercept, coefficientSpeed, coefficientHR, standardError, rSquared, sampleCount] = ...
        fitVO2LinearRegression(vo2, speed, heartRate)
    mask = isfinite(vo2) & isfinite(speed) & isfinite(heartRate);
    y = vo2(mask);
    xSpeed = speed(mask);
    xHR = heartRate(mask);
    sampleCount = numel(y);
    intercept = NaN;
    coefficientSpeed = NaN;
    coefficientHR = NaN;
    standardError = NaN;
    rSquared = NaN;

    if sampleCount < 4
        return;
    end

    design = [ones(sampleCount, 1), xSpeed, xHR];
    if rank(design) < size(design, 2)
        return;
    end

    beta = design \ y;
    predicted = design * beta;
    residual = y - predicted;
    sse = sum(residual.^2);
    sst = sum((y - mean(y)).^2);
    degreesOfFreedom = sampleCount - size(design, 2);

    if degreesOfFreedom > 0
        standardError = sqrt(sse / degreesOfFreedom);
    end
    if sst > 0
        rSquared = 1 - sse / sst;
    end

    intercept = beta(1);
    coefficientSpeed = beta(2);
    coefficientHR = beta(3);
end

function allSamples = buildAllSamplesTable(results)
    athlete = strings(0, 1);
    timeSeconds = [];
    timeMinutes = [];
    vo2Smooth = [];
    vo2Relative = [];
    heartRate = [];
    heartRateRelative = [];
    speedKmh = [];
    speedRelative = [];
    gradePercent = [];
    age = [];
    height = [];
    weight = [];
    restingHR = [];
    maximumHR = [];
    maximumVO2 = [];

    for k = 1:numel(results)
        signal = results(k).DatiSegnale;
        if isempty(signal) || height(signal) == 0
            continue;
        end

        hrRest = results(k).HRrest;
        hrMax = results(k).HRmax;
        vo2Max = results(k).VO2max_smooth;
        speedMax = max(signal.Speed_kmh, [], 'omitnan');

        vo2Rel = signal.VO2kg_smooth / vo2Max;
        if isfinite(hrRest) && isfinite(hrMax) && hrMax > hrRest
            hrRel = (signal.HR - hrRest) / (hrMax - hrRest);
        else
            hrRel = signal.HR / hrMax;
        end
        speedRel = signal.Speed_kmh / speedMax;

        hrRel = max(0, min(1.5, hrRel));
        speedRel = max(0, min(1.5, speedRel));

        mask = isfinite(signal.time_s) & isfinite(vo2Rel) & isfinite(hrRel) & ...
            isfinite(speedRel) & isfinite(results(k).Age) & ...
            isfinite(results(k).Weight_kg) & isfinite(results(k).Height_cm) & ...
            signal.Speed_kmh > 0;
        n = sum(mask);
        if n == 0
            continue;
        end

        athlete = [athlete; repmat(string(results(k).SubjectLabel), n, 1)];
        timeSeconds = [timeSeconds; signal.time_s(mask)];
        timeMinutes = [timeMinutes; signal.time_min(mask)];
        vo2Smooth = [vo2Smooth; signal.VO2kg_smooth(mask)];
        vo2Relative = [vo2Relative; vo2Rel(mask)];
        heartRate = [heartRate; signal.HR(mask)];
        heartRateRelative = [heartRateRelative; hrRel(mask)];
        speedKmh = [speedKmh; signal.Speed_kmh(mask)];
        speedRelative = [speedRelative; speedRel(mask)];
        gradePercent = [gradePercent; signal.Grade_pct(mask)];
        age = [age; repmat(results(k).Age, n, 1)];
        height = [height; repmat(results(k).Height_cm, n, 1)];
        weight = [weight; repmat(results(k).Weight_kg, n, 1)];
        restingHR = [restingHR; repmat(results(k).HRrest, n, 1)];
        maximumHR = [maximumHR; repmat(results(k).HRmax, n, 1)];
        maximumVO2 = [maximumVO2; repmat(results(k).VO2max_smooth, n, 1)];
    end

    allSamples = table(categorical(athlete), timeSeconds, timeMinutes, vo2Smooth, ...
        vo2Relative, heartRate, heartRateRelative, speedKmh, speedRelative, ...
        gradePercent, age, height, weight, restingHR, maximumHR, maximumVO2, ...
        'VariableNames', {'Athlete', 'time_s', 'time_min', 'VO2kg_smooth', ...
        'VO2_rel', 'HR', 'HR_rel', 'Speed_kmh', 'Speed_rel', 'Grade_pct', ...
        'Age', 'Height_cm', 'Weight_kg', 'HRrest', 'HRmax', 'VO2max_smooth'});
end

function [rmse, mae, rSquared] = regressionMetrics(trueValues, predictedValues)
    mask = isfinite(trueValues) & isfinite(predictedValues);
    trueValues = trueValues(mask);
    predictedValues = predictedValues(mask);
    rmse = NaN;
    mae = NaN;
    rSquared = NaN;
    if isempty(trueValues)
        return;
    end
    errorValues = predictedValues - trueValues;
    rmse = sqrt(mean(errorValues.^2, 'omitnan'));
    mae = mean(abs(errorValues), 'omitnan');
    ssResidual = sum((trueValues - predictedValues).^2);
    ssTotal = sum((trueValues - mean(trueValues)).^2);
    if ssTotal > 0
        rSquared = 1 - ssResidual / ssTotal;
    end
end

function predicted = predictFixedEffectsVO2(model, tableData)
    beta = fixedEffects(model);
    names = string(model.CoefficientNames(:));
    predicted = zeros(height(tableData), 1);
    for i = 1:numel(beta)
        name = names(i);
        if name == "(Intercept)" || name == "Intercept"
            predicted = predicted + beta(i);
        elseif name == "Speed_rel"
            predicted = predicted + beta(i) .* tableData.Speed_rel;
        elseif name == "HR_rel"
            predicted = predicted + beta(i) .* tableData.HR_rel;
        elseif name == "Grade_pct"
            predicted = predicted + beta(i) .* tableData.Grade_pct;
        elseif name == "Age"
            predicted = predicted + beta(i) .* tableData.Age;
        elseif name == "Weight_kg"
            predicted = predicted + beta(i) .* tableData.Weight_kg;
        elseif name == "Height_cm"
            predicted = predicted + beta(i) .* tableData.Height_cm;
        elseif name == "Speed_rel:HR_rel" || name == "HR_rel:Speed_rel"
            predicted = predicted + beta(i) .* (tableData.Speed_rel .* tableData.HR_rel);
        end
    end
end

function validation = runLOAOValidation(allSamples)
    athletes = categories(allSamples.Athlete);
    athleteName = strings(numel(athletes), 1);
    nTrain = nan(numel(athletes), 1);
    nTest = nan(numel(athletes), 1);
    rmse = nan(numel(athletes), 1);
    mae = nan(numel(athletes), 1);
    r2 = nan(numel(athletes), 1);
    meanError = nan(numel(athletes), 1);
    meanAbsolutePercentageError = nan(numel(athletes), 1);
    formulaUsed = strings(numel(athletes), 1);

    for a = 1:numel(athletes)
        testMask = allSamples.Athlete == athletes{a};
        trainMask = ~testMask;
        trainTable = allSamples(trainMask, :);
        testTable = allSamples(testMask, :);
        athleteName(a) = string(athletes{a});
        nTrain(a) = height(trainTable);
        nTest(a) = height(testTable);

        if height(trainTable) < 10 || height(testTable) < 3 || ...
                numel(unique(trainTable.Athlete)) < 2
            continue;
        end

        formula = 'VO2_rel ~ Speed_rel + HR_rel + (1 | Athlete)';
        model = fitlme(trainTable, formula);
        prediction = predictFixedEffectsVO2(model, testTable);
        [rmse(a), mae(a), r2(a)] = regressionMetrics(testTable.VO2_rel, prediction);
        meanError(a) = mean(prediction - testTable.VO2_rel, 'omitnan');

        valid = isfinite(testTable.VO2_rel) & isfinite(prediction) & ...
            abs(testTable.VO2_rel) > 0;
        if any(valid)
            meanAbsolutePercentageError(a) = mean(abs( ...
                (prediction(valid) - testTable.VO2_rel(valid)) ./ ...
                testTable.VO2_rel(valid)) * 100, 'omitnan');
        end
        formulaUsed(a) = string(formula);
    end

    validation = table(athleteName, nTrain, nTest, rmse, mae, r2, meanError, ...
        meanAbsolutePercentageError, formulaUsed, ...
        'VariableNames', {'Athlete', 'N_train', 'N_test', 'RMSE', 'MAE', ...
        'R2_pred', 'MeanError', 'MeanAbsPercentError', 'FormulaUsed'});
end