%% RAMPKAR PAPER
% Figure 2 - Cross-correlation RAMPKAR2-HYlight and RAMPKAR2-Perceval

% basic paths
addpath('\\albecklab.mcb.ucdavis.edu\data\Code\DatalocHandler\', ...
    '\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis\', ...
    '\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace\')

iBasePath = '\\albecklab.mcb.ucdavis.edu\data\imageData\RAMPKAR\';
pBasePath = '\\albecklab.mcb.ucdavis.edu\data\Processed Data\RAMPKAR\';

%% Load the data for RAMPKAR-HYlight for crosscorr

% List of replicate subpaths (relative to pBasePath)
replicatePaths = {
    '2023-09-26 10A RAMPKAR2 HYLIGHT\2023-09-26 10A RAMPKAR2 HYLIGHT_Processed.mat' % has MTK + oligo  treatment
    '2023-09-25 10A RAMPKAR2 HYLIGHT\2023-09-25 10A RAMPKAR2 HYLIGHT_Processed.mat' % has MTK + oligo  treatment
    '2025-06-10_MCF10A_R2HmO_R2P\2025-06-10_MCF10A_R2HmO_R2P_Processed.mat'% has MTK + oligo  treatment
    '2025-05-22_MCF10A_R2HmO_R2P\2025-05-22_MCF10A_R2HmO_R2P_Processed.mat' % this one is behaving very weirdly after pwrat
    '2025-05-14_MCF10A_R2HmO_R2P\2025-05-14_MCF10A_R2HmO_R2P_Processed.mat'};

% Preallocate
f2Data = cell(1, numel(replicatePaths));

% Load in a loop
for i = 1:numel(replicatePaths)
    f2D = load([pBasePath, replicatePaths{i}]);
    f2Data{i} = f2D.dataloc;
end

%% Get pulses and subset the data, look at the lag values
% Common pattern to remove from treatment names and create simplified names
pat = optionalPattern(" ") + "at hour " + optionalPattern("+"|"-") + digitsPattern(1) + optionalPattern(" and");

%% Loop for Replicates 1–2 (Markhus')

for i = 1:2
    deltaR2_HY = convertPulseToDataframe(f2Data{i}, {'HYLIGHT','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 3, ...
        'crosscorr', {'HYLIGHT_RAMPKAR2'}, ...
        'exclude', {'IPA'}, ...
        'tmaxaftertx', 12);

    deltaR2_HY.tx = erase(deltaR2_HY.treatment, pat);
    deltaR2_HY.tx = erase(deltaR2_HY.tx, '1 INSIM ');
    deltaR2all_HY{i} = deltaR2_HY;

end

%% Loop for Replicates 3–5 (Marion's)

for i = 3:5
    deltaR2_HY = convertPulseToDataframe(f2Data{i}, {'HYLIGHT','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 5, ...
        'crosscorr', {'HYLIGHT_RAMPKAR2'}, ...
        'tmaxaftertx', 12);

    deltaR2_HY.tx = erase(deltaR2_HY.treatment, pat);
    deltaR2all_HY{i} = deltaR2_HY;
end

for i = [3 5] % currently, pulsanalysis errors out on RAMPKAR2 for x = 4, maybe because of data quality ??
    deltaR2_Per = convertPulseToDataframe(f2Data{i}, {'PERCEVAL','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 5, ...
        'crosscorr', {'PERCEVAL_RAMPKAR2'}, ...
        'tmaxaftertx', 12);

    deltaR2_Per.tx = erase(deltaR2_Per.treatment, pat);
    deltaR2all_Per{i} = deltaR2_Per;

end
%% Pool data
% Concatenate the replicates
dataH = vertcat(deltaR2all_HY{:});
dataP = vertcat(deltaR2all_Per{:});

%% Check out the distribution of x_corr values and lags in all conditions
% HYLIGHT vs RAMPKAR2

hCorr = deltaR2_Per(abs(deltaR2_Per.PERCEVAL_RAMPKAR2_maxXcorr)>0.5,:);
highCorrs = grpstats(hCorr,"tx",["mean","median","sem"],"DataVars",["PERCEVAL_RAMPKAR2_maxXcorrLagInMin","PERCEVAL_RAMPKAR2_maxXcorr"])

binRange = [-150, 150];
edges = linspace(binRange(1), binRange(2), 100); % 100 bins means 51 edges

figure;
subplot(2,2,1)
histogram(dataH.HYLIGHT_RAMPKAR2_maxXcorr,100, Normalization='count');
hold on;
histogram(dataH.HYLIGHT_Scrambled_RAMPKAR2_maxXcorr,100, Normalization='count');
alpha(.4);
xlabel('max X corr');
ylabel('count (HYlight)');
subplot(2,2,2)
histogram(dataH.HYLIGHT_RAMPKAR2_maxXcorrLagInMin,edges, Normalization='count');
hold on;
histogram(dataH.HYLIGHT_Scrambled_RAMPKAR2_maxXcorrLagInMin,edges, Normalization='count');
alpha(.4);
xlim([-150 150]);
xlabel('max X corr lag (min)');
ylabel('count (HYlight)');
legend('actual','scrambled')

% Perceval vs RAMPKAR2

subplot(2,2,3)
histogram(dataP.PERCEVAL_RAMPKAR2_maxXcorr,100, Normalization='count');
hold on;
histogram(dataP.PERCEVAL_Scrambled_RAMPKAR2_maxXcorr,100, Normalization='count');
alpha(.4);
xlabel('max X corr');
ylabel('count (Perceval)')

subplot(2,2,4)
histogram(dataP.PERCEVAL_RAMPKAR2_maxXcorrLagInMin, edges, Normalization='count');
hold on;
histogram(dataP.PERCEVAL_Scrambled_RAMPKAR2_maxXcorrLagInMin, edges, Normalization='count');
alpha(.4);
xlim([-150 150]);
xlabel('max X corr lag (min)');
ylabel('count (Perceval)');
legend('actual','scrambled')
sgtitle('Sensor vs RAMPKAR2: Whole timeseries');

% Ok optimal lag is at 0 min for both Perceval and HYlight

%% Force a lag of 0 minutes for Hylight
% Markhus' datasets
for i = 1:2
    deltaR2_HY = convertPulseToDataframe(f2Data{i}, {'HYLIGHT','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 3, ...
        'crosscorr', {'HYLIGHT_RAMPKAR2'}, ...
        'exclude', {'IPA','NOINSIM'}, ...
        'tmaxaftertx', 12, 'forcecorrlag',0);

    deltaR2_HY.tx = erase(deltaR2_HY.treatment, pat);
    deltaR2_HY.tx = erase(deltaR2_HY.tx, '1 INSIM ');

    keepTx = {'1 Vehicle','3uM Oligo','10uM MK8722 3uM Oligo'};
    deltaR2Hsub = deltaR2_HY(matches(deltaR2_HY.tx, keepTx), :);
    deltaR2Hsub.tx = categorical(deltaR2Hsub.tx, keepTx);

    deltaR2all_HY{i} = deltaR2_HY;
    deltaR2Hsub_all{i} = deltaR2Hsub;
end

% My datasets, which treatments to keep in subset of data

keepTx = {'17mM Gluc 1 Vehicle 1 Vehicle', ...
    '17mM Gluc 3uM Oligo', ...
    '17mM Gluc 3uM Oligo 10uM MK8722'};

for i = 3:5
    deltaR2_HY = convertPulseToDataframe(f2Data{i}, {'HYLIGHT','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 5, ...
        'crosscorr', {'HYLIGHT_RAMPKAR2'}, ...
        'tmaxaftertx', 12, 'forcecorrlag',0);

    deltaR2_HY.tx = erase(deltaR2_HY.treatment, pat);
    deltaR2Hsub = deltaR2_HY(matches(deltaR2_HY.tx, keepTx), :);

    % Store results
    deltaR2all_HY{i} = deltaR2_HY;
    deltaR2Hsub_all{i} = deltaR2Hsub;
end

for i = [3 5] % currently, pulsanalysis errors out on RAMPKAR2 for x = 4, maybe because of data quality ??
    deltaR2_Per = convertPulseToDataframe(f2Data{i}, {'PERCEVAL','RAMPKAR2'}, ...
        'pulsepars', {'Delta_Mean','Mean_Before'}, ...
        'aftertx', 1, 'tmaxback', 5, ...
        'crosscorr', {'PERCEVAL_RAMPKAR2'}, ...
        'tmaxaftertx', 12, 'forcecorrlag',0);

    deltaR2_Per.tx = erase(deltaR2_Per.treatment, pat);
    deltaR2Psub = deltaR2_Per(matches(deltaR2_Per.tx, keepTx), :);
    deltaR2Psub_all{i} = deltaR2Psub;
    deltaR2all_Per{i} = deltaR2_Per;
end

%% Pool data
% Concatenate the replicates
dataHlag = vertcat(deltaR2all_HY{:});
dataPlag = vertcat(deltaR2all_Per{:});

%% Plot histograms again

% HYLIGHT vs RAMPKAR2

binRange = [-150, 150];
edges = linspace(binRange(1), binRange(2), 100); % 100 bins means 51 edges

figure;
subplot(2,2,1)
histogram(dataHlag.HYLIGHT_RAMPKAR2_maxXcorr,100, Normalization='count');
hold on;
histogram(dataHlag.HYLIGHT_Scrambled_RAMPKAR2_maxXcorr,100, Normalization='count');
alpha(.4);
xlabel('max X corr');
ylabel('count (HYlight)');
subplot(2,2,2)
histogram(dataHlag.HYLIGHT_RAMPKAR2_maxXcorrLagInMin,edges, Normalization='count');
hold on;
histogram(dataHlag.HYLIGHT_Scrambled_RAMPKAR2_maxXcorrLagInMin,edges, Normalization='count');
alpha(.4);
xlim([-150 150]);
xlabel('max X corr lag (min)');
ylabel('count (HYlight)');
legend('actual','scrambled')

% Perceval vs RAMPKAR2

subplot(2,2,3)
histogram(dataPlag.PERCEVAL_RAMPKAR2_maxXcorr,100, Normalization='count');
hold on;
histogram(dataPlag.PERCEVAL_Scrambled_RAMPKAR2_maxXcorr,100, Normalization='count');
alpha(.4);
xlabel('max X corr');
ylabel('count (Perceval)')

subplot(2,2,4)
histogram(dataPlag.PERCEVAL_RAMPKAR2_maxXcorrLagInMin, edges, Normalization='count');
hold on;
histogram(dataPlag.PERCEVAL_Scrambled_RAMPKAR2_maxXcorrLagInMin, edges, Normalization='count');
alpha(.4);
xlim([-150 150]);
xlabel('max X corr lag (min)');
ylabel('count (Perceval)');
legend('actual','scrambled')
sgtitle('Sensor vs RAMPKAR2: Whole timeseries, forced lag 0 min');


%% Violin plots of correlation at optimal lag
% concatenate the replicates
% CHECK IF THIS IS FORCED LAG DATA
dataHsub = vertcat(deltaR2Hsub_all{:}); 
dataPsub = vertcat(deltaR2Psub_all{:}); 

% clean up legend
dataHsub.tx(dataHsub.tx == "17mM Gluc 3uM Oligo") = "3uM Oligo"; % make my annotations match Markhus'
dataHsub.tx(dataHsub.tx == "17mM Gluc 1 Vehicle 1 Vehicle") = "Vehicle";
dataHsub.tx(dataHsub.tx == "1 Vehicle") = "Vehicle";
dataHsub.tx(dataHsub.tx == "17mM Gluc 3uM Oligo 10uM MK8722") = "10uM MK8722 3uM Oligo";
% clean up legend
dataPsub.tx(dataPsub.tx == "17mM Gluc 3uM Oligo") = "3uM Oligo";
dataPsub.tx(dataPsub.tx == "17mM Gluc 1 Vehicle 1 Vehicle") = "Vehicle";
dataPsub.tx(dataPsub.tx == "1 Vehicle") = "Vehicle";
dataPsub.tx(dataPsub.tx == "17mM Gluc 3uM Oligo 10uM MK8722") = "10uM MK8722 3uM Oligo";
%order
dataHsub.tx = categorical(dataHsub.tx,{'Vehicle','3uM Oligo','10uM MK8722 3uM Oligo'});
dataPsub.tx = categorical(dataPsub.tx,{'Vehicle','3uM Oligo','10uM MK8722 3uM Oligo'});


%% Plot

figure;
subplot(2,1,1);
violinplot(dataHsub.HYLIGHT_RAMPKAR2_maxXcorr,dataHsub.tx,'ShowData',true,'ViolinAlpha',0.08, 'EdgeColor',[0.2,0.2,0.2])
ylim([-1,1]); ylabel('Corr AMPK vs FBP')
nCellsH = numel(unique(dataHsub.cellID));
text(0.85, 1, sprintf('%d cells total', nCellsH), 'Units','normalized', 'FontSize',9, 'VerticalAlignment','top');
subplot(2,1,2);
violinplot(dataPsub.PERCEVAL_RAMPKAR2_maxXcorr,dataPsub.tx,'ShowData',true,'ViolinAlpha',0.08, 'EdgeColor',[0.2,0.2,0.2])
ylim([-1,1]); ylabel('Corr AMPK vs ATP/ADP')
nCellsH = numel(unique(dataPsub.cellID));
text(0.85, 1, sprintf('%d cells total', nCellsH), 'Units','normalized', 'FontSize',9, 'VerticalAlignment','top');





