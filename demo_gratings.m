%% Demo code to visualise responses to drifting gratings

clear;

addpath('\\zserver.cortexlab.net\Code\2photonPipeline')
addpath('\\zserver.cortexlab.net\Code\Matteobox')
addpath('\\zserver.cortexlab.net\Code\Neuropil Correction');
make_db_FR201; 

iExp = 1;
iPlane = db(iExp).plane;

%%

% function that returns the directory (root) and processed data file
% (refF) from s2p. Edit the function with your local data directory
[root, refF, ~] = getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID), iPlane); 

if exist(fullfile(root, refF), 'file')
    
    load(fullfile(root, refF), 'dat');
    isCell = logical([dat.stat.iscell]);
    red = [dat.stat.redcell];
    isRed = red(isCell);
    
    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts);
    Fcell = dat.Fcell{db(iExp).expID}(isCell, :); % nN*nT
    Fneu = dat.FcellNeu{db(iExp).expID}(isCell, :);% nN*nT
    dF = estimateNeuropil_LFR(Fcell,Fneu); % correct neuropil and then subtracts 5th prctile of trace to equalize eventual drifts across exps
    
end

%%

% smooth and zscore the F signals
dF = gaussFilt(dF', 1)';
dF = zscore(dF');

%timestamp fluorescence signals
nFrames = numel(dF);
planeFrames = db(iExp).plane:info.nPlanes:(nFrames*info.nPlanes);
frameTimes = ppbox.getFrameTimes(info, planeFrames);
frameRate = (1/mean(diff(frameTimes)));

%load and timestamp stimulus matrix
stimTimes = ppbox.getStimTimes(info);
stimSequence = ppbox.getStimSequence_LFR(info);
stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);

if isfield(db(iExp), 'stimList') 
    stimSet = db(iExp).stimList;
    stimLabels{iExp} = stimSequence.labels(stimSet);
else
    stimSet = 1:numel(stimSequence.labels);
    stimLabels{iExp} = stimSequence.labels(stimSet);
end

[responses, aveResp, seResp, kernelTime, stimDur] = ...
    getStimulusSweepsLFR(dF, stimTimes, stimMatrix,frameRate); % responses is (nroi, nStim, nResp, nT)

respWin = [0 3];

[resPeak, aveResPeak, seResPeak] = ...
    gratingOnResp(responses, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

plotSweepResp(responses, kernelTime, stimDur);


%% fit tuning curve

[~,  nStim, nRep, ~] = size(responses);

aveResp = nanmean(responses, 3);
seResp = nanstd(responses, 1,3)/sqrt(nRep);
aveResPeak = mean(resPeak, 3);
seResPeak = std(resPeak, 1, 3)/sqrt(nRep);

figure;
oris = 0:30:330;
for iStim = 1: nStim-1
    
    subplot(1, nStim+2 , iStim)
    
    shadePlot(kernelTime, squeeze(aveResp(1,iStim,:, :)), squeeze(seResp(1,iStim,:, :)), [0 0 0])
    
    xlim([min(kernelTime),max(kernelTime)])
    ylim([min(aveResp(:)), max(responses(:))])
    ylim([-0.5 5])
    set(gca, 'Xtick', [], 'YTick', [],'visible', 'off')
    formatAxes
end
hold on
plot([-1, -1], [1, 2] , 'k')

subplot(1, nStim+2, [nStim: nStim+2])

toFit= makeVec(resPeak(1, 1:end-1, :))';
[tunePars, ~] = fitori(repmat(oris, 1,nRep), toFit);%, [], [NaN, NaN, 0, NaN, NaN]);
% [tunePars, ~] = fitori(oris, makeVec(aveAllResPeak(1, 1:end-1))');%, [], [NaN, NaN, 0, NaN, NaN]);
oriFit = oritune(tunePars, 0:5:330);
errorbar(0:30:360, aveResPeak, seResPeak, 'ok'); hold on
plot(0:5:330, oriFit, '-r', 'LineWidth', 2)
formatAxes
set(gca, 'YColor', 'none','YTick', [], 'XTick', [0:90:360], 'XTickLabel', {0:90:270, 'Blank'});
xlim([-10 , 370])
ylim([min(aveResPeak) - mean(seResPeak), max(aveResPeak) + max(seResPeak)])

try
    title([db(1).mouse_name(1:7), ' ', db(1).mouse_name(9:end)])
catch
    title([db(1).mouse_name])
end

figure;
oriFit = oritune(tunePars, 0:5:360);
polarplot(pi*(0:5:360)/(180), (-min(oriFit)+ oriFit)/(max(oriFit)-min(oriFit)))
formatAxes