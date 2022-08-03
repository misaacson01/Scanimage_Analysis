
%currently outputs some simple summary data, e.g.
%number of yellow, red, total active neurons
%number of neurons linked between before/after datasets
%transients-per-minute of every roi in a big matrix: 
%   [before-spontaneous, after-spontaneous, before-stim, after-stim]
%same data as the matrix, but in a scatterplot

%% settings
analysis_folder = 'C:\Users\misaa\Desktop\2022-07-03 09-18 22_5_8_9_ly6g_test';
fps = 3.4;

%enter roi info for all rois here: [roi#before, roi#after, celltype]
%(celltype = 0 for unknown, 1 for yellow, 2 for red)
roi_info = [16 nan 0;...
    10 28 1;...
    19 16 0;...
    18 54 0;...
    9 nan 0;...
    21 nan 0;...
    11h nan 0;...
    6 nan 0;...
    1 3 0;...
    2 10 0;...
    4 nan 0;...
    15 5 0;...
    8 2 0;...
    5 nan 0];
    

%%load files 
%load files for time 1
unmixed_1 = read_file(fullfile(analysis_folder,'1','celltypes_unmixed.tif'));
load(fullfile(analysis_folder,'1','suite2p','suite2p','plane0','Fall.mat'));
F_1 = F;
spks_1 = spks;
Fneu_1 = Fneu;
stat_1 = stat;
iscell_1 = iscell;
cells50 = find(iscell_1(:,1));
num_rois = length(cells50);
colors = linspecer(num_rois);
o = randperm(num_rois);
colors_1 = colors(o,:);
        
%load files for time 2
unmixed_2 = read_file(fullfile(analysis_folder,'2','celltypes_unmixed.tif'));
load(fullfile(analysis_folder,'2','suite2p','suite2p','plane0','Fall.mat'));
F_2 = F;
spks_2 = spks;
Fneu_2 = Fneu;
stat_2 = stat;
iscell_2 = iscell;
cells50 = find(iscell_2(:,1));
num_rois = length(cells50);
colors = linspecer(num_rois);
o = randperm(num_rois);
colors_2 = colors(o,:);

load(fullfile(analysis_folder,'exp_metadata.mat'));

%for using all cells from suite2p
num_rois_1 = sum(iscell_1(:,1)==1);
num_rois_2 = sum(iscell_2(:,1)==1);
roi_info = nan(max([num_rois_1 num_rois_2]),3);
roi_info(1:num_rois_1,1) = sort(find(iscell_1(:,1)));
roi_info(1:num_rois_2,2) = sort(find(iscell_2(:,1)));
roi_info(:,3) = 0;

%check for duplicates
for t = 1:2
    tmp = roi_info(:,t);
    tmp(isnan(tmp)) = [];
    assert(length(tmp)==length(unique(tmp)),['redundant roi(s) found in roi_info (column ' num2str(t) ')']);
end


%% analyze data
%summary data of roi info
num_active = sum(~isnan(roi_info(:,1:2)))
num_active_yellow = sum(~isnan(roi_info(:,1:2)) & repmat(roi_info(:,3)==1,[1 2]))
num_active_red = sum(~isnan(roi_info(:,1:2)) & repmat(roi_info(:,3)==2,[1 2]))
num_linked = sum(sum(~isnan(roi_info(:,1:2)),2)==2)
num_linked_yellow = sum(sum(~isnan(roi_info(:,1:2)),2)==2 & roi_info(:,3)==1)
num_linked_red = sum(sum(~isnan(roi_info(:,1:2)),2)==2 & roi_info(:,3)==2)

num_rois = size(roi_info,1);
%pre-allocate variables
roi_spikespermin = nan(num_rois,4);
roi_transientspermin = nan(num_rois,4);
roi_kontransientspermin = nan(num_rois,2);
numRepetitions = vstruct(1).num_reps;
repetitionsToAverage = 2:numRepetitions;
stimOnFrames = 7:13;
numConditions = vstruct(1).num_trials;
numFramesPerTrial = mdata(t).num_stim_frames;
stimStart = 1+mdata(1).num_spont_frames;
spontMinutes = mdata(1).num_spont_frames/(fps*60);
stimMinutes = (mdata(t).num_stim_frames*mdata(t).num_stims)/(fps*60);
roi_stimspks = nan([numRepetitions, numConditions, mdata(1).num_stim_frames]);
roi_stimtransients = nan([numRepetitions, numConditions, mdata(1).num_stim_frames]);
roi_DSI = nan(num_rois,2);
roi_OSI = nan(num_rois,2);
for t = 1:2
    %get frame indices for each trial, organized by [repetition, conditions, frame]
    stimFrames = vstruct(t).order-1; %each row is 0:numConditions-1
    stimFrames = stimFrames + repmat(numConditions*(0:numRepetitions-1)',[1 numConditions]); %value for trial1 rep2 is now numConds+1, eg
    stimFrames = stimFrames*numFramesPerTrial; %values are now frame numbers, relative to the first trial of each repetition
    stimFrames = stimFrames + stimStart; %values are the actual frame numbers when each trial started
    stimFrames = repmat(stimFrames, [1 1 numFramesPerTrial]);
    stimFrames = stimFrames + repmat(permute((0:numFramesPerTrial-1),[1 3 2]),[numRepetitions, numConditions, 1]);

    if t==1
        spks = spks_1;
        F = F_1;
    else
        spks = spks_2;
        F = F_2;
    end
    rois = roi_info(:,t);
    for r = 1:num_rois
        roi = rois(r);
        if ~isnan(roi)
            %spikes/min
            roi_spikespermin(r,t) = sum(spks(roi,1:mdata(t).num_spont_frames))/spontMinutes; %spontaneous 
            roi_spikespermin(r,2+t) = sum(spks(roi,1+mdata(t).num_spont_frames:end))/stimMinutes; %stim-evoked
            
            %transients/min
            roi_transientspermin(r,t) = sum(spks(roi,1:mdata(t).num_spont_frames)>0)/spontMinutes; %spontaneous 
            roi_transientspermin(r,2+t) = sum(spks(roi,1+mdata(t).num_spont_frames:end)>0)/stimMinutes; %stim-evoked 
            
            %kontransients/min
            tmp = F(roi,1:mdata(t).num_spont_frames);
            tmpmean = prctile(tmp,20);
            tmpstd = std(tmp);
            tmp = tmp>(tmpmean+3*tmpstd);
            tmp = diff(tmp)==1;
            roi_kontransientspermin(r,t) = sum(tmp)/spontMinutes;
            
            %organize spikes, F, Fneu by [repetition, trial, frame, roi]
            roi_stimspks = spks(sub2ind(size(F),r*ones(size(stimFrames)),stimFrames));
            roi_stimtransients = spks(sub2ind(size(F),r*ones(size(stimFrames)),stimFrames))>0;
            
            %average spike data for repetitions
            sumStimSpks = sum(roi_stimspks(:,:,stimOnFrames),3); %sum the spikes from 0-2s
            meanSumStimSpks = mean(sumStimSpks(repetitionsToAverage,:),1); %average repetitions
            
            %OSI (spks)
            tmpSpksO = mean([meanSumStimSpks(1:6); meanSumStimSpks(7:12)]);
            [maxSpk,maxI] = max(tmpSpksO);
            oppI = 1+mod((maxI+3)-1,6);
            oppSpk = tmpSpksO(oppI);
            roi_OSI(r,t) = (maxSpk-oppSpk)/(maxSpk+oppSpk);
            
            %DSI (spks)
            [maxSpk,maxI] = max(meanSumStimSpks(1:12)); %condition which cause the most spikes
            oppI = 1+mod((maxI+6)-1,12); %condition on the opposite angle
            oppSpk = meanSumStimSpks(oppI);
            roi_DSI(r,t) = (maxSpk-oppSpk)/(maxSpk+oppSpk);
        end
    end
end
%%plot before/after silent/normal/hyperactive
bins = zeros(2,3); %[time, bin]
for t = 1:2
    spm = roi_kontransientspermin(:,t);
    spm(isnan(spm)) = [];
    bins(t,1) = bins(t,1) + sum(spm<=0.25);
    bins(t,2) = bins(t,2) + sum(spm>0.25 & spm<=4);
    bins(t,3) = bins(t,3) + sum(spm>4);
end

%%
colors = linspecer(2);
plot_size = [1000 250];
figure('Position',[100 100 plot_size]);
subplot(1,3,1);
y = roi_kontransientspermin(:,2);
y(isnan(y)) = [];
y = y';
h = histogram(y,(0:10));
h.FaceColor = colors(2,:);
h.FaceAlpha = 0.6;
hold on
y = roi_kontransientspermin(:,1);
y(isnan(y)) = [];
y = y';
h = histogram(y,(0:10));
h.FaceColor = colors(1,:);
h.FaceAlpha = 0.6;
hold on
ylim([0 48])
xlabel('transients/min');
ylabel('# of ROIs');
title('spontaneous activity distribution');

subplot(1,3,2);
% x = [1 2];
y = bins';
bar(y);
% b = gca;
% b.FaceColor = 'flat';
% b(1).CData = repmat(colors(1,:),[3 1]);
% b(2).CData = repmat(colors(2,:),[3 1]);
xlim([0 4])
ylim([0 48]);
xticklabels({'silent','normal','hyperactive'})
ylabel('# ROIs');
% ax = gca;
% ax.Position = [0.13 0.11 0.12 0.815];
title('activity classification')


%% plot before/after cumulative distribution plots
colors = linspecer(2);
plot_size = [1000 250];
figure('Position',[100 100 plot_size]);
subplot(1,4,1);
x = [1 2];
y = num_active;
b = bar(x,y,0.5);
b.FaceColor = 'flat';
b.CData = colors;
xlim([0.5 2.5])
ylim([0 150]);
xticklabels({'baseline','treated'})
ylabel('# active ROIs');
ax = gca;
ax.Position = [0.13 0.11 0.12 0.815];
title('active cells')
%spontaneous activity (transients per min)
subplot(1,4,2);
max_tpm = max(roi_transientspermin,[],'all');
x = 0:0.1:max_tpm;
for d = 1:2
    y = zeros(size(x));
    roi_inds = find(~isnan(roi_info(:,d)));
    yvals = floor(10*sort(roi_transientspermin(roi_inds,d)));
    for i = yvals
        y(i+1) = y(i+1)+1;
    end
    y = y/10;
    y = cumsum(y);
    y = y/max(y);
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
% legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('transients/min')
ylabel('cumulative probability')
xlim([30 50])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
title('spontaneous activity')

%spontaneous activity (spikes per min)
% subplot(1,4,2);
% max_tpm = max(roi_spikespermin,[],'all');
% x = 0:max_tpm;
% for d = 1:2
%     y = zeros(size(x));
%     roi_inds = find(~isnan(roi_info(:,d)));
%     yvals = floor(sort(roi_spikespermin(roi_inds,d)/1000));
%     for i = yvals
%         y(i+1) = y(i+1)+1;
%     end
%     y = cumsum(y);
%     y = y/max(y);
%     plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
%     hold on
% end
% % legend({'baseline','treated'},'Location','NorthEastOutside')
% xlabel('spikes/min')
% ylabel('cumulative probability')
% xlim([0 200])
% ylim([0 1])
% yticks([0 0.2 0.4 0.6 0.8 1])
% title('spontaneous activity')

%orientation selectivity
% figure('Position',[100 100 plot_size]);
subplot(1,4,3)
x = 0:0.01:1;
for d = 1:2
    y = zeros(1,101);
    roi_inds = find(~isnan(roi_info(:,d)));
    yvals = floor(100*sort(roi_OSI(roi_inds,d)));
    for i = yvals
        y(i+1) = y(i+1)+1;
    end
    y = y/100;
    y = cumsum(y);
    y = y/max(y);
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
% legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('OSI')
% ylabel('cumulative probability')
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.2 0.4 0.6 0.8 1])
title('orientation selectivity')

%direction selectivity
subplot(1,4,4)
% figure('Position',[100 100 plot_size]);
x = 0:0.01:1;
for d = 1:2
    y = zeros(1,101);
    roi_inds = find(~isnan(roi_info(:,d)));
    yvals = floor(100*sort(roi_DSI(roi_inds,d)));
    for i = yvals
        y(i+1) = y(i+1)+1;
    end
    y = y/100;
    y = cumsum(y);
    y = y/max(y);
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
% legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('DSI')
% ylabel('cumulative probability')
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.2 0.4 0.6 0.8 1])
title('direction selectivity')

