
%currently outputs some simple summary data, e.g.
%number of yellow, red, total active neurons
%number of neurons linked between before/after datasets
%transients-per-minute of every roi in a big matrix: 
%   [before-spontaneous, after-spontaneous, before-stim, after-stim]
%same data as the matrix, but in a scatterplot

%% settings
analysis_folder = 'C:\Users\misaa\Desktop\2022-05-15 07-19 ly6g_test1';
fps = 3.4;

%enter roi info for all rois here: [roi#, celltype]
%(celltype = 0 for unknown, 1 for yellow, 2 for red)
roi_info = [16 0;...
    10 1;...
    19 0;...
    18 0;...
    9 0;...
    21 0;...
    11 0;...
    6 0;...
    1 0;...
    2 0;...
    4 0;...
    15 0;...
    8 0;...
    5 0];
    

%% load files 
%load files for time 1
load(fullfile(analysis_folder,'suite2p','suite2p','plane0','Fall.mat'));
load(fullfile(analysis_folder,'exp_metadata.mat'));
before_spont_frames = 1:1021;
before_stim_frames = 1:1638;
before_stim_frames = before_stim_frames+before_spont_frames(end);
after_spont_frames = 1:1021;
after_spont_frames = after_spont_frames+before_stim_frames(end);
after_stim_frames = 1:1638;
after_stim_frames = after_stim_frames+after_spont_frames(end);


F_1 = F(:,[before_spont_frames before_stim_frames]);
spks_1 = spks(:,[before_spont_frames before_stim_frames]);
iscell_1 = iscell;
cells50 = find(iscell_1(:,1));
num_rois = length(cells50);
colors = linspecer(num_rois);
o = randperm(num_rois);
colors_1 = colors(o,:);
        
%load files for time 2
F_2 = F(:,[after_spont_frames after_stim_frames]);
spks_2 = spks(:,[after_spont_frames after_stim_frames]);
iscell_2 = iscell;
cells50 = find(iscell_2(:,2));
% num_rois = length(cells50);
% colors = linspecer(num_rois);
% o = randperm(num_rois);
colors_2 = colors(o,:);

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
numRepetitions = 6;
repetitionsToAverage = 2:numRepetitions;
stimOnFrames = 7:13;
numConditions = 13;
numFramesPerTrial = 21;
numSpontFrames = 1021;
stimStart = 1022;
spontMinutes = 1021/(fps*60);
roi_stimspks = nan([numRepetitions, numConditions, numFramesPerTrial]);
roi_stimtransients = nan([numRepetitions, numConditions, numFramesPerTrial]);
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
            roi_spikespermin(r,t) = sum(spks(roi,1:numSpontFrames))/spontMinutes; %spontaneous 
            roi_spikespermin(r,2+t) = sum(spks(roi,1+numSpontFrames:end))/spontMinutes; %stim-evoked
            
            %transients/min
            roi_transientspermin(r,t) = sum(spks(roi,1:numSpontFrames)>0)/spontMinutes; %spontaneous 
            roi_transientspermin(r,2+t) = sum(spks(roi,1+numSpontFrames:end)>0)/spontMinutes; %stim-evoked 
            
            %kontransients/min
            tmp = F(roi,1:numSpontFrames);
            tmpmean = prctile(tmp,10);
            tmpstd = std(tmp);
            tmp = tmp>(tmpmean+2.5*tmpstd);
            tmp = diff(tmp)==1;
            roi_kontransientspermin(r,t) = sum(tmp)/spontMinutes;
            
            %organize spikes, F, Fneu by [repetition, trial, frame, roi]
            roi_stimspks = spks(sub2ind(size(spks),r*ones(size(stimFrames)),stimFrames));
            roi_stimtransients = spks(sub2ind(size(spks),r*ones(size(stimFrames)),stimFrames))>0;
            
            %average spike data for repetitions
            sumStimSpks = sum(roi_stimtransients(:,:,stimOnFrames),3); %sum the spikes from 0-2s
            meanSumStimSpks = mean(sumStimSpks(repetitionsToAverage,:),1); %average repetitions
            
            %OSI (spks)
            tmpSpksO = mean([meanSumStimSpks(1:6); meanSumStimSpks(7:12)]);
            [maxSpk,maxI] = max(tmpSpksO);
            oppI = 1+mod((maxI+3)-1,6);
            oppSpk = tmpSpksO(oppI);
            roi_OSI(r,t) = (maxSpk-oppSpk)/(maxSpk+oppSpk);
            if t==1 && r==57 %r=2,21,57
                best_OSI = meanSumStimSpks(1:12);
                best_OSI_val = roi_OSI(r,t)
            end
            
            %DSI (spks)
            [maxSpk,maxI] = max(meanSumStimSpks(1:12)); %condition which cause the most spikes
            oppI = 1+mod((maxI+6)-1,12); %condition on the opposite angle
            oppSpk = meanSumStimSpks(oppI);
            roi_DSI(r,t) = (maxSpk-oppSpk)/(maxSpk+oppSpk);
            if t==1 && r==2 %t=1 r=2,3,7,10; t=2 
                best_DSI = meanSumStimSpks(1:12);
                best_DSI_val = roi_DSI(r,t)
            end
        end
    end
end
%%plot before/after silent/normal/hyperactive
bins = zeros(2,3); %[time, bin]
% bins(1,1) = 45; %# missing in time 1
% bins(2,1) = 7; %# missing in time 2
for t = 1:2
    spm = roi_kontransientspermin(:,t);
    spm(isnan(spm)) = [];
    bins(t,1) = bins(t,1) + sum(spm<0.25);
    bins(t,2) = bins(t,2) + sum(spm>=0.25 & spm<4);
    bins(t,3) = bins(t,3) + sum(spm>=4);
end
bins

%% create polar plots for OSI, DSI
figure_size = [300 250];
figure('Position',[100 100 figure_size]);
theta = deg2rad(0:30:360);
rho = [best_OSI best_OSI(1)];
rho = rho/max(rho);
polarplot(theta,rho,'LineWidth',1)
ax = gca;
ax.RTick = [0.5 1];
ax.RLim = [0 1.25];
ax.ThetaDir = 'clockwise';
title('orientation selectivity')

figure_size = [300 250];
figure('Position',[500 100 figure_size]);
theta = deg2rad(0:30:360);
rho = [best_DSI best_DSI(1)];
rho = rho/max(rho);
polarplot(theta,rho,'LineWidth',1)
ax = gca;
ax.RTick = [0.5 1];
ax.RLim = [0 1.25];
ax.ThetaDir = 'clockwise';
title('direction selectivity')


%% plot before/after cumulative distribution plots
colors = linspecer(2);
%spontaneous activity (spikes per min)
plot_size = [1000 250];
figure('Position',[100 100 plot_size]);
subplot(1,3,1);
max_tpm = max(roi_spikespermin,[],'all');
x = 0:max_tpm;
for d = 1:2
    y = zeros(size(x));
    roi_inds = find(~isnan(roi_info(:,d)));
    yvals = floor(sort(roi_spikespermin(roi_inds,d)/1000));
    for i = yvals
        y(i+1) = y(i+1)+1;
    end
    y = cumsum(y);
    y = y/max(y);
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
% legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('spikes/min')
ylabel('cumulative probability')
xlim([0 200])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
title('spontaneous activity')

%spontaneous activity transients/min
% plot_size = [800 250];
% figure('Position',[100 100 plot_size]);
% subplot(1,3,1);
% title('spontaneous activity')
% max_tpm = max(roi_transientspermin,[],'all');
% x = 0:0.1:max_tpm;
% for d = 1:2
%     y = zeros(size(x));
%     roi_inds = find(~isnan(roi_info(:,d)));
%     yvals = floor(10*sort(roi_transientspermin(roi_inds,d)));
%     for i = yvals
%         y(i+1) = y(i+1)+1;
%     end
%     y = y/10;
%     y = cumsum(y);
%     y = y/max(y);
%     plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
%     hold on
% end
% % legend({'baseline','treated'},'Location','NorthEastOutside')
% xlabel('transients/min')
% ylabel('cumulative probability')
% xlim([20 40])
% ylim([0 1])
% yticks([0 0.2 0.4 0.6 0.8 1])

%orientation selectivity
% figure('Position',[100 100 plot_size]);
subplot(1,3,2)
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
subplot(1,3,3)
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


%% plot data
plot_size = [400 250];
colors = linspecer(2);

%plot before/after spontaneous and stimulated activity scatterplot
figure('Position',[100 100 plot_size]);
title('neural excitability')
for d = 1:4
    scatter(1:num_rois',roi_transientspermin(:,d),20,'filled')
    hold on
end
legend({'before-spont','after-spont','before-stim','after-stim'},'Location','NorthEastOutside')
xlabel('roi')
ylabel('transients/min')

%plot before/after spontaneous activity sorted lineplot
figure('Position',[100 100 plot_size]);
title('neural excitability')
for d = 1:2
    roi_inds = find(~isnan(roi_info(:,d)));
    nr = length(roi_inds);
    x = (0:(nr-1))/(nr-1);
    y = sort(roi_transientspermin(roi_inds,d));
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('roi (sorted)')
ylabel('transients/min')
xlim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])

%plot before/after OSI sorted lineplot
figure('Position',[100 100 plot_size]);
title('orientation selectivity')
for d = 1:2
    roi_inds = find(~isnan(roi_info(:,d)));
    nr = length(roi_inds);
    x = (0:(nr-1))/(nr-1);
    y = sort(roi_OSI(roi_inds,d));
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('roi (sorted)')
ylabel('OSI')
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.2 0.4 0.6 0.8 1])

%plot before/after DSI sorted lineplot
figure('Position',[100 100 plot_size]);
title('direction selectivity')
for d = 1:2
    roi_inds = find(~isnan(roi_info(:,d)));
    nr = length(roi_inds);
    x = (0:(nr-1))/(nr-1);
    y = sort(roi_DSI(roi_inds,d));
    plot(x,y,'-k','LineWidth',1.5,'Color',colors(d,:));
    hold on
end
legend({'baseline','treated'},'Location','NorthEastOutside')
xlabel('roi (sorted)')
ylabel('DSI')
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.2 0.4 0.6 0.8 1])

