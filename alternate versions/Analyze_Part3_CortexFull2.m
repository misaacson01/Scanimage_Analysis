%currently calculates some simple summary data, e.g.
%number of total active neurons
%transients-per-minute of every cell
%OSI, DSI of every cell
%plots highly active spontaneous cells
%plots highly selective cells

%% settings
% analysis_folder = 'C:\Users\misaa\Desktop\2022-10-17 14-50 22_08_30_31 full20x_M_APP_ly6g_mouse_R_spot_2';
analysis_folder = "C:\Users\misaa\Desktop\2022-10-10 16-13 22_06_30_07_01 full20x_F_APP_ly6g_mouse_L";
savename = 'BA1016.mat';
fps = 3.4;

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
num_frames = size(F,2);
        
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

%for loading cell linking data
%
%
%

%summary of roi data
roi_summary_1 = nan(size(iscell_1,1),6); %[iscell1 tpm OSI DSI color roi2]
roi_summary_2 = nan(size(iscell_2,1),6); %[iscell2 tpm OSI DSI color roi1]
roi_summary_1(:,1) = iscell_1(:,1); 
roi_summary_2(:,1) = iscell_2(:,1); 

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
repetitionsToAverage = 2:numRepetitions; %skip first, it tends to have an outsized effect
stimOnFrames = 8:17;
baselineFrames = 1:6;
numConditions = vstruct(1).num_trials;
numFramesPerTrial = mdata(1).num_stim_frames;
stimStart = 1+mdata(1).num_spont_frames;
spontMinutes = mdata(1).num_spont_frames/(fps*60);
stimMinutes = (mdata(1).num_stim_frames*mdata(1).num_stims)/(fps*60);
% roi_stimspks = nan([numRepetitions, numConditions, mdata(1).num_stim_frames]);
% roi_stimtransients = nan([numRepetitions, numConditions, mdata(1).num_stim_frames]);
% roi_stimkontransients = nan([numRepetitions, numConditions, mdata(1).num_stim_frames]);
roi_DSI = nan(num_rois,2);
roi_OSI = nan(num_rois,2);


%% un-randomize stimulus data
%create frame indices for each trial, organized by [repetition, conditions, frame]
trial_order = repmat(0:(numConditions-1),[numRepetitions 1]);
trial_order = trial_order + repmat(numConditions*(0:numRepetitions-1)',[1 numConditions]); %value for trial1 rep2 is now numConds+1, eg
start_of_trial_frame = trial_order*numFramesPerTrial; %values are now frame numbers, relative to the first trial of each repetition
%start_of_trial_frame = the first frame # for each repetition and
%condition, if the data were NOT randomized

roi_stimF_1 = nan(num_rois,numRepetitions,numConditions,numFramesPerTrial);
roi_stimF_2 = nan(num_rois,numRepetitions,numConditions,numFramesPerTrial);
for t = 1:2 %loop for both
    %sort by increasing condition
    [~, sorted_col_order] = sort(vstruct(t).order,2); %get the randomized order of this dataset
    sorted_row_order = repmat((1:numRepetitions)',[1 numConditions]); %repetition # for each trial
    sorted_trial_start_frame = start_of_trial_frame(sub2ind([numRepetitions numConditions],sorted_row_order,sorted_col_order));
    sorted_trial_start_frame = sorted_trial_start_frame + stimStart; %values are the actual frame numbers when each trial started
    %sorted_trial_start_frame = the first frame # for each repetition and
    %condition, of data that WAS randomized
    
    stimFrames = repmat(sorted_trial_start_frame, [1 1 numFramesPerTrial]); %3rd dimension will store all frames for each trial
    stimFrames = stimFrames + repmat(permute((0:numFramesPerTrial-1),[1 3 2]),[numRepetitions, numConditions, 1]);
    %stimFrames = similar to sorted_trial_start_frame, except that it
    %extends in the 3rd dimension to all frames of that trial

    %get suite2p data for this dataset
    if t==1
        spks = spks_1;
        F = F_1;
        Fneu = Fneu_1;
        roi_summary = roi_summary_1;
    else
        spks = spks_2;
        F = F_2;
        Fneu = Fneu_2;
        roi_summary = roi_summary_2;
    end
    FmFneu = F-(0.7*Fneu);
    
    %for konnerth-like measurement of "transients":
    spontF = F(:,100:mdata(t).num_spont_frames);
    estBaselineF = prctile(spontF,50,2);
    shiftF = spontF-estBaselineF;
    shiftF(shiftF>0) = nan;
    shiftF = [shiftF -shiftF];
    estNoiseF = std(shiftF,0,2,'omitnan');
    estThresholdF = estBaselineF + (3*estNoiseF);
    aboveThresholdF = F>repmat(estThresholdF,[1 num_frames]);
    transientsF = aboveThresholdF & logical([ones(size(F,1),1) diff(aboveThresholdF,1,2)==1]); %turn spans of 1's into a single 1
            
    %loop for every roi
    rois = roi_info(:,t); %all rois classified as cells in suite2p
    rois2 = find(roi_summary(:,1)==1);
    for r = 1:num_rois
        roi = rois(r); %get current roi #
        if ~isnan(roi)
            %spontaneous spikes/min
%             roi_spikespermin(r,t) = sum(spks(roi,1:mdata(t).num_spont_frames))/spontMinutes; %spontaneous 
%             roi_spikespermin(r,2+t) = sum(spks(roi,1+mdata(t).num_spont_frames:end))/stimMinutes; %stim-evoked
            
            %spontaneous spike transients/min
%             roi_transientspermin(r,t) = sum(spks(roi,1:mdata(t).num_spont_frames)>0)/spontMinutes; %spontaneous 
%             roi_transientspermin(r,2+t) = sum(spks(roi,1+mdata(t).num_spont_frames:end)>0)/stimMinutes; %stim-evoked 
            
            %spontaneous konnerth-like transients/min
            tpm = sum(transientsF(roi,1:mdata(t).num_spont_frames))/spontMinutes;
            roi_kontransientspermin(r,t) = tpm;
            roi_kontransientspermin(r,2+t) = sum(transientsF(roi,1+mdata(t).num_spont_frames:end))/stimMinutes;
            roi_summary(roi,2) = tpm;
            
            %organize spikes, transients, F, Fneu by [repetition, trial, frame, roi]
%             roi_stimspks = spks(sub2ind(size(F),roi*ones(size(stimFrames)),stimFrames));
%             roi_stimtransients = spks(sub2ind(size(F),roi*ones(size(stimFrames)),stimFrames))>0;
%             roi_stimkontransients = transientsF(sub2ind(size(F),roi*ones(size(stimFrames)),stimFrames));
            roi_stimF = FmFneu(sub2ind(size(F),roi*ones(size(stimFrames)),stimFrames));
            if t==1
                roi_stimF_1(r,:,:,:) = roi_stimF;
            else
                roi_stimF_2(r,:,:,:) = roi_stimF;
            end
            
            %average Fluorescence data for repetitions
            StimActivity = 100*((mean(roi_stimF(:,:,stimOnFrames),3)./mean(roi_stimF(:,:,baselineFrames),3))-1); % %mean dF/F from 0-2s over baseline (-2-0s)
            meanStimActivity = mean(StimActivity(repetitionsToAverage,:),1); %average repetitions
            
            %OSI
            tmpSpksO = mean([meanStimActivity(1:6); meanStimActivity(7:12)]);
            [maxActivity,maxI] = max(tmpSpksO);
            oppI = 1+mod((maxI+3)-1,6);
            oppActivity = tmpSpksO(oppI);
            OSI = (maxActivity-oppActivity)/(maxActivity+oppActivity);
            OSI = min([max([0 OSI]) 1]);
            roi_OSI(r,t) = OSI;
            roi_summary(roi,3) = OSI;
            
            %DSI
            [maxActivity,maxI] = max(meanStimActivity(1:12)); %condition which cause the most spikes
            oppI = 1+mod((maxI+6)-1,12); %condition on the opposite angle
            oppActivity = meanStimActivity(oppI);
            DSI = (maxActivity-oppActivity)/(maxActivity+oppActivity);
            DSI = min([max([0 DSI]) 1]);
            roi_DSI(r,t) = DSI;
            roi_summary(roi,4) = DSI;
        end
    end
    if t==1
        roi_summary_1 = roi_summary;
    else
        roi_summary_2 = roi_summary;
    end
end

% %% find best cells, plot the dF/F for all repetitions
% figure('Position',[100 100 500 400])
% 
% %for best active rois, plot spontaneous dF/F
% timeX = (0.5:mdata(1).num_spont_frames)/fps;
% spontF = F_2(:,1:mdata(1).num_spont_frames) - 0.7*Fneu_2(:,1:mdata(1).num_spont_frames);
% spontF = 100*((spontF./repmat(prctile(spontF(:,100:end),50,2),[1 mdata(1).num_spont_frames]))-1);
% spontF = movmean(spontF,3,2);
% %pick most active rois
% t=2;
% [sorted_roi_active, sorted_order] = sort(roi_kontransientspermin(:,t),'descend');
% % handpickedrois = [3 4 5 6 7   13 14 17 18 19    20 21 27 28 29    110]; %best active
% % handpickedrois = [149 153 142 143 144]; %good quiet
% handpickedrois = [3 4 149 5 6    153 110 13 142 14    17 18 19 143 20    21 144 27 28 29 7];
% best_active = sorted_order(handpickedrois)';
% % best_active = best_active(randperm(16));
% spw = 14;
% for roi = 1:16
%     r = best_active(roi);
%     subplot(16,spw,1 + spw*(roi-1));
%     plot(timeX,zeros(size(timeX)),'--k','Color',[0.7 0.7 0.7],'LineWidth',0.5)
%     hold on
%     plot(timeX,spontF(r,:),'b','LineWidth',1)
%     xlim([0 18])
%     ylim([0 150])
%     ax = gca;
%     ax.Clipping = 'off';
%     axis off
%     text(-20,0, num2str(roi))
%     subplot(16,spw,spw + spw*(roi-1));
%     text(0.5,0.2,sprintf('%3.1f',(roi_kontransientspermin(r,2))),'FontSize',8);
%     axis off
%     if roi==1
%         text(-0.8,1.6,'transients/min','FontSize',8,'FontWeight','bold');
%     end
% end
% %%
% figure('Position',[100 100 500 400])
% % t = 2;
% % [sorted_roi_DSI, sorted_order] = sort(roi_DSI(:,t),'descend');
% % handpickedrois = [2 3 6 9 10   13 15 17 18 20   21 22 23 26 27]; %best DSI
% % best_DSI = sorted_order(handpickedrois)';
% % 
% % [sorted_roi_OSI, sorted_order] = sort(roi_OSI(:,t),'descend');
% % handpickedrois = [2 4 7 10 14 16]; %best OSI
% % best_OSI = sorted_order(handpickedrois)';
% % 
% % %best ROIs:
% % best_both = unique(sort([best_DSI best_OSI]));
% % best_both = best_both(randperm(length(best_both)));
% 
%  %manually selected for 2022-10-17 14-50 22_08_30_31 full20x_M_APP_ly6g_mouse_R_spot_2
%  
% subplot(1,2,2)
% best_both = [77 16 34 53 4     32 1 10 69 48   115 3 26 14 132   70 ];
% timeX = ((1:numFramesPerTrial)/fps)-2;
% spw = 14;
% figure()
% for roi = 1:16
%     r = best_both(roi);
%     for cond = 1:12
%         meanF = squeeze(mean(roi_stimF_2(r,repetitionsToAverage,cond,:),2));
%         baselineF = mean(roi_stimF_2(r,repetitionsToAverage,:,baselineFrames),'all');
%         meanDF = 100*((meanF/baselineF)-1);
%         subplot(16,spw,cond + spw*(roi-1));
%         plot([-2 4],[0 0],'--k','Color',[0.7 0.7 0.7],'LineWidth',0.5)
%         hold on
%         patch([0 2 2 0],[0 0 150 150],'k','EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.1)
%         plot(timeX,meanDF,'b','LineWidth',1)
%         xlim([-2 4])
%         ylim([0 150])
%         axis off
%         ax = gca;
%         ax.Clipping = 'off';
%         if roi==16 && cond==1
%             line([-4 -4],[0 100])
%         end
%         if roi==1            
%             ax = gca;
%             xlim([-2 4])
%             ylim([0 150])
% 
%             angles = 0:30:330;
%             angle = deg2rad(angles(cond));
%             arrowX = [1-2.5*cos(angle) 1+2.5*cos(angle)];
%             arrowY = [250+75*sin(angle) 250-75*sin(angle)];
% 
%             axisXRange = diff(ax.XLim);
%             axisYRange = diff(ax.YLim);
%             axisXStart = ax.XLim(1);
%             axisYStart = ax.YLim(1);
%             normXRange = ax.Position(3);
%             normYRange = ax.Position(4);
%             normXStart = ax.Position(1);
%             normYStart = ax.Position(2);
%             arrowX = normXRange*((arrowX - axisXStart)/axisXRange) + normXStart;
%             arrowY = normYRange*((arrowY - axisYStart)/axisYRange) + normYStart;
%             
%             ar = annotation('arrow',arrowX,arrowY);
%             ar.Color = 'blue';
%             ar.LineWidth = 1;
%             ar.HeadLength = 5;
%             ar.HeadWidth = 7;
%             ar.LineWidth = 1.5;
%         end
%     end
%     %OSI
%     subplot(16,spw,13 + spw*(roi-1));
%     text(0.5,0.2,sprintf('%0.2f',(roi_OSI(r,2))),'FontSize',8);
%     axis off
%     if roi==1
%         text(0.5,1.6,'OSI','FontSize',8,'FontWeight','bold');
%     end
%     %DSI
%     subplot(16,spw,14 + spw*(roi-1));
%     text(0.5,0.2,sprintf('%0.2f',(roi_DSI(r,2))),'FontSize',8);
%     axis off
%     if roi==1
%         text(0.5,1.6,'DSI','FontSize',8,'FontWeight','bold');
%     end
% end
% 

%% plot before/after silent/normal/hyperactive
bins = zeros(2,3); %[time, bin]
for t = 1:2
    spm = roi_kontransientspermin(:,t);
    spm(isnan(spm)) = [];
    bins(t,1) = bins(t,1) + sum(spm<=0.25);
    bins(t,2) = bins(t,2) + sum(spm>0.25 & spm<=4);
    bins(t,3) = bins(t,3) + sum(spm>4);
end
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


%% plot before/after plots
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
colors = linspecer(2);
subplot(1,4,2);
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
     
           
%% load roi info (for 2)
load(fullfile(analysis_folder,'roi_info.mat'))
ri = roi_info_2;

linked_cells_2 = ~isnan(ri(:,1)) & ~isnan(ri(:,3));
num_linked = sum(linked_cells_2);
linked_cellnums = ri(find(linked_cells_2),:);
linked_cellnums = fliplr(linked_cellnums); %now before is the 1st column, after is 3rd, color is middle

tpm_ba_red = [];
DSI_ba_red = [];
OSI_ba_red = [];

tpm_ba_yel = [];
DSI_ba_yel = [];
OSI_ba_yel = [];

tpm_ba_bth = [];
DSI_ba_bth = [];
OSI_ba_bth = [];

tpm_ba_all = [];
DSI_ba_all = [];
OSI_ba_all = [];

for i = 1:length(linked_cellnums)
    before_cellnum = linked_cellnums(i,1);
    before_iscellnum = find(roi_info(:,1)==before_cellnum);
    after_cellnum = linked_cellnums(i,3);
    after_iscellnum = find(roi_info(:,2)==after_cellnum);
    color = linked_cellnums(i,2);
    
    DSI_ba = [roi_DSI(before_iscellnum,1) roi_DSI(after_iscellnum,2)];
    OSI_ba = [roi_OSI(before_iscellnum,1) roi_OSI(after_iscellnum,2)];
    tpm_ba = [roi_kontransientspermin(before_iscellnum,1) roi_kontransientspermin(after_iscellnum,2)];
    
    tpm_ba_all = [tpm_ba_all;tpm_ba];
    DSI_ba_all = [DSI_ba_all;DSI_ba];
    OSI_ba_all = [OSI_ba_all;OSI_ba];
    
    if color==1
        tpm_ba_yel = [tpm_ba_yel;tpm_ba];
        DSI_ba_yel = [DSI_ba_yel;DSI_ba];
        OSI_ba_yel = [OSI_ba_yel;OSI_ba];
    elseif color==2
        tpm_ba_red = [tpm_ba_red;tpm_ba];
        DSI_ba_red = [DSI_ba_red;DSI_ba];
        OSI_ba_red = [OSI_ba_red;OSI_ba];
    elseif color==3
        tpm_ba_bth = [tpm_ba_bth;tpm_ba];
        DSI_ba_bth = [DSI_ba_bth;DSI_ba];
        OSI_ba_bth = [OSI_ba_bth;OSI_ba];
    end
end
ba_all = tpm_ba_all;
ba_all(:,:,2) = OSI_ba_all;
ba_all(:,:,3) = DSI_ba_all;

ba_yel = tpm_ba_yel;
ba_yel(:,:,2) = OSI_ba_yel;
ba_yel(:,:,3) = DSI_ba_yel;

ba_red = tpm_ba_red;
ba_red(:,:,2) = OSI_ba_red;
ba_red(:,:,3) = DSI_ba_red;

save(['C:\Users\misaa\Desktop\' savename ],'ba_all','ba_yel','ba_red')