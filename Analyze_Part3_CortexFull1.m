
%currently outputs some simple summary data, e.g.
%number of yellow, red, total active neurons
%number of neurons linked between before/after datasets
%transients-per-minute of every roi in a big matrix: 
%   [before-spontaneous, after-spontaneous, before-stim, after-stim]
%same data as the matrix, but in a scatterplot

%% settings
analysis_folder = 'C:\Users\misaa\Desktop\2022-06-15 06-12 22_5_8_9_ly6g_test';
fps = 3.6;

%enter roi info for all rois here: [roi#before, roi#after, celltype]
%(celltype = 0 for unknown, 1 for yellow, 2 for red)
roi_info = [16 nan 0;...
    10 28 1;...
    19 16 0;...
    18 54 0;...
    9 nan 0;...
    21 nan 0;...
    11 nan 0;...
    6 nan 0;...
    1 3 0;...
    2 10 0;...
    4 nan 0;...
    15 5 0;...
    8 2 0;...
    5 nan 0];
    
%check for duplicates
for t = 1:2
    tmp = roi_info(:,t);
    tmp(isnan(tmp)) = [];
    assert(length(tmp)==length(unique(tmp)),['redundant roi(s) found in roi_info (column ' num2str(t) ')']);
end


%% load files 
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
cells50 = find(iscell_2(:,2));
num_rois = length(cells50);
colors = linspecer(num_rois);
o = randperm(num_rois);
colors_2 = colors(o,:);

load(fullfile(analysis_folder,'exp_metadata.mat'));


%% analyze data
%summary data of roi info
num_active = sum(~isnan(roi_info(:,1:2)))
num_active_yellow = sum(~isnan(roi_info(:,1:2)) & repmat(roi_info(:,3)==1,[1 2]))
num_active_red = sum(~isnan(roi_info(:,1:2)) & repmat(roi_info(:,3)==2,[1 2]))
num_linked = sum(sum(~isnan(roi_info(:,1:2)),2)==2)
num_linked_yellow = sum(sum(~isnan(roi_info(:,1:2)),2)==2 & roi_info(:,3)==1)
num_linked_red = sum(sum(~isnan(roi_info(:,1:2)),2)==2 & roi_info(:,3)==2)

num_rois = size(roi_info,1);
roi_spikespermin = nan(num_rois,4);
roi_transientspermin = nan(num_rois,4);
for t = 1:2
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
            roi_spikespermin(r,t) = 60*sum(spks(roi,1:mdata(t).num_spont_frames))/fps; %spontaneous 
            roi_spikespermin(r,2+t) = 60*sum(spks(roi,1+mdata(t).num_spont_frames:end))/fps; %stim-evoked
            
            %transients/min
            roi_transientspermin(r,t) = 60*sum(spks(roi,1:mdata(t).num_spont_frames)>0)/fps; %spontaneous 
            roi_transientspermin(r,2+t) = 60*sum(spks(roi,1+mdata(t).num_spont_frames:end)>0)/fps; %spontaneous 
            
            %OSI

            %DSI
        end
    end
end

roi_transientspermin

title('spontaneous activity')
scatter(repmat((1:num_rois)',[1 4]),roi_transientspermin,20,'filled')
legend({'before-spont','after-spont','before-stim','after-stim'},'Location','NorthEastOutside')
xlabel('roi')
ylabel('transients/min')