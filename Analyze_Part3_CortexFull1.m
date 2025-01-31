%% settings
analysis_folder = "D:\SN Lab\Visual Cortex\Raw Data\2022-10-10 16-13 22_06_30_07_01 full20x_F_APP_ly6g_mouse_L_3331";


%% load files and metadata
% load experiment metadata
filename = fullfile(analysis_folder,'exp_metadata.mat');
assert(isfile(filename),'exp_metadata.mat file note found in analysis folder')
load(filename);

% load roi_info
filename = fullfile(analysis_folder,'roi_info.mat');
assert(isfile(filename),'roi_info.mat file note found in analysis folder')
load(filename);
if exist('roi_info','var')
    single_roi_info = true;
    num_timepoints = size(roi_info,3);
    num_rois = nan(1,num_timepoints);
    for t = 1:num_timepoints
        num_rois(t) = sum(~isnan(roi_info(:,1,t)));
    end
elseif exist('roi_info_1','var')
    single_roi_info = false;
    t = 0;
    endreached = 0;
    while endreached==0
        t = t+1;
        roi_info_name = ['roi_info_' num2str(t)];
        if exist(roi_info_name,'var')
            num_timepoints = t;
            num_rois(t) = size(eval(roi_info_name),1);
        else
            endreached = 1;
        end
    end        
else
    error('unexpected variables in roi_info.mat')
end

%load suite2p settings to get fps
filename = fullfile(analysis_folder,num2str(1),'suite2p','suite2p','plane0','Fall.mat');
assert(isfile(filename),['suite2p Fall.mat file note found in analysis folder for timepoint ' num2str(1)])
load(filename);
fps = ops.fs;

% experiment settings
numRepetitions = vstruct(1).num_reps;
repetitionsToAverage = 2:numRepetitions; %skip first, it tends to have an outsized effect
stimOnFrames = 8:17;
baselineFrames = 1:6;
spontFrames = 100:mdata(1).num_spont_frames;
numConditions = vstruct(1).num_trials;
numFramesPerTrial = mdata(1).num_stim_frames;
stimStart = 1+mdata(1).num_spont_frames;
spontMinutes = mdata(1).num_spont_frames/(fps*60);
stimMinutes = (mdata(1).num_stim_frames*mdata(1).num_stims)/(fps*60);

column_names = {'1: cell type (user classified)';...
                '2: center x';...
                '3: center y';...
                '4: spikes-per-minute';...
                '5: OSI';...
                '6: DSI';...
                '7: gcamp brightness';...
                '8: eyfp brightness';...
                '9: mruby2 brightness';...
                '10: transients-per-minute-f';...
                '11: transients-per-minute-s';...
                '12: roi # (of tracked cells)';...
                '13: average Fisher information'};
            
            
%% load and process cell data for all timepoints
cell_data = nan(max(num_rois),12,num_timepoints); %[num_rois, num_datatypes, num_timepoints]
for t = 1:num_timepoints
    % load suite2p (for activity, OSI, DSI, x, y, [masks]
    filename = fullfile(analysis_folder,num2str(t),'suite2p','suite2p','plane0','Fall.mat');
    assert(isfile(filename),['suite2p Fall.mat file note found in analysis folder for timepoint ' num2str(t)])
    load(filename);
    FmFneu = F-(0.7*Fneu);
    
    % load unmixed colors (use masks for brightness)
    filename = fullfile(analysis_folder,num2str(t),'celltypes_unmixed.tif');
    assert(isfile(filename),['celltypes_unmixed.tif file note found in analysis folder for timepoint ' num2str(t)])
    unmixed = read_file(filename);
    
    % analyze spontaneous activity 
    num_frames = length(spontFrames);

    %konnerth-like measurement of "transients")
    spontF = FmFneu(:,spontFrames);
    estBaselineF = prctile(spontF,50,2); %%%%%%%%%%%improve this?
    shiftF = spontF-repmat(estBaselineF,[1 num_frames]);    
    noiseF = shiftF;
    noiseF(noiseF>0) = nan;
    noiseF = [noiseF -noiseF];
    noiseStd = std(noiseF,0,2,'omitnan');
    estThresholdF = estBaselineF + (3*noiseStd);
    aboveThresholdF = spontF>repmat(estThresholdF,[1 num_frames]);
    transientsF = aboveThresholdF & logical([ones(size(F,1),1) diff(aboveThresholdF,1,2)==1]); %turn spans of 1's into a single 1
    transientsperminute_f = sum(transientsF,2)/spontMinutes;

    %spike-based "transients"
    transientsperminute_s = sum(spks(:,spontFrames)>0,2)/spontMinutes;

    %spikes-per-minute
    spikesperminute = sum(spks(:,spontFrames),2)/spontMinutes;
    cell_data(1:num_rois(t),4,t) = spikesperminute;
    cell_data(1:num_rois(t),10,t) = transientsperminute_f;
    cell_data(1:num_rois(t),11,t) = transientsperminute_s;
    
    % analyze stimulus-evoked activity data
    stimulus_data = nan(num_rois(t),numFramesPerTrial,numRepetitions,numConditions);
    for r = 1:numRepetitions
        for c = 1:numConditions
            %get info about the current experimental condition
            cond = vstruct(t).order(r,c);
            start_frame = stimStart + ((c-1)*numFramesPerTrial) + ((r-1)*numFramesPerTrial*numConditions);
            frame_inds = start_frame:(start_frame+numFramesPerTrial-1);
            
            %get fluorescence traces for all ROIs for this condition
            stimulus_data(:,:,r,cond) = FmFneu(:,frame_inds);
        end
    end
    stimDF = 100*((mean(stimulus_data(:,stimOnFrames,:,:),2)./mean(stimulus_data(:,baselineFrames,:,:),2))-1); %mean dF/F from 0-2s over baseline (-2-0s)
    stimDF = permute(stimDF,[1 3 4 2]); %[roi, rep, cond]
    meanStimDF = mean(stimDF(:,repetitionsToAverage,:),2,'omitnan');
    meanStimDF = permute(meanStimDF,[1 3 2]); %[roi, cond]
    stdStimDF = std(stimDF(:,repetitionsToAverage,:),0,2,'omitnan');
    stdStimDF = permute(stdStimDF,[1 3 2]); %[roi, cond]
    angles = deg2rad(0:30:330);
    
    for roi = 1:num_rois(t)
        meanActivity = meanStimDF(roi,1:12);
        
        %OSI
        OSI_complex = sum(meanActivity.*(exp(2*1i*angles)))./sum(abs(meanActivity));
        OSI = abs(OSI_complex);
        cell_data(roi,5,t) = OSI;

        %DSI
        DSI_complex = sum(meanActivity.*(exp(1i*angles)))./sum(abs(meanActivity));
        DSI = abs(DSI_complex);
        cell_data(roi,6,t) = DSI;
        
        %get cell x,y positions and brightnesses
        xpix = stat{roi}.xpix+1;
        ypix = stat{roi}.ypix+1;
        
        cell_data(roi,2,t) = mean(xpix); %col 2: center x
        cell_data(roi,3,t) = mean(ypix); %col 3: center y
        
        cell_data(roi,7,t) = mean(unmixed(sub2ind(size(unmixed),ypix,xpix,1*ones(size(xpix)))),'omitnan'); % col 7: mean gcamp brightness
        cell_data(roi,8,t) = mean(unmixed(sub2ind(size(unmixed),ypix,xpix,2*ones(size(xpix)))),'omitnan'); % col 8: mean eyfp brightness
        cell_data(roi,9,t) = mean(unmixed(sub2ind(size(unmixed),ypix,xpix,3*ones(size(xpix)))),'omitnan'); % col 9: mean mruby2 brightness
        
        %get fisher information score
        dMu_ds = diff(meanStimDF(roi,1:12));  % Numerical derivative of mean response w.r.t. stimulus
        dMu_ds(12) = meanStimDF(roi,1)-meanStimDF(roi,12); %circle back to the start
        FI = (dMu_ds.^2) ./ stdStimDF(roi,1:12);  % FI per stimulus
        FI = mean(FI); % average across stimulus to get average Fisher Information
        cell_data(roi,12,t) = FI;

        %get user classification of cell type
        if single_roi_info
            cell_data(roi,1,t) = roi_info(roi,2,t);
        else
            roi_info_name = ['roi_info_' num2str(t)];
            roi_info_t = eval(roi_info_name);
            cell_data(roi,1,t) = roi_info_t(roi,2);
        end
    end
end


%% organize data for cells tracked across all timepoints
% num_tracked_cells = 0;
% tracked_cell_data = nan(num_tracked_cells,10,num_timepoints);

for t = 1:num_timepoints
    if single_roi_info
        roi_info_t = roi_info(:,[1 3],t);
    else
        roi_info_name = ['roi_info_' num2str(t)];
        roi_info_t = eval(roi_info_name);
        roi_info_t = roi_info_t(:,[1 3]);
    end
    
    t_tracked_rois = find(all(~isnan(roi_info_t),2));
    roi_info_t(t_tracked_rois,:);
    size(t_tracked_rois);
    
    if t==1
        rt1 = roi_info_t(t_tracked_rois,:);
    else
        rt2 = roi_info_t(t_tracked_rois,:);
        %%%%%%%%%%% for now, only use rt2 until roi_info problem is fixed
        
        %check for duplicates in t1 rois
        if length(unique(rt2(:,2)))<length(unique(rt2(:,1)))
            [~, unique_inds, ~] = unique(rt2(:,2)); %indices of unique numbers
            %%%%%%%%%%%% for now, only use the 1st one
            unique_inds = sort(unique_inds);
            rt2 = rt2(unique_inds,:);
        end
        tracked_rois = fliplr(rt2);
        num_tracked_rois = size(rt2,1);
    end
end

% sort tracked rois by increasing combined roi #
[~, order] = sort(sum(tracked_rois,2));
tracked_rois = tracked_rois(order,:);

%loop through number of pairs
tracked_cell_data = nan(num_tracked_rois,13,num_timepoints); %[num_rois, num_datatypes, num_timepoints]
for r = 1:num_tracked_rois
    %for each pair, populate tracked data with cell data
    for t = 1:num_timepoints
        roi = tracked_rois(r,t);
        tracked_cell_data(r,[1:11 13],t) = cell_data(roi,:,t);
        tracked_cell_data(r,12,t) = roi;
    end
end


%% save data
filename = fullfile(analysis_folder,'cell_data.mat');
save(filename,'cell_data','tracked_cell_data','column_names');


