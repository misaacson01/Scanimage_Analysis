
%set analysis folder
analysis_folder = 'C:\Users\misaa\Desktop\Awake Cortex Tests\2025-05-01 14-14 2024-09-24 APP Control';


%% load files and metadata
% load experiment metadata
filename = fullfile(analysis_folder,'exp_metadata.mat'); 
assert(isfile(filename),'exp_metadata.mat file note found in analysis folder')
load(filename); %load mdata and vstruct

%load mouse metadata
filename = fullfile(analysis_folder,'mouse_metadata.csv'); 
assert(isfile(filename),'mouse_metadata.csv file note found in analysis folder')
mouse_metadata = readtable(fullfile(analysis_folder,'mouse_metadata.csv'));

%load plaque data
filename = fullfile(analysis_folder,'plaque_results.mat');
assert(isfile(filename),'plaque_results.mat file note found in analysis folder')
load(filename); %load results

%load suite2p data
filename = fullfile(analysis_folder,num2str(1),'suite2p','suite2p','plane0','Fall.mat');
assert(isfile(filename),'suite2p Fall.mat file note found in analysis folder of timepoint 1')
load(filename); %load F, Fneu, iscell, ops, redcell, spks, stat

%experiment settings
spontaneousSecondsToIgnore = 5;
fps = ops.fs;
numFramesToIgnore = spontaneousSecondsToIgnore*fps;
numRepetitions = vstruct(1).num_reps;
repetitionsToAverage = 2:numRepetitions; %skip first, it tends to have an outsized effect
stimOnFrames = 8:17;
baselineFrames = 1:6;
numSpontFrames = mdata(1).num_spont_frames;
numConditions = vstruct(1).num_trials;
numFramesPerTrial = mdata(1).num_stim_frames;
numStimFrames = mdata(1).num_stim_frames*mdata(1).num_stims;
stimulus_angles = deg2rad(0:30:330);
preSpontFrames = round(numFramesToIgnore+1):numSpontFrames; %get spontaneous frames, but ignore the first 5 seconds
postSpontFrames = round(numSpontFrames+numStimFrames+numFramesToIgnore+1):((2*numSpontFrames)+numStimFrames);
preStimFrames = (1+numSpontFrames):(numSpontFrames+numStimFrames);
postStimFrames = ((2*numSpontFrames)+numStimFrames+1):(2*(numSpontFrames+numStimFrames));
spontMinutes = length(preSpontFrames)/(fps*60);
stimMinutes = length(preStimFrames)/(fps*60);
cellRows = find(iscell(:,1));
numCells = length(cellRows);

%organization for long data format
column_names = {'roi #';...                         %1 (e.g. 1)
                'mouse name';...                    %2 (e.g. 1R2L_4294)
                'sex';...                           %3 (e.g. f)
                'order';...                         %4 (e.g. control-ly6g)
                'treatment';...                     %5 (e.g. control)
                'timepoint';...                     %6 (e.g. pre)
                'genotype';...                      %7 (e.g. app)
                'cell type';...                     %8 (e.g. excitatory)
                'center x';...                      %9
                'center y';...                      %10
                'distance to nearest plaque';...    %11
                'distance to nearest stall';...     %12
                'gcamp brightness';...              %13
                'spont spikes per minute';...       %14
                'spont transients per minute';...   %15
                'stim spikes per minute';...        %16
                'stim transients per minute';...    %17
                'OSI';...                           %18
                'preferred orientation';...         %19
                'orientation deviation';...         %20
                'DSI';...                           %21
                'preferred direction';...           %22
                'direction deviation';...           %23
                'average Fisher information';...    %24   
                'FI within stim variance';...       %25
                'FI between stim variance'};        %26   
variable_types = {'double';...  % 1 roi #
    'string';...                % 2 mouse name
    'string';...
    'string';...
    'string';...
    'string';...
    'string';...
    'string';...                % 8 cell type
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...                % 16 
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...
    'double';...                % 24
    'double';...                % 25
    'double';...                % 26
    };
numDatatypes = length(column_names);


%% process cell data for both timepoints
data_table = table('Size',[numCells*2, numDatatypes],'VariableTypes',variable_types,'VariableNames',column_names);
corr_matrix = nan([numCells numCells 4]); %[cell, cell, {prespont prestim postspont poststim}]
FmFneu = F - 0.7*Fneu;

%%% process spike data for all other measurements
for t = 1:2
    if t==1
        spontFrames = preSpontFrames;
        stimFrames = preStimFrames;
        dataRows = 1:numCells;
        timepoint = 'pre';
    else
        spontFrames = postSpontFrames;
        stimFrames = postStimFrames;
        dataRows = (1+numCells):(2*numCells);
        timepoint = 'post';
    end

    %%% calculate correlation matrices
    spontFmFneu = FmFneu(cellRows,spontFrames);
    stimFmFneu = FmFneu(cellRows,stimFrames);
    for i = 1:numCells %loop through all cells
        for j = (i+1):numCells  %loop through upper triangle
            %calculate the correlation coefficient between the two cells
            corr_matrix(i,j,1 + 2*(t-1)) = corr(spontFmFneu(i,:)', spontFmFneu(j,:)');
            corr_matrix(i,j,2*t) = corr(stimFmFneu(i,:)', stimFmFneu(j,:)');
        end
    end

    %%% add metadata
    data_table{dataRows,1} = (1:numCells)';
    data_table(dataRows,2) = mouse_metadata(1,'mouse_name');
    data_table(dataRows,3) = mouse_metadata(1,'sex');
    data_table(dataRows,4) = mouse_metadata(1,'order');
    data_table(dataRows,5) = mouse_metadata(1,'treatment');
    data_table(dataRows,6) = repmat({timepoint},[numCells 1]);
    data_table(dataRows,7) = mouse_metadata(1,'genotype');

    %%% add cell type (to be added later)
    data_table(dataRows,8) = repmat({'unknown'},[numCells 1]);

    %now loop for every cell
    for c = 1:numCells
        cell = cellRows(c);

        %%% add cell position
        cell_x = double(stat{cell}.med(2));
        cell_y = double(stat{cell}.med(1));
        data_table{dataRows(c),9} = cell_x;
        data_table{dataRows(c),10} = cell_y;

        %%% add cell distance to nearest stall/plaque (if there are any)
        data_table{dataRows(c),11} = results.plaque_dist(cell_y,cell_x); %distance to nearest plaque
        num_stall = table2array(mouse_metadata(end,'num_stall'));
        if num_stall~=size(mouse_metadata,1)
            d = nan;
        else
            stall_x = table2array(mouse_metadata(:,'stall_x'));
            stall_y = table2array(mouse_metadata(:,'stall_y'));
            stall_z = table2array(mouse_metadata(:,'stall_z'));
            dx = cell_x - stall_x;
            dy = cell_y - stall_y;
            dz = results.z_pos - stall_z;
            d = sqrt(dx.^2 + dy.^2 + dz.^2); %distance to all detected stalls
            d = min(d); 
        end
        data_table{dataRows(c),12} = d; %distance to nearest stall

        %%% add gcamp brightness and overall avtivity level
        data_table{dataRows(c),13} = median(stat{cell}.lam); %median baseline gcamp brightness?
        data_table{dataRows(c),14} = sum(spks(cell,spontFrames),2)/spontMinutes; %spontaneous spikes per minute
        data_table{dataRows(c),15} = sum(spks(cell,spontFrames)>0,2)/spontMinutes; %spontaneous transients per minute
        data_table{dataRows(c),16} = sum(spks(cell,stimFrames),2)/spontMinutes; %stimulated spikes per minute
        data_table{dataRows(c),17} = sum(spks(cell,stimFrames)>0,2)/spontMinutes; %stimulated transients per minute

        %%% organize stimulus-evoked spike data
        stimulus_data = nan(numFramesPerTrial,numRepetitions,numConditions);
        for rep = 1:numRepetitions
            for cond = 1:numConditions
                trial = cond + (rep-1)*numConditions;
                frame_inds = (stimFrames(1) + (trial-1)*numFramesPerTrial):(stimFrames(1) - 1 + trial*numFramesPerTrial);

                %get info about the current experimental condition
                actual_cond = vstruct(t).order(rep,cond);

                %get spikes traces for all ROIs for this condition
                stimulus_data(:,rep,cond) = spks(cell,frame_inds);
            end
        end
        stim_spks = mean(stimulus_data(stimOnFrames,:,:),1); %increase in spike rate (spikes/frame) over baseline
        stim_spks = permute(stim_spks,[2 3 1]); %[rep, cond] %mean spike rate increase
        mean_stim_spks = mean(stim_spks(repetitionsToAverage,1:12),1,'omitnan'); %[1 cond]
        std_stim_spks = std(stim_spks(repetitionsToAverage,1:12),0,1,'omitnan'); %[1 cond]

        %%% calculate orientation selectivity
        OSI_complex = sum(mean_stim_spks.*(exp(2*1i*stimulus_angles)))./sum(abs(mean_stim_spks));
        data_table{dataRows(c),18} = abs(OSI_complex); %OSI: magnitude of vector sum (displacement/distance)
        data_table{dataRows(c),19} = mod(rad2deg(angle(OSI_complex)), 180); % preferred orientation, in deg (angle of vector sum)
        %calculate circular standard deviation
        R = abs(sum(mean_stim_spks.*exp(2*1i*stimulus_angles)))/sum(abs(mean_stim_spks));
        circular_std = sqrt(-2*log(R)); %in radians for doubled angles
        circular_std_orientation = circular_std/2; %convert back to orientation space
        data_table{dataRows(c),20} = rad2deg(circular_std_orientation);
    
        %%% calculate direction selectivity
        DSI_complex = sum(mean_stim_spks.*(exp(1i*stimulus_angles)))./sum(abs(mean_stim_spks));
        data_table{dataRows(c),21} = abs(DSI_complex); %DSI
        data_table{dataRows(c),22} = mod(rad2deg(angle(DSI_complex)), 360);
        R = abs(sum(mean_stim_spks.*exp(1i*stimulus_angles)))/sum(abs(mean_stim_spks));
        circular_std_direction = sqrt(-2*log(R)); % standard deviation in radians
        data_table{dataRows(c),23} = rad2deg(circular_std_direction);

        %%% calculate fischer information
        dMu_ds = diff(mean_stim_spks);  % Numerical derivative of mean response w.r.t. stimulus
        dMu_ds(12) = mean_stim_spks(1)-mean_stim_spks(12); %circle back to the start
        FI = (dMu_ds.^2)./std_stim_spks;  % FI per stimulus
        FI = mean(FI); % average across stimulus to get average Fisher Information
        data_table{dataRows(c),24} = FI; 
        data_table{dataRows(c),25} = mean(std_stim_spks); %variance component
        data_table{dataRows(c),26} = mean(dMu_ds.^2); %activity component
    end
end


%% create images of capillary stall distance
if num_stall==size(mouse_metadata,1) %if there are any stalls
    %load stabilized ch2
    stall_timepoints = table2array(mouse_metadata(:,'stall_timepoint'));
    % stall_x = table2array(mouse_metadata(:,'stall_x'));
    % stall_y = table2array(mouse_metadata(:,'stall_y'));
    % stall_z = table2array(mouse_metadata(:,'stall_z'));
    stall_start_frame = table2array(mouse_metadata(:,'stall_start_frame'));
    stall_stop_frame = table2array(mouse_metadata(:,'stall_stop_frame'));

    %create distace map
    isize = mdata(1).imagesize.iseries;
    stall_dist = nan(isize);
    [xmap,ymap] = meshgrid(1:isize,1:isize); %x/y-maps based on isometric volume
    pixelSize = getPixelSize('3', datetime('7/1/2018'), 'Zeiss 20x', 2); %0.3907 um/pixel @ 1024x1024 & 1.5zoom
    mpp_xy = pixelSize(1)*(1024/isize(1)); %microns per pixel, x/y axis
    mpp_z = 1; %microns per pixel, z axis
    for i = 1:num_stall
        dx = (xmap - stall_x(i))*mpp_xy;
        dy = (ymap - stall_y(i))*mpp_xy;
        dz = (results.z_pos - stall_z(i))*mpp_z;
        d = sqrt(dx.^2 + dy.^2 + dz.^2); %distance to all detected stalls
        tmp = stall_dist;
        tmp(:,:,2) = d;
        stall_dist = min(tmp,[],3,'omitnan');
    end
    figure()
    imshow(stall_dist/80)
    cm = hot;
    colormap(flipud(cm))
    cb = colorbar;
    cb.Ticks = 0:0.25:1;
    cb.TickLabels = {'0','20','40','60','80+'};
    ax = gca;
    % ax.YDir = 'normal';
    % axis off
    saveas(gcf,fullfile(analysis_folder,'stall_distance.svg'))
    
    %%%create example image series of stall and gcamp
    %load ch2, ch3 of timepoint 1
    i = 1;
    filename = fullfile(analysis_folder,num2str(i),'MC functional','CH2_functional.tif');
    gcamp_I = read_file(filename);
    filename = fullfile(analysis_folder,num2str(i),'MC functional','CH4_functional.tif');
    rgb = zeros([isize 3]);
    vessel_I = read_file(filename);
    frames = [stall_start_frame(i)-4 stall_start_frame(i) stall_start_frame(i)+4 stall_stop_frame(i)];
    for f = 1:4
        tmp = double(vessel_I(:,:,frames(f)));
        tmp = tmp-prctile(tmp(:),1);
        tmp = tmp/prctile(tmp(:),99);
        rgb(:,:,1) = tmp;
        tmp = double(gcamp_I(:,:,frames(f)));
        tmp = tmp-prctile(tmp(:),1);
        tmp = tmp/prctile(tmp(:),99);
        % rgb(:,:,2) = tmp;
        figure()
        imshow(rgb)
        saveas(gcf,fullfile(analysis_folder,['stall_frame_' num2str(f) '.svg']))
    end
    
end


%% save data
filename = fullfile(analysis_folder,'data_table.mat');
save(filename,'data_table');


