

%% settings
analysis_name = '22_5_9_ly6g_test';
exp_dirs{1} = 'C:\Users\misaa\Desktop\2022-05-08 fulltest_beforeLy6g_22_5_8'; %e.g. baseline
exp_dirs{2} = 'C:\Users\misaa\Desktop\2022-05-09 fulltest_afterLy6g_22_5_9'; %e.g. post-injection


%% pre-process imaging
%change settings if desired
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels
opts_tiff.append = false; opts_tiff.big = true; opts_tiff.message = false;

%create directories
num_dirs = length(exp_dirs);
base_dir = fileparts(exp_dirs{1});
save_dir = fullfile(base_dir,[datestr(now,'yyyy-mm-dd HH-MM ') analysis_name]);
assert(~isfolder(save_dir),'Analysis folder already exists');
mkdir(save_dir);
for d = 1:num_dirs
    dstr = num2str(d);
    mkdir(fullfile(save_dir,dstr));
    mkdir(fullfile(save_dir,dstr,'suite2p'));
end

%load metadata for each experiment
for d = 1:num_dirs
    mdataname = ls(fullfile(exp_dirs{d},'*metadata.mat'));
    vsname = ls(fullfile(exp_dirs{d},'*vs.mat'));
    load(fullfile(exp_dirs{d},mdataname),'metadata');
    mdata(d) = metadata;
    load(fullfile(exp_dirs{d},vsname),'vs');
    vstruct(d) = vs;
end

%set unmixing coefficients
gc = [0    0.3  0.6  0    0    0.05 0.05 0    ]; %gcamp6s
yc = [0    0.05 0.3  0    0    0.1  0.5  0.05 ]; %eyfp
rc = [0    0    0    0.04 0    0.04 0.38 0.54 ]; %mruby2
tc = [0    0    0    1    0    0    0    0    ]; %texas red (maybe some in cham ch4?)
nc = [0.2  0.1  0.1  0.1  0.2  0.1  0.1  0.1  ]; %pmt noise
unmixingCoeffs = [gc', yc', rc', tc', nc'];


%pre-process imaging for each directory
for d = 2:num_dirs
    dstr = num2str(d);
    fprintf('Pre-processing complete.\n\n');
    
    %concatentate all functional imaging together into "functional" folder
    fprintf('concatenating files...')
    filenames = dir(fullfile(exp_dirs{d},'stim','*.tif'));
    ns = length(filenames);
    filenames(ns+1) = dir(fullfile(exp_dirs{d},'spont*.tif'));
    functional = read_file(fullfile(filenames(ns+1).folder,filenames(ns+1).name));
    [iheight, iwidth, nspf] = size(functional);
    nspf = nspf/4; %number of spontaneous frames
    tmp = read_file(fullfile(filenames(1).folder,filenames(1).name));
    nsf = size(tmp,3)/4; %number of stim frames
    functional(:,:,(1+(nspf*4)):((nspf*4)+(ns*nsf*4))) = nan; %pre-allocate space for stim frames
    for f = 1:ns
        tmp = read_file(fullfile(filenames(f).folder,filenames(f).name));
        functional(:,:,1+(nspf*4)+((f-1)*nsf*4):(nspf*4)+(f*nsf*4)) = tmp;
    end
    imageName = fullfile(save_dir,dstr,'functional.tif');
    saveastiff(functional,imageName,opts_tiff);
    fprintf('done.\n')
    
    %motion correct it 
    fprintf('Starting motion correction:\n');
    run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
    
    %save it in a "suite2p" folder (ch2 only?)
    imageName = fullfile(save_dir,dstr,'MC functional','CH2_functional.tif');
    copyfile(imageName,fullfile(save_dir,dstr,'suite2p','CH2_functional.tif'))
    
    %load imaging for color unmixing (create 8ch mixed)
    fprintf('Starting fluorophore unmixing...');
    sats = read_file(fullfile(exp_dirs{d},'1030_00001.tif'));
    ncf = size(sats,3)/4;
    sats = reshape(sats,[iheight iwidth 4 ncf]);
    cham = reshape(functional(:,:,1:(4*ncf)),[iheight iwidth 4 ncf]);
    mixed = mean(cham,4);
    mixed(:,:,5:8) = mean(sats,4);
    saveastiff(mixed,fullfile(save_dir,dstr,'mixed.tif'));
    mixed = uint16(mixed);
    
    %save metadata
    mdata(d).num_stims = ns;
    mdata(d).num_stim_frames = nsf;
    mdata(d).num_spont_frames = nspf;
    mdata(d).num_color_frames = ncf;
    
    %test out histogram matching
    % figure()
    % for c = 1:8
    %    subplot(8,2,1+2*(c-1))
    %    histogram(mixed(:,:,c))
    %    xlim([0 65535])
    % end
    % ref = mixed(:,:,3);
    % scale = 1/(double(max(ref,[],'all'))/65535);
    % ref = uint8(255*scale*(double(ref)/65535));
    % for c = 1:8
    %    mixed(:,:,c) = imhistmatch(uint8(255*(double(mixed(:,:,c))/65535)),ref);
    %    subplot(8,2,2*c)
    %    histogram(mixed(:,:,c))
    %    xlim([0 255])
    % end
    % figure()
    % imshow(ref)
    % figure()
    % imshow(mixed(:,:,3))

    %unmix fluorophores from channels
    [~, ~, numChannels, numFrames, numSlices] = size(mixed);
    mixed = permute(mixed,[3,1,2,4,5]);               % Order matrix [channel, y, x, frame, slice]
    mixed = reshape(mixed,numChannels,[]);                      % Flatten matrix
    numPts = size(mixed,2);
    numFluorophores = size(unmixingCoeffs,2);
    Unmixed = nan(numFluorophores,numPts);
    for i = 1:numPts
        Unmixed(:,i) = cast(lsqnonneg(unmixingCoeffs,double(mixed(:,i))),'int16');  % Use nonnegative least squares to solve system
    end
    Unmixed = reshape(Unmixed,[numFluorophores,iheight,iwidth,numFrames,numSlices]); % Expand matrix to 5D
    Unmixed = permute(Unmixed,[2,3,1,4,5]);               % volume matrix re-ordered to [y, x, channel, frame, slice]
    saveastiff(Unmixed,fullfile(save_dir,dstr,'unmixed.tif'),opts_tiff);
    
    %save png of eyfp/mruby yellow/red image 
    colorIm = zeros([iheight iwidth 3]);
    colorIm(:,:,1) = mean(Unmixed(:,:,[2 3]),3);
    colorIm(:,:,2) = Unmixed(:,:,2)/2;
    imwrite(uint8(100*colorIm/prctile(colorIm,98,'all')),fullfile(save_dir,dstr,'YellowRed.tif'));
    fprintf('done.\n\n');
end
save(fullfile(save_dir,'exp_metadata.mat'),'mdata','vstruct');
fprintf('Pre-processing complete.\n\n');

