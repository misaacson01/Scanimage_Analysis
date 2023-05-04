

%% settings
analysis_name = '22_5_8_9_ly6g_test';
exp_dirs{1} = 'C:\Users\misaa\Desktop\2022-06-30 full20x_APP_before_ly6g_mouse L'; %e.g. baseline
exp_dirs{2} = 'C:\Users\misaa\Desktop\2022-07-01 full20x_APP_after_ly6g_mouse_L'; %e.g. post-injection
register_all_to_first = true;

%% pre-process imaging
%change settings if desired
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels
opts_tiff.append = false; opts_tiff.big = true; opts_tiff.message = false;
[optimizer,metric] = imregconfig('multimodal');

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
    mkdir(fullfile(save_dir,dstr,'suite2p_comp'));
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

%set cell-type unmixing coefficients
gc = [0     0.3   0.6   0     0     0.05  0.05  0    ]; %gcamp6s
yc = [0     0.05  0.3   0     0     0.1   0.5   0.05 ]; %eyfp
rc = [0     0     0     0.04  0     0.04  0.38  0.54 ]; %mruby2
tc = [0     0     0     1     0     0     0     0    ]; %texas red (maybe some in sats ch4?)
nc = [0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125]; %pmt noise
unmixingCoeffs = [gc', yc', rc', tc', nc'];

%new dilution coefficients
gc = [0     0.33  0.24  0.06  0     0.15  0.17  0.05 ]; %gcamp6s (
yc = [0     0     0     0     0     0.15   0.75   0.1 ]; %eyfp
rc = [0     0     0     0     0     0     0.28  0.72 ]; %mruby2
tc = [0     0     0     1     0     0     0     0    ]; %texas red (maybe some in sats ch4?)
nc = [0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125]; %pmt noise
unmixingCoeffs = [gc', yc', rc', tc', nc'];

%set 780 nm unmixing coefficients
%%%%%%%%%%%fix these
% bc = [1    0     0     0    ];
% gc = [0    0.35  0.65  0    ];
% yc = [0    0.15  0.85  0    ];
% tc = [0    0.05  0.03  0.65 ];
% nc = [0.25 0.25  0.25  0.25 ];
% unmixingCoeffs_780 = [bc', gc', yc', tc', nc'];

%set 920 nm unmixing coefficients
% unmixingCoeffs_920 = unmixingCoeffs(1:4,:)/repmat(sum(unmixingCoeffs(1:4,:)),[4 1]);

%% pre-process imaging for each directory
for d = 1:num_dirs
    dstr = num2str(d);
    fprintf(['Starting pre-processing for exp ' dstr ' (of ' num2str(num_dirs) '):\n']);
   
    %%%% functional imaging %%%%
    %concatentate all functional imaging together into "functional" folder
    fprintf('Concatenating files (this can take some time)... ')
    filenames = dir(fullfile(exp_dirs{d},'stim','*.tif'));
    ns = length(filenames);
    filenames(ns+1) = dir(fullfile(exp_dirs{d},'spont*.tif'));
    functional = read_file(fullfile(filenames(ns+1).folder,filenames(ns+1).name));
    [iheight, iwidth, nspf] = size(functional);
    nspf = nspf/4; %number of spontaneous frames
    
    tmp = read_file(fullfile(filenames(1).folder,filenames(1).name));
    nstf = size(tmp,3)/4; %number of stim trial frames
    functional(:,:,(1+(nspf*4)):((nspf*4)+(ns*nstf*4))) = nan; %pre-allocate space for stim frames
    for f = 1:ns
        tmp = read_file(fullfile(filenames(f).folder,filenames(f).name));
        functional(:,:,1+(nspf*4)+((f-1)*nstf*4):(nspf*4)+(f*nstf*4)) = tmp;
    end
    nsf = (size(functional,3)/4) - nspf; %number of stim frames (all trials)
    
    imageName = fullfile(save_dir,dstr,'functional.tif');
    saveastiff(functional,imageName,opts_tiff);
    fprintf('done.\n')
    
    %motion correct it 
    fprintf('Starting motion correction:\n');
    run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
    dstr_first = num2str(d);
    template_first = read_file(fullfile(save_dir,dstr_first,'MC functional','TEMPLATE_functional.tif'));
    %%
    if d>1 && register_all_to_first
        %register functional imaging to first session
        template = read_file(fullfile(save_dir,dstr,'MC functional','TEMPLATE_functional.tif'));
        tform = imregtform(template,template_first,'affine',optimizer,metric);
        ch2_functional = read_file(fullfile(save_dir,dstr,'MC functional','CH2_functional.tif'));
        sameAsInput = affineOutputView([iheight iwidth],tform,'BoundsStyle','SameAsInput');
        for f = 1:size(ch2_functional,3)
            imwarp(ch2_functional(:,:,f),tform,'OutputView',sameAsInput);
        end
        saveastiff(ch2_functional,fullfile(save_dir,dstr,'suite2p','CH2_functional.tif'));
    else
        %copy the functional imaging to a "suite2p" folder (ch2 only for now)
        imageName = fullfile(save_dir,dstr,'MC functional','CH2_functional.tif');
        copyfile(imageName,fullfile(save_dir,dstr,'suite2p','CH2_functional.tif'));
    end
    %%
    %spectral unmixing for excitatory/inhibitory neurons
    %load images for color unmixing (create 8ch mixed)
    fprintf('Starting fluorophore unmixing...');
    sats = read_file(fullfile(exp_dirs{d},'1030_00001.tif'));
    ncf = size(sats,3)/4; %number of color frames
    sats = mean(reshape(sats,[iheight iwidth 4 ncf]),4);
    cham = mean(reshape(functional(:,:,1:(4*ncf)),[iheight iwidth 4 ncf]),4);
    if register_all_to_first
        assert(size(sats,1)==size(template_first,1),'functional imaging and color imaging are not the same size');
        %register satsuma images to first functional template
        tform = imregtform(max(sats(:,:,channel_options.chsh),[],3),template_first,'affine',optimizer,metric);
        sameAsInput = affineOutputView([iheight iwidth],tform,'BoundsStyle','SameAsInput');
        for c = 1:size(sats,3)
            imwarp(sats(:,:,c),tform,'OutputView',sameAsInput);
        end
        %register chameleon images to first functional template
        tform = imregtform(max(cham(:,:,channel_options.chsh),[],3),template_first,'affine',optimizer,metric);
        sameAsInput = affineOutputView([iheight iwidth],tform,'BoundsStyle','SameAsInput');
        for c = 1:size(cham,3)
            imwarp(cham(:,:,c),tform,'OutputView',sameAsInput);
        end
    else
        %register satsuma images to chameleon images
        tform = imregtform(max(sats(:,:,channel_options.chsh),[],3),max(cham(:,:,channel_options.chsh),[],3),'affine',optimizer,metric);
        sameAsInput = affineOutputView([iheight iwidth],tform,'BoundsStyle','SameAsInput');
        for c = 1:size(sats,3)
            sats(:,:,c) = imwarp(sats(:,:,c),tform,'OutputView',sameAsInput);
        end
    end
    mixed = cham;
    mixed(:,:,5:8) = sats;
    mixed = uint16(mixed);
    
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
    saveastiff(Unmixed,fullfile(save_dir,dstr,'celltypes_unmixed.tif'),opts_tiff);
    
    %save green, red, and yellow .tifs separately (added for ROImaker pipeline)
    saveastiff(Unmixed(:,:,1),fullfile(save_dir,dstr,'gcamp_unmixed.tif'),opts_tiff);
    saveastiff(Unmixed(:,:,2),fullfile(save_dir,dstr,'eyfp_unmixed.tif'),opts_tiff);
    saveastiff(Unmixed(:,:,3),fullfile(save_dir,dstr,'mruby_unmixed.tif'),opts_tiff);
    
    %save png of eyfp/mruby yellow/red image 
    colorIm = zeros([iheight iwidth 3]);
    colorIm(:,:,1) = mean(Unmixed(:,:,[2 3]),3);
    colorIm(:,:,2) = Unmixed(:,:,2)/2;
    imwrite(uint8(100*colorIm/prctile(colorIm,98,'all')),fullfile(save_dir,dstr,'YellowRed.tif'));
    
    %combining unmixed red/yellow images into ch2 functional (for suite2p manual ROIs)
    %get the mean ch2 image in green
    ch2_functional = read_file(fullfile(save_dir,dstr,'MC functional','CH2_functional.tif'));
    ch2_functional_mean = mean(ch2_functional,3,'omitnan');
    %get the red/yellow image
    anatomical_mean = mean(Unmixed(:,:,[2 3]),3,'omitnan');

    %run a GUI to manually register
    registerApp = manualRegistration(ch2_functional_mean,anatomical_mean);
    waitfor(registerApp,'editing','off');
    rotation = registerApp.rotation;
    translation = registerApp.translation;
    registerApp.delete
    
    % rotate/translate the red/yellow image
    tform = rigidtform2d(rotation,translation);
    centerOutput = affineOutputView(size(anatomical_mean),tform,"BoundsStyle","CenterOutput");
    anatomical_mean = imwarp(anatomical_mean,tform,"OutputView",centerOutput);
    
    % create a composite image for suite2p
    ch2_functional_composite = (ch2_functional/2) + repmat(anatomical_mean/2,[1 1 num_ch2_frames]);
    imageName = fullfile(save_dir,dstr,'suite2p_comp','CH2_functional_comp.tif');
    saveastiff(ch2_functional_composite,imageName,opts_tiff);

%     %%%% spectral unmixing of stacks for plaques, vessels %%%%
%     %plaques stack
%     tmp = read_file(fullfile(exp_dirs{d},'plaques_00001.tif'));
%     numChannels = 4;
%     [stackheight, stackwidth, numFrames] = size(tmp);
%     numFrames = numFrames/4;
%     numSlices = 1;
%     mixed = uint16(reshape(tmp,[stackheight stackwidth numChannels numFrames]));
%     mixed = permute(mixed,[3,1,2,4,5]);               % Order matrix [channel, y, x, frame, slice]
%     mixed = reshape(mixed,numChannels,[]);                      % Flatten matrix
%     numPts = size(mixed,2);
%     numFluorophores = size(unmixingCoeffs_780,2);
%     Unmixed = nan(numFluorophores,numPts);
%     for i = 1:numPts
%         Unmixed(:,i) = cast(lsqnonneg(unmixingCoeffs_780,double(mixed(:,i))),'int16');  % Use nonnegative least squares to solve system
%     end
%     Unmixed = reshape(Unmixed,[numFluorophores,stackheight,stackwidth,numFrames,numSlices]); % Expand matrix to 5D
%     Unmixed = permute(Unmixed,[2,3,4,1,5]);               % volume matrix re-ordered to [y, x, frame, channel, slice]
%     saveastiff(Unmixed(:,:,:,1),fullfile(save_dir,dstr,'plaques_unmixed.tif'),opts_tiff);
%     %%
%     %vessels stack
%     tmp = read_file(fullfile(exp_dirs{d},'vessels_00001.tif'));
%     numChannels = 4;
%     [stackheight, stackwidth, numFrames] = size(tmp);
%     numFrames = numFrames/4;
%     numSlices = 1;
%     mixed = uint16(reshape(tmp,[stackheight stackwidth numChannels numFrames]));
%     mixed = permute(mixed,[3,1,2,4,5]);               % Order matrix [channel, y, x, frame, slice]
%     mixed = reshape(mixed,numChannels,[]);                      % Flatten matrix
%     numPts = size(mixed,2);
%     numFluorophores = size(unmixingCoeffs_920,2);
%     Unmixed = nan(numFluorophores,numPts);
%     for i = 1:numPts
%         Unmixed(:,i) = cast(lsqnonneg(unmixingCoeffs_920,double(mixed(:,i))),'int16');  % Use nonnegative least squares to solve system
%     end
%     Unmixed = reshape(Unmixed,[numFluorophores,stackheight,stackwidth,numFrames,numSlices]); % Expand matrix to 5D
%     Unmixed = permute(Unmixed,[2,3,4,1,5]);               % volume matrix re-ordered to [y, x, frame, channel, slice]
%     saveastiff(Unmixed(:,:,:,4),fullfile(save_dir,dstr,'vessels_unmixed.tif'),opts_tiff);
%     %%
%     %%%% spectral unmixing of functional gcamp, instead of just using ch2 %%%%
%     numFrames = nspf+nstf;
%     mixed = nan([iheight, iwidth, 4, numFrames]);
%     for c = 1:4
%         imageName = fullfile(save_dir,dstr,'MC functional',['CH' num2str(c) '_functional.tif']);
%         mixed(:,:,c,:) = uint16(read_file(imageName));
%     end
%     mixed = permute(mixed,[3,1,2,4,5]);               % Order matrix [channel, y, x, frame, slice]
%     mixed = reshape(mixed,numChannels,[]);                      % Flatten matrix
%     numPts = size(mixed,2);
%     numFluorophores = size(unmixingCoeffs_920,2);
%     Unmixed = nan(numFluorophores,numPts);
%     for i = 1:numPts
%         Unmixed(:,i) = cast(lsqnonneg(unmixingCoeffs_920,double(mixed(:,i))),'int16');  % Use nonnegative least squares to solve system
%     end
%     Unmixed = reshape(Unmixed,[numFluorophores,iheight,iwidth,numFrames,numSlices]); % Expand matrix to 5D
%     Unmixed = permute(Unmixed,[2,3,4,1,5]);               % volume matrix re-ordered to [y, x, frame, channel, slice]
%     saveastiff(Unmixed(:,:,:,2),fullfile(save_dir,dstr,'gcamp_unmixed.tif'),opts_tiff);
    
    %store metadata
    mdata(d).num_stims = ns;
    mdata(d).num_stim_frames = nstf;
    mdata(d).num_spont_frames = nspf;
    mdata(d).num_color_frames = ncf;
    fprintf('done.\n\n');
    %%
end

%% save metadata
save(fullfile(save_dir,'exp_metadata.mat'),'mdata','vstruct');
fprintf('Pre-processing complete.\n\n');

