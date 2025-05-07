

%% settings
analysis_name = '2024-09-24 APP Control';
exp_dirs{1} = 'C:\Users\misaa\Desktop\Awake Cortex Tests\2024-09-24 full20x_F_APP_before_control_mouse_1R2L_4294'; %e.g. baseline
exp_dirs{2} = 'C:\Users\misaa\Desktop\Awake Cortex Tests\2024-09-25 full20x_F_APP_after_control_mouse_1R2L_4294'; %e.g. post-injection
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
end

%load metadata for each experiment
for d = 1:num_dirs
    mdataname = dir(fullfile(exp_dirs{d},'*metadata.mat'));
    mdataname = mdataname(1).name;
    vsname = dir(fullfile(exp_dirs{d},'*vs.mat'));
    vsname = vsname(1).name;
    load(fullfile(exp_dirs{d},mdataname),'metadata');
    mdata(d) = metadata;
    load(fullfile(exp_dirs{d},vsname),'vs');
    vstruct(d) = vs;
end


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

    %Extract blue channel from plaque stacks
    filename = dir(fullfile(exp_dirs{d},'plaques*.tif'));
    plaques = read_file(fullfile(filename.folder,filename.name));
    blueplaques = plaques(:,:,1:4:end);
    saveastiff(blueplaques,fullfile(save_dir,dstr,'blue_plaques.tif'));
    greenplaques = plaques(:,:,2:4:end);
    saveastiff(greenplaques,fullfile(save_dir,dstr,'green_plaques.tif'));

    %store metadata
    mdata(d).num_stims = ns;
    mdata(d).num_stim_frames = nstf;
    mdata(d).num_spont_frames = nspf;
    fprintf('done.\n\n');
end


%% save metadata
save(fullfile(save_dir,'exp_metadata.mat'),'mdata','vstruct');
fprintf('Pre-processing complete.\n\n');

