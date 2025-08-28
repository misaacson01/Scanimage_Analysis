
%set analysis folder
analysis_folder = 'C:\Users\misaa\Desktop\Awake Cortex Tests\2025-05-01 14-14 2024-09-24 APP Control';

%check if part2 has already been run
if isfile(fullfile(analysis_folder,'plaque_results.mat'))
    answer = questdlg('Part 2 has already been run on this dataset. Overwrite?', ...
                         'Overwite Question', 'yes', 'no (exit)', 'no (exit)');
    if ~strcmp(answer,'yes')
        fprintf('script exited\n')
        return
    end
end
fprintf('running Part2: pre-processing anatomical imaging data\n')


%% load imaging data
wfig = waitbar(0,'loading data');
% load experiment metadata
filename = fullfile(analysis_folder,'exp_metadata.mat');
assert(isfile(filename),'exp_metadata.mat file not found in analysis folder')
load(filename);

%load plaque volume (both blue and green channels) from timepoint 2
blue_plaques = read_file(fullfile(analysis_folder,'2','blue_plaques.tif'));
green_plaques = read_file(fullfile(analysis_folder,'2','green_plaques.tif'));
[h,w,nf] = size(green_plaques);

%load functional imaging plane (ch2 from timepoint 2)
green_functional = read_file(fullfile(analysis_folder,'2','suite2p','CH2_functional.tif'));
green_functional = mean(green_functional(:,:,1:mdata(2).num_spont_frames),3,'omitnan');
[h_f,w_f] = size(green_functional);
fprintf('data loaded\n')

%% register functional plane assuming different zoom factors to find best fit
search_zoom_ind = 2; %start with only seaching for zoom 1.333, since this is the most likely situation
filt_sigma = 3;
zoom_factors = [1 1.333 1.5 1.666 2];
num_zooms = length(zoom_factors);
best_correlations = zeros(1, num_zooms);
best_slices = zeros(1, num_zooms);
best_transforms = cell(1, num_zooms);
%list of zoom factors to test:
%1 -- same zoom as stack
%1.3333x zoom relative to stack (e.g. 2.0/1.5)
%1.5x zoom relative to stack (e.g. 1.5/1.0)
%1.666x zoom relative to stack (e.g. 2.5/1.5)
%2x zoom relative to stack (e.g. 2.0/1.0)

[optimizer, metric] = imregconfig('multimodal'); % Suitable for images from different sensors or modalities
optimizer.MaximumIterations = 300; % Increase iterations for better convergence
optimizer.InitialRadius = 0.009; % Tune optimizer parameters
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;

R = nan(nf,num_zooms);
if nf>180 && nf<220 
    startframe = 50;
    stopframe = 150;
elseif nf>480 && nf<520 
    startframe = 100;
    stopframe = 300;
else
    error('unexpected number of frames')
end

%manually set start/stop frames
% startframe = 204;
% stopframe = 206;

regwarning_id = 'images:imregcorr:weakPeakCorrelation';
warning('off',regwarning_id) %turn off warning for poor registration (many are expected)
for z = search_zoom_ind
    %create images of functional for current zoom factor
    cur_func_image = imresize(green_functional,[h w]/zoom_factors(z));
    cur_func_image = imgaussfilt(cur_func_image,filt_sigma);
    cur_func_image = cur_func_image - prctile(cur_func_image(:),1);
    cur_func_image = cur_func_image/prctile(cur_func_image(:),99);
    for i = startframe:stopframe 
        if isvalid(wfig)
            waitbar(0.9*((i-startframe)/(stopframe-startframe)),wfig,['finding imaging plane within volume (z' num2str(z) ' of ' num2str(num_zooms) ')'])
        else
            wfig = waitbar(0.9*((i-startframe)/(stopframe-startframe)),['finding imaging plane within volume (z' num2str(z) ' of ' num2str(num_zooms) ')']);
        end
        cur_stack_image = double(green_plaques(:,:,i));
        cur_stack_image = imgaussfilt(cur_stack_image,filt_sigma);
        cur_stack_image = cur_stack_image - prctile(cur_stack_image(:),1);
        cur_stack_image = cur_stack_image/prctile(cur_stack_image(:),99);
        cur_stack_image = imhistmatchn(cur_stack_image,cur_func_image);
        registered_func_image = imregister(cur_func_image, cur_stack_image, 'rigid', optimizer, metric);
        R(i,z) = corr2(registered_func_image,cur_stack_image);
    end
end
warning('on',regwarning_id)

%plot correlation by zoom/depth
figure()
plot(R)
legend({'z1.0','z1.33','z1.5','z1.66'})
xlabel('stack frame number')
ylabel('correlation coefficient')
[best_corrs,best_frames] = max(R);
[~,best_zoom_ind] = max(best_corrs);
best_zoom = zoom_factors(best_zoom_ind);
best_frame = best_frames(best_zoom_ind);
title(['best frame: ' num2str(best_frame) ' (zoom ' num2str(best_zoom) ')'])
saveas(gcf,fullfile(analysis_folder,'stack_location_plot.png'))
fprintf('functional plane correlated to plaque volume\n')


%% create overlay for manual verification
if isvalid(wfig)
    waitbar(0.91,wfig,'creating imaging plane/stack overlay for manual verification')
else
    wfig = waitbar(0.91,'creating imaging plane/stack overlay for manual verification');
end

cur_func_image = imresize(green_functional,[h w]/best_zoom);
cur_func_image = cur_func_image - prctile(cur_func_image(:),1);
cur_func_image = cur_func_image/prctile(cur_func_image(:),99);

cur_stack_image = double(green_plaques(:,:,best_frame));
cur_stack_image = cur_stack_image - prctile(cur_stack_image(:),1);
cur_stack_image = cur_stack_image/prctile(cur_stack_image(:),99);
cur_stack_image = imhistmatchn(cur_stack_image,cur_func_image);

ref_stack = imref2d(size(cur_stack_image));
ref_func = imref2d(size(cur_func_image));
tform = imregtform(cur_func_image, ref_func, cur_stack_image, ref_stack, 'rigid', optimizer, metric);
registered_func_image = imwarp(cur_func_image, tform, 'OutputView', ref_stack);
%registered_func_image = imregister(cur_func_image, cur_stack_image, 'rigid', optimizer, metric);

comb_rgb = registered_func_image;
comb_rgb(:,:,3) = cur_stack_image;
stack_rgb = repmat(cur_stack_image,[1 1 3]);
stack_rgb(:,:,1:2) = 0;

figure()
subplot(1,3,1)
imshow(stack_rgb)
title('stack frame (blue)');
subplot(1,3,2)
imshow(comb_rgb)
title('functional plane (red)');
subplot(1,3,3)
imshow(stack_rgb)
hold on
[h_rf, w_rf] = size(cur_func_image); 
[Xb_original, Yb_original] = meshgrid(1:w_rf, 1:h_rf);
original_points = [Xb_original(:), Yb_original(:)]; % N x 2 matrix of [x, y]
registered_points = transformPointsForward(tform, original_points);
x_pos = reshape(registered_points(:,1), size(Xb_original));
y_pos = reshape(registered_points(:,2), size(Yb_original));
x_pos = imresize(x_pos,[h_f w_f],"bilinear");
y_pos = imresize(y_pos,[h_f w_f],"bilinear");
z_pos = best_frame;
scatter(x_pos,y_pos,3,'green');
title('functional plane position')
saveas(gcf,fullfile(analysis_folder,'stack_location_overlay.png'))

results.x_pos = x_pos;
results.y_pos = y_pos;
results.z_pos = z_pos;
fprintf('best fit displayed\n')


%% perform plaque binarization
if isvalid(wfig)
    waitbar(0.91,wfig,'binarizing plaque volume')
else
    wfig = waitbar(0.91,'binarizing plaque volume');
end
switch best_zoom
    case 1
        stack_zoom = 1.5;
    case 1.333
        stack_zoom = 1.5;
    case 1.5
        stack_zoom = 1.0;
    case 1.666
        stack_zoom = 1.5;
    case 2
        stack_zoom = 1.0;
end
%find scale factor and subtract non-blue labelling (GCaMP) from blue
filt_sigma = 5;
pixelSize = getPixelSize('3', datetime('7/1/2018'), 'Zeiss 20x', stack_zoom); %0.3907 um/pixel @ 1024x1024 & 1.5zoom
mpp_xy = pixelSize(1)*(1024/w); %microns per pixel, x/y axis
mpp_z = 1; %microns per pixel, z axis

%subtract green gcamp labels from blue volume assuming linear scale factor
p = polyfit(single(green_plaques(:)), single(blue_plaques(:)), 1); % 1st-order polynomial fit (linear)
alpha = p(1); % slope is the scaling factor
plaques = single(blue_plaques) - alpha * single(green_plaques);

%make stack isometric
plaques = imresize(plaques,round([h w]*(mpp_xy/mpp_z)));
plaques = max(plaques, 0); % clip negative values
plaques = imgaussfilt3(plaques,filt_sigma);

%binarize plaques using 5*sd threshold above bottom 95%
prc95 = prctile(plaques(:),95);
bottom95 = plaques(plaques<prc95);
mean_bottom95 = mean(bottom95(:));
std_bottom95 = std(bottom95(:));
threshold = mean_bottom95+(5*std_bottom95);
binary_vol = plaques>threshold;

%find connected components and remove smallest
cc = bwconncomp(binary_vol, 26); % 26-connectivity for 3D
stats = regionprops3(cc, 'Volume', 'Centroid');
min_plaque_radius = 5/mpp_z; %minimum radius, in pixels
min_plaque_volume = (4/3)*pi*(min_plaque_radius^3);
results.rawstats = stats;
if isempty(stats)
    disp('no plaques detected');
    plaques_bin = false(size(corrected_vol));
else
    plaques_bin = ismember(labelmatrix(cc), find(stats.Volume>=min_plaque_volume));
end
stats = stats(stats.Volume >= min_plaque_volume, :);
results.stats = stats;
results.min_plaque_radius = min_plaque_radius;
results.min_plaque_volume = min_plaque_volume;
results.plaques_bin = plaques_bin;

%visualize plaque binarization
volshow(plaques_bin)
plaques_b = plaques/(mean_bottom95+(10*std_bottom95));
plaques_rb = plaques_b/5;
plaques_rb(plaques_bin) = 1;
ready_to_save = true;
if isfile(fullfile(analysis_folder,'plaques_bin_check.tif'))
    delete(fullfile(analysis_folder,'plaques_bin_check.tif'));
end
try
    saveastiff([plaques_b plaques_rb],fullfile(analysis_folder,'plaques_bin_check.tif'));
catch
    warning('couldn''t save "plaques_bin_check.tif"; it might already exist and matlab can''t delete it')
end
fprintf('plaque volume binarized\n')


%% calculate distance to nearest plaque for every pixel of the cell plane
[h_iso, w_iso] = size(plaques_bin,[1 2]);
plaque_dist = nan(h_f,w_f);
zdist = mpp_z:mpp_z:(nf*mpp_z);
zdist = zdist - (best_frame*mpp_z);
zdist = permute(zdist,[1 3 2]);
zdist = repmat(zdist,[h_iso w_iso 1]);

tmp = double(plaques_bin);
tmp(tmp==0) = nan;
minzdist = abs(zdist).*tmp;
minzdist = min(minzdist,[],3); %z-distance to closest plaque
[xmap,ymap] = meshgrid(linspace(1,w,w_iso),linspace(1,h,h_iso)); %x/y-maps based on isometric volume

for x_ind = 1:length(x_pos)
    if isvalid(wfig)
        waitbar(0.92 + 0.08*(x_ind/length(x_pos)),wfig,'calculating distance to plaques')
    else
        wfig = waitbar(0.92 + 0.08*(x_ind/length(x_pos)),'calculating distance to plaques');
    end
    for y_ind = 1:length(y_pos)
    % 
        %create 3D map of distance to position
        cur_x_pos = x_pos(y_ind,x_ind);
        cur_y_pos = y_pos(y_ind,x_ind);

        %calculate distance to all plaques
        dist = sqrt((xmap-cur_x_pos).^2 + (ymap-cur_y_pos).^2 + minzdist.^2);
        plaque_dist(y_ind,x_ind) = min(dist,[],'all','omitnan');
    end
end
% plaque_dist(plaque_dist<0) = 0;
results.plaque_dist = plaque_dist;
fprintf('distance to nearest plaque calculated\n')

%visualize plaque distance maps
figure()
imshow(plaque_dist/80)
cm = hot;
colormap(flipud(cm))
cb = colorbar;
cb.Ticks = 0:0.25:1;
cb.TickLabels = {'0','20','40','60','80+'};
ax = gca;
ax.YDir = 'normal';
axis on
saveas(gcf,fullfile(analysis_folder,'plaque_distance.png'))


%% save data
delete(wfig)
save(fullfile(analysis_folder,'plaque_results.mat'),'results')
fprintf('script finished\n')
