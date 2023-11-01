%run this script after suite2p analysis has finished for both before/after
%timepoints
%
%After this script has finished, a CellLinker.tif file will be saved in the
%analysis folder, which you can then open up in imagej to manually identify
%cells between timepoints, and (ideally) identify their cell type
%
%If we stick with this kind of manual analysis, then we may want to replace
%imagej with a matlab gui that makes it easier to save cell linking info

%% settings
analysis_folder = 'C:\Users\misaa\Desktop\2022-07-03 09-18 22_5_8_9_ly6g_test';
image_size = [500 1000]; %[height, width] of output image (in pixels)

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
        
%create tif stack of before/after:
%% 1. mean functional image
mean_functional_1 = read_file(fullfile(analysis_folder,'1','MC functional','TEMPLATE_functional.tif'));
mean_functional_2 = read_file(fullfile(analysis_folder,'2','MC functional','TEMPLATE_functional.tif'));
[h, w] = size(mean_functional_1);
mean_functional = uint16(2*[mean_functional_1 mean_functional_2]);

%resize figure to match aspect ratio of original image
aspect_ratio = w/h;
if aspect_ratio<1
    image_size(2) = round(image_size(2)*aspect_ratio);
elseif aspect_ratio>1
    image_size(1) = round(image_size(1)/aspect_ratio);
end

figure('Position',[100 100 fliplr(image_size)]);
imshow(imresize(mean_functional,image_size));
axis off
box off
a = gca;
a.Position = [0 0 1 1];
F = getframe(gcf);
[X, ~] = frame2im(F);
Xall = imresize(X,image_size);
close(gcf)


%% 2. unmixed yellow
color_im_1 = read_file(fullfile(analysis_folder,'1','celltypes_unmixed.tif'));
color_im_2 = read_file(fullfile(analysis_folder,'2','celltypes_unmixed.tif'));

figure('Position',[100 100 fliplr(image_size)]);
yellow_1 = color_im_1(:,:,2);
yellow_1 = 400*(double(yellow_1)/double(max(yellow_1,[],'all')));
yellow_2 = color_im_2(:,:,2);
yellow_2 = 400*(double(yellow_2)/double(max(yellow_2,[],'all')));
color_im = uint8(0.5*[yellow_1 yellow_2]);
color_im(:,:,2) = color_im;
color_im(:,:,3) = uint8(0);
imshow(imresize(color_im,image_size));
axis off
box off
a = gca;
a.Position = [0 0 1 1];
F = getframe(gcf);
[X, ~] = frame2im(F);
Xall(:,:,:,2) = imresize(X,image_size);
close(gcf)


%% 3. unmixed red
color_im_1 = read_file(fullfile(analysis_folder,'1','celltypes_unmixed.tif'));
color_im_2 = read_file(fullfile(analysis_folder,'2','celltypes_unmixed.tif'));

figure('Position',[100 100 fliplr(image_size)]);
red_1 = color_im_1(:,:,3);
red_1 = 400*(double(red_1)/double(max(red_1,[],'all')));
red_2 = color_im_2(:,:,3);
red_2 = 400*(double(red_2)/double(max(red_2,[],'all')));
color_im = uint8([red_1 red_2]);
color_im(:,:,[2 3]) = uint8(0);
imshow(imresize(color_im,image_size));
axis off
box off
a = gca;
a.Position = [0 0 1 1];
F = getframe(gcf);
[X, ~] = frame2im(F);
Xall(:,:,:,3) = imresize(X,image_size);
close(gcf)


%% 4. ROIs
figure('Position',[100 100 fliplr(image_size)]);

roiImg = zeros(h,w,3,2);
for t = 1:2
    if t==1
        stat = stat_1;
        cells50 = find(iscell_1(:,1));
        num_rois = length(cells50);
        colors = colors_1;
    else
        stat = stat_2;
        cells50 = find(iscell_2(:,1));
        num_rois = length(cells50);
        colors = colors_2;
    end
    for r = 1:num_rois
        roi = cells50(r);
        for c = 1:3
            rowsub = stat{roi}.ypix+1;
            colsub = stat{roi}.xpix+1;
            colorsub = c*ones(size(rowsub));
            tsub = t*ones(size(rowsub));
            brightnessvec = stat{roi}.lam;
            brightnessvec = brightnessvec/max(brightnessvec);
            roiImg(sub2ind([h w 3 2],rowsub,colsub,colorsub,tsub)) = colors(r,c)*brightnessvec;
        end
    end
end
roiImg = uint8(255*[roiImg(:,:,:,1) roiImg(:,:,:,2)]);
image(imresize(roiImg,image_size));
axis off
box off
a = gca;
a.Position = [0 0 1 1];
F = getframe(gcf);
[X, ~] = frame2im(F);
Xall(:,:,:,4) = imresize(X,image_size);
close(gcf)


%% 5. ROI #s
figure('Position',[100 100 fliplr(image_size)]);

roiImg = uint8(zeros(h,2*w,3));
scalefactor = image_size(1)/256;
image(imresize(roiImg,image_size));
for t = 1:2
    if t==1
        stat = stat_1;
        cells50 = find(iscell_1(:,1));
        num_rois = length(cells50);
        colors = colors_1;
        xadd = 0;
    else
        stat = stat_2;
        cells50 = find(iscell_2(:,1));
        num_rois = length(cells50);
        colors = colors_2;
        xadd = scalefactor*256;
    end
    for r = 1:num_rois
        roi = cells50(r);
        x = xadd + scalefactor*stat{roi}.med(2);
        y = scalefactor*stat{roi}.med(1);
        text(double(x),double(y),num2str(roi),'Color',colors(r,:));
    end
end
axis off
box off
a = gca;
a.Position = [0 0 1 1];
F = getframe(gcf);
[X, ~] = frame2im(F);
Xall(:,:,:,5) = imresize(X,image_size);
close(gcf)


%% save the stack
savename = fullfile(analysis_folder,'CellLinkingSideBySide.tif');
options.color = true;
options.compress = 'no';
options.message = false;
options.append = false;
options.overwrite = true;
options.big = false;
saveastiff(Xall,savename,options)

savename = fullfile(analysis_folder,'CellLinkingOverlayed.tif');
numFrames = size(Xall,4);
w = floor(image_size(2)/2);
XallS = Xall(:,1:w,:,:);
XallS(:,:,:,(numFrames+1):(2*numFrames)) = Xall(:,w+1:2*w,:,:);
saveastiff(XallS,savename,options)


%% run the cell linker gui with the stack name
%TODO