%% CLEAR
close all;
clear all;

%% 1. OP Map Data Import & Preprocessing
% jpg: ? X ? X 3 uint8 format
DATA_MAP = imread('F1_cut.jpg');

% White Space Filling (using 0)
DATA_MAP_reshape = reshape(DATA_MAP, [size(DATA_MAP,1)*size(DATA_MAP,2) 3]);
thr = 190;
white_mask = (DATA_MAP_reshape(:,1)>thr) & (DATA_MAP_reshape(:,2)>thr) & (DATA_MAP_reshape(:,3)>thr);
DATA_MAP_reshape(white_mask,1) = 0;
DATA_MAP_reshape(white_mask,2) = 0;
DATA_MAP_reshape(white_mask,3) = 0;

DATA_MAP = reshape(DATA_MAP_reshape, [size(DATA_MAP,1), size(DATA_MAP,2) 3]);
% figure;
% subplot(1,2,1); image(DATA_MAP); title("Raw image"); axis image;
% subplot(1,2,2); image(DATA_MAP); title("Blank space removal"); axis image;

% RGB - OPMAP Conversion
res = 500;
OP_MAP = v20181008_RGB_OP_Mapping(res);

% Gaussian Denoising
OP_MAP_F = angle(imgaussfilt(real(exp(2i*OP_MAP)),2.5)+1i*imgaussfilt(imag(exp(2i*OP_MAP)),2.5))/2;
NEW_OP_MAP = OP_MAP_F;
OP_MAP_F(OP_MAP_F<0) = OP_MAP_F(OP_MAP_F<0)+pi;

% figure;
% subplot(2,2,1); image(DATA_MAP); axis image; title("RGB data image");
% subplot(2,2,2); scatter3(RGB(1,:),RGB(2,:),RGB(3,:),10,ORI,'filled'); colormap(hsv); colorbar; title("RGB curve, mapped onto OP");
% xlabel("R"); ylabel("G"); zlabel("B");
% subplot(2,2,3); imagesc(OP_MAP); axis image; colormap(hsv); colorbar; title("Orientation image");
% subplot(2,2,4); imagesc(OP_MAP_F); axis image; colormap(hsv); colorbar; title("Filtered OP image");

OP_MAP = NEW_OP_MAP;

% Pinwheel Identification
Result = v20181008_Pinwheel_Identification(OP_MAP);
Pwl_p = cell2mat(Result(1));
Pwl_n = cell2mat(Result(2));

%% 2. Local Region Selection
figure('position', [100, 100, 1400, 500]); c = subplot(1,2,1);
imagesc(OP_MAP_F); hold on;
plot(Pwl_p(1,:),Pwl_p(2,:),'wo','MarkerSize',4, 'LineWidth',1);
plot(Pwl_n(1,:),Pwl_n(2,:),'ko','MarkerSize',4, 'LineWidth',1);
title('Pinwheel analysis'); colormap(c,hsv); caxis([0 pi]); axis image;
c = colorbar; c.Label.String = 'Preferred orientation (radian)';

LOCAL_SIZE = 50;
sy = 10:LOCAL_SIZE:size(OP_MAP,1); sx = 10:LOCAL_SIZE:size(OP_MAP,2);
valid = zeros([length(sy) length(sx)]);
for ii = 1:length(sy)
    for jj = 1:length(sx)
        pos1 = sx(jj); pos2 = sy(ii);
        if ~isnan(OP_MAP(pos2,pos1))
            plot(pos1, pos2, "kx",'MarkerSize',10, 'LineWidth',4);
            valid(ii,jj) = 1;
        end
    end
end

%% 3. Local Column Spacing and Pinwheel Density Calculation
% Assume that measured V1 area is a convex polygon
Pwl_density = ones([length(sy) length(sx)])*-1;
area_cnt = 0;
for ii = 1:length(sy)-1
    for jj = 1:length(sx)-1
        if valid(ii,jj) && valid(ii+1,jj) && valid(ii,jj+1) && valid(ii+1,jj+1)
            area_cnt = area_cnt+1;
            local_map = OP_MAP(sy(ii):sy(ii+1),sx(jj):sx(jj+1));
            Y = xcorr2(exp(2i*local_map));
            Y = abs(Y)/max(max(abs(Y)))^7; Y = abs(Y)/max(max(abs(Y)));
            
            Pwl_cnt =...
                sum(Pwl_p(1,:) >= sx(jj) & Pwl_p(2,:) >= sy(ii) & Pwl_p(1,:) < sx(jj+1) & Pwl_p(2,:) < sy(ii+1))+...
                sum(Pwl_n(1,:) >= sx(jj) & Pwl_n(2,:) >= sy(ii) & Pwl_n(1,:) < sx(jj+1) & Pwl_n(2,:) < sy(ii+1));
            
            PSEUDO_COLUMN_SPACING = LOCAL_SIZE;
            Pwl_density(ii,jj) = Pwl_cnt/LOCAL_SIZE^2*PSEUDO_COLUMN_SPACING^2;
            
%             figure('position', [0, 0, 1000, 500]);
%             c = subplot(1,2,1); imagesc(local_map); axis xy image; colormap(c,hsv); colorbar;
%             subplot(1,2,2); imagesc(Y); axis xy image; colorbar;
        end
    end
end

d = subplot(1,2,2); imagesc(Pwl_density);
title("Local pinwheel counts"); axis image; colorbar; colormap(d,[0 0 0;hot]);

% OP_MAP_NAN0 = OP_MAP(sy(2):sy(end-1),sx(4):sx(end-3));
% OP_MAP_NAN0(isnan(OP_MAP_NAN0)) = 0;
% Y = abs(xcorr2(exp(2i*(OP_MAP_NAN0))));
% Y = Y/max(max(Y));
% figure;
% imagesc(Y); axis image; colorbar;

%% 4. Pinwheel Density vs. Region
