%% CLEAR
close all;
clear all;

%% OP Map Data Import & Preprocessing
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

% Gaussian Noise Reduction
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

% figure;
% subplot(1,2,1); image(DATA_MAP); axis image; title("RGB data image");
% subplot(1,2,2); imagesc(OP_MAP_F); hold on;
% plot(Pwl_p(1,:),Pwl_p(2,:),'wo','MarkerSize',4, 'LineWidth',1);
% plot(Pwl_n(1,:),Pwl_n(2,:),'ko','MarkerSize',4, 'LineWidth',1);
% title('Pinwheel analysis'); colormap(hsv); caxis([0 pi]); axis image;
% c = colorbar; c.Label.String = 'Preferred orientation (radian)';

%% Local Region Analysis
LOCAL_SIZE = 40;
C_spacing = ones(size(OP_MAP,1), size(OP_MAP,2))*nan;
Pwl_num = ones(size(OP_MAP,1), size(OP_MAP,2))*nan;
Pwl_density = ones(size(OP_MAP,1), size(OP_MAP,2))*nan;
for ii = 1:size(OP_MAP,1)
    for jj = 1:size(OP_MAP,2)
        l1 = max(ii-LOCAL_SIZE,1);
        u1 = min(ii+LOCAL_SIZE-1,size(OP_MAP,1));
        l2 = max(jj-LOCAL_SIZE,1);
        u2 = min(jj+LOCAL_SIZE-1,size(OP_MAP,2));
        % Assume that measured V1 area is a convex polygon
        if ~sum(isnan([OP_MAP(ii,jj), OP_MAP(l1,l2), OP_MAP(l1,u2), OP_MAP(u1,l2), OP_MAP(u1,u2)]))
            LOCAL_MAP = OP_MAP(l1:u1,l2:u2);
            % Local column spacing
            C_spacing(ii,jj) = v20181017_Column_Spacing(LOCAL_MAP);
            if isnan(C_spacing(ii,jj)) || C_spacing(ii,jj) == 0
                C_spacing(ii,jj) = C_spacing(ii-1,jj-1);
                % figure; imagesc(LOCAL_MAP); axis image; colormap(hsv); colorbar;
                % figure; imagesc(real(xcorr2(exp(2i*LOCAL_MAP),exp(2i*LOCAL_MAP)))); error('Column spacing not found');
            end
            % Pinwheel number in local region
            Pwl_num(ii,jj) = sum((Pwl_p(1,:)>=l2) & (Pwl_p(1,:)<=u2) & (Pwl_p(2,:)>=l1) & (Pwl_p(2,:)<=u1)) ...
                + sum((Pwl_n(1,:)>=l2) & (Pwl_n(1,:)<=u2) & (Pwl_n(2,:)>=l1) & (Pwl_n(2,:)<=u1));
            % Local pinwheel density per hypercolumn
            Pwl_density(ii,jj) = Pwl_num(ii,jj)/(size(LOCAL_MAP,1)*size(LOCAL_MAP,2))*C_spacing(ii,jj)^2;
        end
        fprintf("%.3g percent complete...\n",((ii-1)*size(OP_MAP,2)+jj)/(size(OP_MAP,1)*size(OP_MAP,2))*100);
    end
end

figure('position', [10, 10, 1400, 900]);
c = subplot(2,2,1); imagesc(OP_MAP_F); hold on;
plot(Pwl_p(1,:),Pwl_p(2,:),'wo','MarkerSize',4, 'LineWidth',1);
plot(Pwl_n(1,:),Pwl_n(2,:),'ko','MarkerSize',4, 'LineWidth',1);
title('Pinwheel analysis'); colormap(c,hsv); caxis([0 pi]); axis image;
c = colorbar; c.Label.String = 'Preferred orientation (radian)';
d = subplot(2,2,2); imagesc(Pwl_num);
title("Local pinwheel counts"); axis image; colorbar; colormap(d,[0 0 0;hot]);
e = subplot(2,2,3); imagesc(C_spacing);
title("Local column spacing"); axis image; colorbar; colormap(e,[0 0 0;hot]); caxis([0 100]);
f = subplot(2,2,4); imagesc(Pwl_density);
title("Local pinwheel density/hypercolumn"); axis image; colorbar; colormap(f,[0 0 0;hot]); caxis([0 8]);

disp(nanmean(nanmean(Pwl_density)));

%% Thresholding
thr = pi; THR = Pwl_density;
U = double(Pwl_density>=thr); U(U==0) = nan;
U = Pwl_density.*(U); U_mean = nanmean(nanmean(U));
D = double(Pwl_density<thr); D(D==0) = nan;
D = Pwl_density.*(D); D_mean = nanmean(nanmean(D));
THR = zeros(size(OP_MAP,1), size(OP_MAP,2));
THR = THR+~isnan(U).*U_mean;
THR = THR+~isnan(D).*D_mean;

U_size = sum(sum(~isnan(U)));
D_size = sum(sum(~isnan(D)));
UD_ratio = [U_size D_size]/(U_size+D_size);
UD_ratio_pi = [(pi-D_mean)/(U_mean-D_mean) (U_mean-pi)/(U_mean-D_mean)];
Mean = nanmean(nanmean(Pwl_density));

disp(UD_ratio);
disp(UD_ratio_pi);

figure; imagesc(THR); axis image;
colormap(jet); colorbar;
title(['Upper density: ',num2str(U_mean),', Lower density: ',num2str(D_mean),', Mean: ',num2str(Mean)]);