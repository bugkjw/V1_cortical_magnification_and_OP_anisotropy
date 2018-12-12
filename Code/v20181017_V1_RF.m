function [Result] = v20181017_V1_RF(rsize, r_resol, Mx, d, POS, V1_N_pos)
OFF_pos = cell2mat(POS(1));
ON_pos = cell2mat(POS(2));
V1_sample_pos = cell2mat(POS(3));

%% 3. Computing LGN RF and weight
% Defined in retinal space.
x_rspace = -rsize:r_resol:rsize; x_rsize = length(x_rspace);
y_rspace = -rsize*Mx:r_resol:rsize*Mx; y_rsize = length(y_rspace);

% Weight assignment
pool_sig = 0.5*d;
dist_mag = [1 1]; % sqrt(2)*[Mx My]/((My^2+Mx^2)^0.5);
ON_weight = zeros(V1_N_pos,size(ON_pos,2));
OFF_weight = zeros(V1_N_pos,size(OFF_pos,2));

for ii = 1:V1_N_pos
    pos_temp = V1_sample_pos(:,ii);
    
    ON_dist = sum(dist_mag*((ON_pos-...
        [ones(1,size(ON_pos,2))*pos_temp(1);
        ones(1,size(ON_pos,2))*pos_temp(2)]).^2),1).^0.5; % 1 X size(ON_pos,2)
    ON_weight(ii,:) = exp(-ON_dist.^2/pool_sig^2);
    
    OFF_dist = sum(dist_mag*((OFF_pos-...
        [ones(1,size(OFF_pos,2))*pos_temp(1);
        ones(1,size(OFF_pos,2))*pos_temp(2)]).^2),1).^0.5; % 1 X size(OFF_pos,2)
    OFF_weight(ii,:) = exp(-OFF_dist.^2/pool_sig^2);
end

% RF in visual (retinal) space
sig_cen = 0.2*d; sig_sur = 1.0*d; Ncoef = 5;
ON_layers = zeros(y_rsize,x_rsize,size(ON_pos,2));
OFF_layers = zeros(y_rsize,x_rsize,size(OFF_pos,2));

for ii = 1:size(ON_pos,2)
    pos_temp = ON_pos(:,ii); xpos = pos_temp(1); ypos = pos_temp(2);
    [xx,yy] = meshgrid(x_rspace, y_rspace);
    RFcen = 1/(2*pi*sig_cen^2)*exp(-((xx-xpos).^2+(yy-ypos).^2)/(2*sig_cen^2));
    RFsur = -Ncoef/(2*pi*sig_sur^2)*exp(-((xx-xpos).^2+(yy-ypos).^2)/(2*sig_sur^2));
    layer = RFcen+RFsur; layer = layer./max(max(abs(layer)));
    ON_layers(:,:,ii) = layer;
end
for ii = 1:size(OFF_pos,2)
    pos_temp = OFF_pos(:,ii); xpos = pos_temp(1); ypos = pos_temp(2);
    [xx,yy] = meshgrid(x_rspace, y_rspace);
    RFcen = 1/(2*pi*sig_cen^2)*exp(-((xx-xpos).^2+(yy-ypos).^2)/(2*sig_cen^2));
    RFsur = -Ncoef/(2*pi*sig_sur^2)*exp(-((xx-xpos).^2+(yy-ypos).^2)/(2*sig_sur^2));
    layer = -(RFcen+RFsur); layer = layer./max(max(abs(layer)));
    OFF_layers(:,:,ii) = layer;
end

% %%
% figure; title('Pooling range of one particular V1 cell'); hold on;
% scatter(OFF_pos(1,:),OFF_pos(2,:),10,OFF_weight(2000,:),'filled')
% scatter(ON_pos(1,:),ON_pos(2,:),10,ON_weight(2000,:),'filled')
% axis image
% 
% figure; suptitle('RGC receptive field');
% tt = -3*d:0.01:3*d; [xx,yy] = meshgrid(-3*d:0.1:3*d,-3*d:0.1:3*d);
% subplot(1,2,1);
% RFcen = 1/(2*pi*sig_cen^2)*exp(-(tt.^2)/(2*sig_cen^2));
% RFsur = -Ncoef/(2*pi*sig_sur^2)*exp(-(tt.^2)/(2*sig_sur^2));
% layer = (RFcen+RFsur); plot(tt,layer);
% xlabel('Distance, in units of d_{OFF}');
% xticks(-3*d:d:3*d); xticklabels(-3:3); xlim([-3*d 3*d]);
% subplot(1,2,2);
% RFcen = 1/(2*pi*sig_cen^2)*exp(-(xx.^2+yy.^2)/(2*sig_cen^2));
% RFsur = -Ncoef/(2*pi*sig_sur^2)*exp(-(xx.^2+yy.^2)/(2*sig_sur^2));
% layer = (RFcen+RFsur); imagesc(layer); colormap(jet); colorbar;
% axis image off
% x0 = 100; y0 = 100; width = 1000; height = 400;
% set(gcf,'units','points','position',[x0,y0,width,height]);

%% 4. Computing V1 RF
ON_layers_2d = reshape(ON_layers,[y_rsize*x_rsize,size(ON_pos,2)]);
OFF_layers_2d = reshape(OFF_layers,[y_rsize*x_rsize,size(OFF_pos,2)]);

ON_weight_temp = ON_weight;
ON_weight_temp(ON_weight_temp<1e-6) = 0;
OFF_weight_temp = OFF_weight;
OFF_weight_temp(OFF_weight_temp<1e-6) = 0;

V1_layers_2d = ON_weight_temp*ON_layers_2d'+OFF_weight_temp*OFF_layers_2d';
V1_layers = reshape(V1_layers_2d',[y_rsize x_rsize V1_N_pos]);

Result = V1_layers;

end