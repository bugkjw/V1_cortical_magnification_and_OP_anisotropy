function [Result] = v20181017_spacing(OR_map, rr_sig)
%% 8. Acquiring mean column spacing: NN and projection
% Autocorrelation method
autocorr = xcorr2(exp(2i*OR_map),exp(2i*OR_map));

uniform = real(autocorr);
contrast = (uniform)/mean(mean(uniform)).^8;
p = FastPeakFind(contrast);
xpos = p(1:2:end)'; ypos = p(2:2:end)';
pos = [xpos;ypos];
NN6 = [];
for ii = 1:size(pos,2)
    nndist = [];
    for jj = 1:size(pos,2)
        if ii~=jj
            nndist = [nndist, sum((pos(:,ii)-pos(:,jj)).^2).^0.5];
        end
    end
    temp = sort(nndist);
    NN6 = [NN6 mean(temp(1:min(6,length(temp))))];
end
s_mean_NN = mean(NN6);

% figure;
% ax1 = subplot(1,2,1); hold on;
% imagesc(OR_map);
% title(['Simulated orientation map, Gaussian std = ' num2str(rr_sig)]);
% colormap(ax1,hsv); caxis([-pi/2 pi/2]);
% axis xy image
% c = colorbar; c.Label.String = 'Preferred orientation (radian)';
% 
% subplot(1,2,2); hold on;
% imagesc(real(autocorr));
% axis xy image;
% plot(p(1:2:end),p(2:2:end),'r+');
% title(['Map autocorrelation, Gaussian std = ' num2str(rr_sig)]);
% x0 = 100; y0 = 100; width = 1000; height = 400;
% set(gcf,'units','points','position',[x0,y0,width,height]);

figure;
findpeaks(sum(contrast,1));
[~,locs] = findpeaks(sum(contrast,1));
s_mean_proj = mean(locs(2:end)-locs(1:end-1));
close;

Result = [s_mean_NN s_mean_proj];

end