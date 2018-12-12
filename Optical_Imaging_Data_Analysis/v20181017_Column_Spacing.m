function [Result] = v20181017_Column_Spacing(LOCAL_MAP)
%% 8. Acquiring mean column spacing: NN
% Autocorrelation method
autocorr = xcorr2(exp(2i*LOCAL_MAP),exp(2i*LOCAL_MAP));

uniform = real(autocorr);
contrast = (uniform)/mean(mean(uniform)).^6;
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
% title('OP Map Data');
% colormap(ax1,hsv); caxis([-pi/2 pi/2]);
% axis xy image
% c = colorbar; c.Label.String = 'Preferred orientation (radian)';
% 
% subplot(1,2,2); hold on;
% imagesc(real(autocorr));
% axis xy image;
% plot(p(1:2:end),p(2:2:end),'r+');
% title('Map Autocorrelation');
% x0 = 100; y0 = 100; width = 1000; height = 400;
% set(gcf,'units','points','position',[x0,y0,width,height]);

Result = s_mean_NN;

end