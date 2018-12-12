function [Result] = v20181017_OP_measure(V1_N_pos, V1_layers, rsize, r_resol, Mx, POS)
V1_sample_pos = cell2mat(POS(3));
V1_pos = cell2mat(POS(4));

%% 5. Cortical map measurement
Pref_ori = zeros(V1_N_pos,1);
yf = size(V1_layers,1); xf = size(V1_layers,2);
[xx,yy] = meshgrid((1:xf)-round(xf/2),(1:xf)-round(xf/2));
[~,rr] = cart2pol(xx,yy);
map_angle = angle(xx+1i*yy);
map_angle_complex = exp(2i*map_angle);

for ii = 1:V1_N_pos
    temp_pos = round(1+[rsize;rsize*Mx]/r_resol+V1_sample_pos(:,ii)/r_resol);

    winSize = xf;
    %temp_x = temp_pos(1);
    temp_y = temp_pos(2);
    if temp_y<=winSize/2
        Local_pool = V1_layers(1:xf,1:xf,ii);
    elseif temp_y>=yf-winSize/2
        Local_pool = V1_layers(yf-winSize+1:yf,1:xf,ii);
    elseif mod(winSize,2)==0
        Local_pool = V1_layers...
            (temp_y-winSize/2:temp_y+winSize/2-1,1:xf,ii);
    else %if mod(wiSize,2)==1
        Local_pool = V1_layers...
            (temp_y-floor(winSize/2):temp_y+floor(winSize/2),1:xf,ii);
    end
    
    Y = fftshift(fft2(Local_pool));
    
    temp = abs(Y).*map_angle_complex.*rr;
    
    mu = sum(temp(:))/sum(abs(Y(:)));
    pref_angle_temp = angle(mu)/2;
    disp(['Measuring... ' num2str(round(ii/V1_N_pos*100)) '% ' num2str(pref_angle_temp) 'rad']);
    Pref_ori(ii) = pref_angle_temp+pi/2;
end

Pref_ori(Pref_ori>pi/2) = Pref_ori(Pref_ori>pi/2)-pi;

% figure; subplot(1,2,1);
% title('Sampled OPs');
% scatter(V1_sample_pos(1,:),V1_sample_pos(2,:),10,Pref_ori,'filled');
% axis xy image off;
% colormap(hsv); colorbar;
% 
% subplot(1,2,2);
% scatter(V1_pos(1,:), V1_pos(2,:), 10, Pref_ori, 'filled');
% title('OP map, unfiltered'); colorbar; colormap(hsv); axis xy image;

Result = Pref_ori;

end