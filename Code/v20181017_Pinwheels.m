function [Result] = v20181017_Pinwheels(csize, OR_map)
%% 7. Pinwheel detection
% Basic idea: point of intersection between real and imag zero-magnitude contour of complex orientation map
figure; hold on;
[xximg,yyimg] = meshgrid(1:2*csize+1,1:2*csize+1);

imagesc(OR_map); title('Simulated orientation map');
colormap(hsv); caxis([-pi/2 pi/2])
c = colorbar; c.Label.String = 'Preferred orientation (radian)';

OR_Cmap = exp(2i*OR_map);
OR_Rmap = real(OR_Cmap);
OR_Imap = imag(OR_Cmap);
[C1,~] = contour(xximg,yyimg,OR_Rmap);
[C2,~] = contour(xximg,yyimg,OR_Imap);
axis xy image off;

[x1,y1,z1] = C2xyz(C1);
[x2,y2,z2] = C2xyz(C2);

x1(z1~=0)=[]; y1(z1~=0)=[]; z1(z1~=0)=[];
x2(z2~=0)=[]; y2(z2~=0)=[]; z2(z2~=0)=[];
Pwl = []; Pwl_p = []; Pwl_n = [];
for ii = 1:length(z1)
    x1temp = cell2mat(x1(ii));
    y1temp = cell2mat(y1(ii));
    for jj = 1:length(z2)
        x2temp = cell2mat(x2(jj));
        y2temp = cell2mat(y2(jj));
        Ptemp = InterX([x1temp;y1temp],[x2temp;y2temp]);
        Pwl = [Pwl Ptemp];
    end
end
for Ptemp = Pwl
    img_temp = OR_map...
        (max(1,round(Ptemp(2))-6):min(size(OR_map,1),round(Ptemp(2))+6),...
        max(1,round(Ptemp(1))-6):min(size(OR_map,2),round(Ptemp(1))+6));
    cnt = 0;
    for kk = 1:size(img_temp,1)-1
        if img_temp(kk+1,1)-img_temp(kk,1)>=0
            cnt = cnt+1;
        else
            cnt = cnt-1;
        end
    end
    for kk = 1:size(img_temp,2)-1
        if img_temp(end,kk+1)-img_temp(end,kk)>=0
            cnt = cnt+1;
        else
            cnt = cnt-1;
        end
    end
    for kk = size(img_temp,1):1
        if img_temp(kk-1,end)-img_temp(kk,end)>=0
            cnt = cnt+1;
        else
            cnt = cnt-1;
        end
    end
    for kk = size(img_temp,2):1
        if img_temp(1,kk-1)-img_temp(1,kk)>=0
            cnt = cnt+1;
        else
            cnt = cnt-1;
        end
    end
    %             figure; hold on;
    %             imagesc(img_temp); axis xy image off; colormap(hsv); colorbar;
    if cnt>0
        Pwl_p = [Pwl_p Ptemp];
        %                 plot(size(img_temp,2)/2,size(img_temp,1)/2,'ko');
    else
        Pwl_n = [Pwl_n Ptemp];
        %                 plot(size(img_temp,2)/2,size(img_temp,1)/2,'wo');
    end
end

close;
% figure; hold on;
% imagesc(OR_map);
% plot(Pwl_p(1,:),Pwl_p(2,:),'wo','LineWidth',2);
% plot(Pwl_n(1,:),Pwl_n(2,:),'ko','LineWidth',2);
% title('Simulated orientation map with pinwheels');
% colormap(hsv); caxis([-pi/2 pi/2])
% axis xy image off
% c = colorbar; c.Label.String = 'Preferred orientation (radian)';

Result = {Pwl_p Pwl_n};

end