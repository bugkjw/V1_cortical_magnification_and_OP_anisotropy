function [Result] = v20181017_filt_Gaussian(r_resol, c_resol, csize, rr_sig, V1_N_pos, POS, Pref_ori)
V1_pos = cell2mat(POS(4));

%% 6. Map filtering (Gaussian)
[xx,yy] = meshgrid(-csize:csize,-csize:csize);
xy = xx+1i*yy;
OR_Cmap = zeros(length(-csize:csize),length(-csize:csize));

for ii = 1:V1_N_pos
    disp(['Filtering with Gaussian... ' num2str(round(ii/V1_N_pos*100)) '%']);
    pos_temp = V1_pos(1,ii)+1i*V1_pos(2,ii); %% pos_temp = x_pos + 1i*y_pos
    rr_temp = abs(xy-pos_temp)*r_resol/c_resol;
    gaus_temp = 1/rr_sig^2/(2*pi)*exp(-(rr_temp/rr_sig).^2/2);
    
    ori_temp = exp(2i*Pref_ori(ii));
    OR_Cmap = OR_Cmap + gaus_temp.*ori_temp;
end
OR_map = angle(OR_Cmap)/2;
OR_map(OR_map>pi/2) = OR_map(OR_map>pi/2)-pi;

Result = OR_map;

end