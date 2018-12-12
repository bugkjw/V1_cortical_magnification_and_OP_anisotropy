function [Result] = v20181017_sampling(rsize, r_resol, Mx, csize, c_resol)
%% 2. Cortical Sampling
% Sample points in retinal space
r_sample_x = -rsize:r_resol:rsize;
r_sample_y = -rsize*Mx:r_resol*Mx:rsize*Mx;

% Corresponding points in cortical space
c_sample_x = -csize:c_resol:csize;
c_sample_y = -csize:c_resol:csize;

V1_sample_pos = [];
V1_pos = [];
for ii = 1:min(length(r_sample_y),length(c_sample_y))
    for jj = 1:min(length(r_sample_y),length(c_sample_y))
        V1_sample_pos = [V1_sample_pos [r_sample_x(ii);r_sample_y(jj)]];
        V1_pos = [V1_pos [c_sample_x(ii); c_sample_y(jj)]];
    end
end

figure; subplot(1,2,1);
scatter(V1_sample_pos(1,:),V1_sample_pos(2,:),10);
axis image; title('OP sample points');

subplot(1,2,2);
scatter(V1_pos(1,:),V1_pos(2,:),10);
axis image; title('V1 cell pos, cortical space');

Result = {V1_sample_pos, V1_pos};

end