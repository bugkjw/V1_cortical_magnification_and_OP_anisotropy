%% Simulating orientation maps as moire interference patterns between RGC mosaics
clear all;
close all;

% Retinal space definition
rsize = 700; r_resol = 10;
rspace = -rsize:r_resol:rsize;
r_msize = length(rspace);
% Cortical space definition
csize = 100; c_resol = 2;
cspace = -csize:c_resol:csize;
c_msize = length(cspace);

% RGC mosaic
alpha = 0.1;
d = 30; % Grid spacing for OFF mosaic
noise_factor = 0.12; %Noise std in respect to d, 0 to 0.12
phi = 10*pi/180; theta = 5*pi/180;
Result = v20181017_RGC_mosaic(alpha, d, noise_factor,phi,theta, rsize*3);
OFF_pos = cell2mat(Result(1));
ON_pos = cell2mat(Result(2));

% Map smoothing std
rr_sig = 30;

N = 100;
Pwl_density = zeros(3,10,N);

for trial = 1:N
    for ee = 1:11
        My = 1; % mmc/deg
        aniso_index = 1+(ee-1)/10; % Anisotropy index = My/Mx, should be >1 !!
        Mx = My*aniso_index; % mmc/deg
        
        OFF_pos_temp = []; ON_pos_temp = [];
        for ii = 1:size(OFF_pos,2)
            pos_temp = OFF_pos(:,ii);
            if pos_temp(1)>-rsize && pos_temp(1)<rsize && pos_temp(2)>-rsize*Mx && pos_temp(2)<rsize*Mx
                OFF_pos_temp = [OFF_pos_temp pos_temp];
            end
        end
        for ii = 1:size(ON_pos,2)
            pos_temp = ON_pos(:,ii);
            if pos_temp(1)>-rsize && pos_temp(1)<rsize && pos_temp(2)>-rsize*Mx && pos_temp(2)<rsize*Mx
                ON_pos_temp = [ON_pos_temp pos_temp];
            end
        end
        
        Result = v20181017_sampling(rsize, r_resol, Mx, csize, c_resol);
        V1_sample_pos = cell2mat(Result(1)); V1_pos = cell2mat(Result(2));
        
        V1_N_pos = size(V1_pos,2);
        POS = {OFF_pos_temp, ON_pos_temp, V1_sample_pos, V1_pos};
        
        V1_layers = v20181017_V1_RF(rsize, r_resol, Mx, d, POS, V1_N_pos);
        
        Pref_ori = v20181017_OP_measure(V1_N_pos, V1_layers, rsize, r_resol, Mx, POS);
        
        OR_map = v20181017_filt_Gaussian(r_resol, c_resol, csize, rr_sig, V1_N_pos, POS, Pref_ori);
        
        Result = v20181017_Pinwheels(csize, OR_map);
        Pwl_p = cell2mat(Result(1)); Pwl_n = cell2mat(Result(2));
        
        s_mean = v20181017_spacing(OR_map, rr_sig);
        s_mean_NN = s_mean(1); s_mean_proj = s_mean(2);
        
        %% Mean pinwheel density
        if rr_sig>=40 && aniso_index>=1.5
            s_mean = s_mean_proj;
        else
            s_mean = s_mean_NN;
        end
        Pwl_density(1,ee,trial) = size(Pwl_p,2)/(2*csize+1)^2*s_mean^2;
        Pwl_density(2,ee,trial) = size(Pwl_n,2)/(2*csize+1)^2*s_mean^2;
        Pwl_density(3,ee,trial) = (size(Pwl_p,2)+size(Pwl_n,2))/(2*csize+1)^2*s_mean^2;
        disp(['Pinwheel(+) density: ' num2str(Pwl_density(1,ee,trial)) ' per hypercolumn']);
        disp(['Pinwheel(-) density: ' num2str(Pwl_density(2,ee,trial)) ' per hypercolumn']);
        disp(['Pinwheel density: ' num2str(Pwl_density(3,ee,trial)) ' per hypercolumn']);
        
    end
end

figure; hold on;
xx = 1+(0:10)/10;
plot(xx, Pwl_density(1,:,1),'r'); plot(xx, Pwl_density(2,:,1),'b'); plot(xx, Pwl_density(3,:,1),'k');
xlabel('Aniotropy index'); ylabel('Pinwheel density');
legend('+','-','Total'); axis tight; ylim([0 6]);