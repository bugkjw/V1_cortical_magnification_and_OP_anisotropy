function [Result] = v20181017_RGC_mosaic(alpha, d, noise_factor,phi,theta, rsize)
%% 1. RGC mosaic generator
OFF_pos = []; ON_pos = []; % RF centers

rotMat_whole = [cos(phi) sin(phi); -sin(phi) cos(phi)];
rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];

hexMat = 1/2*[1 1;sqrt(3) -sqrt(3)];
displacement = 0;
maxij = 200;
for i = -maxij:maxij
    for j = -maxij:maxij
        Lij = rotMat_whole*hexMat*[i;j]; % Hexagonal grid
        nij = randn(2,1)*noise_factor*d;
        c_OFF = d*Lij+nij; c_ON = (1+alpha)*d*rotMat*Lij+nij+displacement;
        if c_OFF(1)>-rsize && c_OFF(1)<rsize && c_OFF(2)>-rsize && c_OFF(2)<rsize
            OFF_pos = [OFF_pos c_OFF];
        end
        if c_ON(1)>-rsize && c_ON(1)<rsize && c_ON(2)>-rsize && c_ON(2)<rsize
            ON_pos = [ON_pos c_ON];
        end
    end
end

scaling_factor = (1+alpha)/sqrt(alpha^2+2*(1-cos(theta))*(1+alpha));
figure; suptitle('Total retinal space');
hold on; plot(OFF_pos(1,:),OFF_pos(2,:),'.r'); plot(ON_pos(1,:),ON_pos(2,:),'.b');
xlim([-rsize rsize]); ylim([-rsize rsize]); axis image;
title(['ON/OFF RGC mosaic, S = ', num2str(scaling_factor)]);
pbaspect([1 1 1]);

Result = {OFF_pos, ON_pos};

end