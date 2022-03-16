function[piezo] = new_chris_piezo(sz, v_scale)

addpath 'C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Usefull stuff'
% sz = 128;
% v_scale = 3.4*1e-7;

M = 4; % rads
std_nor = v_scale;

% stack = CStack(num2cell(norm_vels));

xs = linspace(-1,1,sz);
ys = xs;
R = 0.7;
[xx,yy] = meshgrid(xs,ys);

theta = atan2(yy,xx);
rr = sqrt(xx.^2 + yy.^2);

piezo = rr < R;
piezo = cast(piezo, 'double');

intensity = linspace(-1, 1, 20); 
radiusN = [0, 0.2625, 0.2625 + 0.21875, 0.2625 + 0.21875 + 0.21875];
radiusN = [0, 0.275758, 0.275758 + 0.212121, 0.275758 + 0.212121 + 0.212121];

angs = 17;
    thethaN = linspace(-pi, pi, angs);
    thethaN_inner = linspace(-pi, pi, 9);
    dif = diff(thethaN_inner);
    dif = dif(1) / 2;
    
    outerTan = [];
    for k = 1 : length(thethaN_inner) - 1
        if mod(k, 2) == 1
            outerTan = [outerTan;  thethaN_inner(k)];
            outerTan = [outerTan;  outerTan(end) + dif];
        else 
            outerTan = [outerTan;  thethaN_inner(k)];
            outerTan = [outerTan;  outerTan(end) + dif];

        end
    end
    outerTan = [outerTan;  outerTan(end) + dif];
    

for i = 1:M-1
    
    rang = linspace(1, angs - 1, angs/2);

    if i == 1
         j = 1;
         while j <= length(thethaN_inner) - 1
             curr_slice = theta > thethaN_inner(j) & theta < thethaN_inner(j + 1) & ...
                 rr > radiusN(i) & rr < radiusN(i + 1);
                 curr_slice = curr_slice * intensity(randi([1 20]));
                 piezo = piezo + curr_slice;
                 j = j + 1;
         end
    else
        
        for j = 1 : angs - 1
%             if ~stack.isempty()
                curr_slice = theta < outerTan(j+1) & theta > outerTan(j) & ...
                                            rr > radiusN(i) & rr < radiusN(i + 1);
                curr_slice = curr_slice * intensity(randi([1 20]));
                piezo = piezo + curr_slice;
%             end
        end 
    end
end


% figure, image(circle,'CDatamapping','scaled'),set(gca,'FontSize',20);
piezo = piezo .* std_nor;
% pcolor(xs, ys, piezo);
title('Piezo Spacial varying Hz rage of [-0.5 0,5]');

% lambda = 532*1e-9;
% vels_range = lambda .* Hz_rage;
% 
% num_of_piezos = 40; % 8 + 16 + 16
% vels = randperm(50, num_of_piezos);
% norm_vels = (max(vels_range) - min(vels_range)) * ((vels - min(vels))...
%     / (max(vels) - min(vels))) + min(vels_range);

% std_nor = std(norm_vels);

end