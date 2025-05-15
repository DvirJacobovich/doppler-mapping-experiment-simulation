% Dvir Jacobovich 2024 - Pr John Howell lab Hebrew University Of Jerusalem

function[circle] = chrs_pezo(sz)

% sz = 128;
M = 4; % rads

xs = linspace(-1,1,sz);
ys = xs;
R = 0.7;
[xx,yy] = meshgrid(xs,ys);

theta = atan2(yy,xx);
rr = sqrt(xx.^2 + yy.^2);

circle = rr < R;
circle = cast(circle, 'double');

% inner_rad = 0.3;

intensity = linspace(-1, 1, 20); 
radiusN = [0, 0.2625, 0.2625 + 0.21875, 0.2625 + 0.21875 + 0.21875];
radiusN = [0, 0.275758, 0.275758 + 0.212121, 0.275758 + 0.212121 + 0.212121];

for i = 1:M-1
    angs = 16;
    thethaN = linspace(-pi, pi, angs);
    %rang = linspace(1, angs, angs/2);
    rang = linspace(1, angs - 1, angs/2);

    if i == 1
%         for j = 1 : rang(~ismember(rang, 2:2:rang(end)))
        for j = rang(1:end)
            final_point = theta < thethaN(mod(j + 3, 15)) & theta > thethaN(mod(j + 1, 15)) & ...
                                        rr > radiusN(i) & rr < radiusN(i + 1);
            final_point = final_point * intensity(randi([1 20])) + 1.5;
            circle = circle + final_point;
        end 
        
    else
        
        for j = 1 : angs - 1
            final_point = theta < thethaN(j+1) & theta > thethaN(j) & ...
                                        rr > radiusN(i) & rr < radiusN(i + 1);
            final_point = final_point * intensity(randi([1 20]));
            circle = circle + final_point;
        end 
    end
end

pcolor(xs, ys, circle);

end
