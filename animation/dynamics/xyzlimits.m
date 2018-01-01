function [xlims, ylims, zlims] =  xyzlimits()

steps = [0:100];
radius_file = 'radius.txt';
rs = dlmread(radius_file);
count = 1;
nparts = size(rs, 1);
xyzmins = zeros(nparts, 3);
xyzmaxs = zeros(nparts, 3);
for step = steps
    xyzs_file = ['xyzs', sprintf('%03d', step), '.txt'];
    disp(xyzs_file);
    xyzs = dlmread(xyzs_file);
    xyzmins(count, :) = min(xyzs - rs);
    xyzmaxs(count, :) = max(xyzs + rs);
    count = count + 1;
end
maxs = num2cell(max(xyzmaxs));
mins = num2cell(min(xyzmins));

[x_max, y_max, z_max] = maxs{:};
[x_min, y_min, z_min] = mins{:};
xlims = [x_min, x_max];
ylims = [y_min, y_max];
zlims = [z_min, z_max];
end
