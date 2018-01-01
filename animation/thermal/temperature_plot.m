function temperature_plot(xyzs_file, radius_file, temperature_file, filename)
h = figure;
set(h, 'Visible', 'off');

xyzs = dlmread(xyzs_file);
r_p = dlmread(radius_file);
ts = dlmread(temperature_file);

partical_x_pos = xyzs(:,1);
partical_y_pos = xyzs(:,2);
partical_z_pos = xyzs(:,3);
temperatures = ts;
colors = vals2colormap(temperatures, 'jet');

len = length(partical_z_pos);

for i = 1:len
    [x, y, z] = ellipsoid(partical_x_pos(i),partical_y_pos(i),...
        partical_z_pos(i),r_p(i),r_p(i),r_p(i),10);

    surf(x,y,z,'FaceColor',colors(i, :),...
        'FaceLighting','gouraud','EdgeColor','none'); hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'xtick',[]);
set(gca, 'ytick', []);
set(gca, 'ztick', []);
axis equal
view([1 1 1]);
camlight
colormap('jet')
colorbar('SouthOutside')
caxis([0, 4000])
if nargin == 4
    saveas(gcf, filename);
end

end

