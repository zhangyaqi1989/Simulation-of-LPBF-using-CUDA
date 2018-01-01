function particle_plot(xyzs_file, radius_file, xlims, ylims, zlims, filename)
h = figure;
set(h, 'Visible', 'off');

xyzs = dlmread(xyzs_file);
r_p = dlmread(radius_file);

partical_x_pos = xyzs(:,1);
partical_y_pos = xyzs(:,2);
partical_z_pos = xyzs(:,3);

len = length(partical_z_pos);

xlim(xlims);
ylim(ylims);
zlim(zlims);

xticks([0:2:10]*1e-4);

for i = 1:len
    [x, y, z] = ellipsoid(partical_x_pos(i),partical_y_pos(i),...
        partical_z_pos(i),r_p(i),r_p(i),r_p(i),10);
    surf( x,y,z,'FaceColor',[0 0 1], 'FaceLighting','gouraud',...
        'EdgeColor','none'); hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'xtick',[]);
set(gca, 'ytick', []);
set(gca, 'ztick', []);
axis equal
view([1 1 1])
camlight
if nargin == 6
    saveas(gcf, filename);
end
end
