radius_file = 'radius.txt';
[xlims, ylims, zlims] = xyzlimits();

steps = [0:100];
count = 1;
for step = steps
    xyzs_file = ['xyzs', sprintf('%03d', step), '.txt'];
    disp(xyzs_file);
    image_file = sprintf('image%03d.png', count);
    particle_plot(xyzs_file, radius_file, xlims, ylims, zlims, image_file);
    count = count + 1;
    close;
end
