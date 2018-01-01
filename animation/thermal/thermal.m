xyzs_file = 'xyzs.txt';
radius_file = 'radius.txt';
temperature_file = 'temperatures.txt';
temperature_plot(xyzs_file, radius_file, temperature_file);


steps = [0:2000:99000];
count = 1;
for step = steps
    temperature_file = ['temperatures_', num2str(step), '.txt'];
    disp(temperature_file);
    image_file = sprintf('image%05d.png', count);
    temperature_plot(xyzs_file, radius_file, temperature_file, image_file);
    count = count + 1;
    close
end