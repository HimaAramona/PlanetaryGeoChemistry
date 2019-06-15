filename = 'output_';
curTime = datestr(now,'mm-dd-yyyy HH-MM'); 
suffix = '.mat';
filename = strcat(filename, curTime, suffix);
disp(filename);