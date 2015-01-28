%% Cleaning workspace & Setup
clear;
clc;

addpath(genpath('./simulation'));
addpath(genpath('./mm'));

if isunix()
    setenv('PATH', [getenv('PATH'), ':../build/bin']);
    if ismac()
        exec = @(cmd) system(sprintf('LD_LIBRARY_PATH="" && %s', cmd));
    else
        exec = @(cmd) system(sprintf('LD_LIBRARY_PATH="" && %s', cmd));
    end
else
    setenv('PATH', [getenv('PATH'), ';..\\build\\bin']);
    exec = @(cmd) system(cmd);
end

%% Simulating simple Poisson steps:
% SimParams class contains all necessary simulation parameters
% Default parameters can be changed by directly manipulating public member 
% variables: e.g. the number of datapoints: sp.N
% The following piece of code generates noisy step data:
sp = SimParams();
sp.N = 1e5;
sd = SimData(sp);
sd = sd.simulatePoissonSteps();
sd = sd.simOptTrapNoise();

mmwrite('noisy_data.mm', sd.data);

%% Getting Lambda_opt and denoise the signal
cmd = ['lambdaopt', ...
       ' --input ', 'noisy_data.mm'];
[status, out] = exec(cmd);
assert(status == 0);
lambda_opt = str2double(out);

cmd = ['denoising', ...
       ' --lambda ', num2str(lambda_opt), ...
       ' --input ', 'noisy_data.mm', ... 
       ' --output ', 'denoised_data.mm'];
[status] = exec(cmd);
assert(status == 0);

%% Building a level grid and running the graph cut
distance = 0.2;
cmd = ['level_generator', ...
       ' --level-distance ', num2str(distance), ...
       ' --input ', 'denoised_data.mm', ...
       ' --output ', 'level_data.mm'];
[status] = exec(cmd);
assert(status == 0);

cmd = ['graph_processing', ...
       ' --levels ', 'level_data.mm', ...
       ' --rho-d ', num2str(0.01), ...
       ' --rho-s ', num2str(0.1), ...
       ' --rho-p ', num2str(0.2), ...
       ' --prior-distance ', num2str(0.2), ...
       ' --input ', 'denoised_data.mm', ...
       ' --output ', 'clustered_data.mm'];
[status] = exec(cmd);
assert(status == 0);

%% Make a plot of the result
figure;
clf;
fontSize = 16;
fontName = 'Helvetica';

hold on;
plot1  = plot(sd.time, sd.data,'.', 'Color', [0.8,0.8,0.8], 'MarkerSize', 5);
plot2 = plot(sd.time, sd.pwcs, '.r', 'MarkerSize', 1, 'LineWidth', 1.5);
plot3 = plot(sd.time, mmread('denoised_data.mm'), '.g', 'MarkerSize', 1, 'LineWidth', 1.5);
plot4 = plot(sd.time, mmread('clustered_data.mm'), '.b', 'MarkerSize', 1, 'LineWidth', 1.5);
hold off;

title('Noisy, denoised and clustered data', ...
      'FontSize', fontSize, ...
      'FontName', fontName);
legend([plot1 plot2 plot3 plot4], ...
        'Noisy full bandwidth data set', ...
        'Pure step signal', ...
        'TVDN result', ...
        'Clustered result', ...
        'FontSize', floor(0.65 * fontSize), ...
        'FontName', fontName, ...
        'Location', 'NorthWest');
xlabel('time/s', 'FontSize', floor(0.85 * fontSize), 'FontName', fontName);
ylabel('bead motion/nm', 'FontSize', floor(0.85 * fontSize), 'FontName', fontName);
