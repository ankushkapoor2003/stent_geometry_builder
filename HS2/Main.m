function Main(x)
% Saving the current wd path
original_dir = pwd;

% Setting up the results directory
result_dir_name = 'HS2_Results';
result_dir_path = fullfile(original_dir, result_dir_name);

% Generating the results directory if not already created
if ~exist(result_dir_path, 'dir')
    mkdir(result_dir_path);
end

% Call the stent generation function
HS2_Generation(x, result_dir_path);
end