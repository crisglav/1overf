function params = load_params(params_path)
% params_path: path to the discover-eeg saved parameters
fid = fopen(params_path);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
params = jsondecode(str);

% Check that the preprocess the data exists and has not been moved
if ~exist(params.PreprocessedDataPath)
    warning('The data might have been moved. Using relative paths.')
    
    params.RawDataPath = '../data/blinded';
    params.PreprocessedDataPath = '../data/blinded/derivatives_v2023_08_18';
end

end

