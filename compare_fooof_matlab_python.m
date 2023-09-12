% Compare that the fooof reimplementation in matlab is the same as in
% python
%
% Original python code by Tom Donoghue et. al in https://github.com/fooof-tools/fooof
%
% Cristina Gil Avila, TUM, 12.09.2023
close all
clear all

% Python fooof files
py_path = '../results/fooof/PFC';
py_files = dir(fullfile(py_path,'*.mat'));

% Matlab fooof files
mat_path = '../results/fooof_matlab';
mat_files = dir(fullfile(mat_path,'*.mat'));

n = length(py_files);


py = cell(1,n);
mat = cell(1,n);

r2_py = nan(1,n);
r2_mat = nan(1,n);
mae_py = nan(1,n);
mae_mat = nan(1,n);
npeaks_py = nan(1,n);
npeaks_mat = nan(1,n);
exp_py = nan(1,n);
exp_mat = nan(1,n);


for iFile =1:n
    
    py_file = py_files(iFile);
    load(fullfile(py_file.folder,py_file.name));
    
    py{iFile}.aperiodic_params = aperiodic_params;
    py{iFile}.peak_params = peak_params;
    py{iFile}.r2 = r_squared;
    py{iFile}.mae = error;
    
    r2_py(iFile) = r_squared;
    mae_py(iFile) = error;
    npeaks_py(iFile) = size(peak_params,1);
    exp_py(iFile) = aperiodic_params(2);

    clear aperiodic_params error gaussian_params peak_params r_squared
    
    
    mat_file = mat_files(iFile);
    load(fullfile(mat_file.folder,mat_file.name));
    
    mat{iFile}.aperiodic_params = self.aperiodic_params;
    mat{iFile}.peak_params = self.peak_params;
    mat{iFile}.r2 = self.r2;
    mat{iFile}.mae = self.error_mae;
    
    r2_mat(iFile) = self.r2;
    mae_mat(iFile) = self.error_mae;
    npeaks_mat(iFile) = size(self.peak_params,1);
    exp_mat(iFile) = self.aperiodic_params(2);
       
    clear self
    
end

% Compare structs
isequal(npeaks_py,npeaks_mat);

tol = 1e-5;
find(not(abs(r2_py-r2_mat) < tol))
find(not(abs(mae_py-mae_mat) < tol))
find(not(abs(exp_py-exp_mat) < tol))

py{7}.peak_params - mat{7}.peak_params