vdatapath = '/rechenmagd3/Experiments/2023_1overf/results/features/vdata/';
files = dir([vdatapath '*.mat']);
if(isempty(gcp('nocreate')))
    parObj = parpool(30);
end
% for iRand = 1:500
%     for iSpec =1:48
        parfor iSubj=1:length(files)
            aux = load(fullfile(files(iSubj).folder,files(iSubj).name));
            fprintf(['Loaded ' files(iSubj).name '\n']);
%             save(fullfile(files(iSubj).folder,files(iSubj).name),'vdata');
%             fprintf(['Saved ' files(iSubj).name '\n']);
%             clear vdata;
        end
%     end
% end