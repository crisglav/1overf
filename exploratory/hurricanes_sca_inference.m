% Hurricanes example
clear all
close all

% Hurricanes csv
filename = 'C:\Users\Mitarbeiter\Downloads\Specification Curve\Results for plotting\specifications_hurricanes.csv'; 
data = readtable(filename);

nSpec = height(data);