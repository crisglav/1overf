function [resTab] = ttestBF2(X,Y)
X = reshape(X, 1, length(X));
Y = reshape(Y, 1, length(Y));
if(isfile('tempTTest.csv'))
   delete('tempTTest.csv')
end
if(isfile('tempTTestRes.csv'))
   delete('tempTTestRes.csv')
end
if(length(X) > length(Y))
   datComb = [X',[Y';nan(length(X)-length(Y),1)]];
elseif(length(X) < length(Y))
    datComb = [[X';nan(length(Y)-length(X),1)],Y']; 
else
    datComb = [X',Y'];
end


tabXY = array2table(datComb, 'VariableNames', {'X', 'Y'});
writetable(tabXY, '/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface/tempTTest.csv')

% execute system command
system('export R_LIBS_USER=/rechenmagd4/Cristina/Rpackages/4.2.1; /sw/R/xenial/4.0.1/bin/Rscript /rechenmagd3/Experiments/2023_1overf/code/Rinterface/ttestBF24Matlab.R');

pause(3);

resTab = readtable('/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface/tempTTestRes.csv');

delete('/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface/tempTTest.csv')
delete('/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface/tempTTestRes.csv')
end