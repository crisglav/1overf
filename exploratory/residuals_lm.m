close all;

rand1 = rand(150,1);
rand2 = rand(150,1);
rand3 = rand(150,1);
age = 30 + 2*rand1 + rand(150,1);

% figure, scatter(age,rand1);
% figure, scatter(age,rand2);

model_rand1 = fitlm(age,rand1);
model_rand2 = fitlm(age,rand2);

figure, plot(model_rand1)
figure, plot(model_rand2)

model_res = fitlm(model_rand1.Residuals.Raw, model_rand2.Residuals.Raw);
figure, plot(model_res);

model_rand1 = fitlm(rand3,rand1);
model_rand2 = fitlm(rand3,rand2);
figure, plot(model_rand1)
figure, plot(model_rand2)
model_res = fitlm(model_rand1.Residuals.Raw, model_rand2.Residuals.Raw);
figure, plot(model_res);