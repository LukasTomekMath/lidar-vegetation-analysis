clear all
clc

%% Data import and prep
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

% sentinel data
dataSent = readmatrix("112x72.csv");
dataMonoSent = readmatrix("mono_11x72_2018.csv");

% point cloud data
% dataRM = readmatrix("RM_curves_rev2.3.csv");
dataRM = readmatrix("RM_rev2.3_onlyVegHigh.csv");
% dataMonoRM = readmatrix("RM_curvesMono_inSVK.csv");
dataMonoRM = readmatrix("RM_mono_onlyVegHigh.csv");

% zaujimave kombinacie (iba vysoka vegetacia)
% 1) 7 16 18 (78.91 %)
% 2) 5 13 20 -> v 3D (91.2 %)
% 3) 7 18 20 -> v 3D (83.53 %) aj 2D
% 4) 2 5 13 20 -> v 3D (92.12 %)
index = @(i) ((i-1)*4 + 1):1:(i*4);

selected = [];
metrics = [ 5 13 20 ];
for i = metrics
	selected = [ selected index(i) ]; %#ok<AGROW> 
end

% data = dataSent;
data = dataRM(:, selected);
% data = [dataSent dataRM(:, selected)];

% dataMono = dataMonoSent;
dataMono = dataMonoRM(:, selected);
% dataMono = [dataMonoSent dataMonoRM(:, selected)];

data = [data; dataMono];

dataMin = min(data);
dataMax = max(data);
dataScaled = (data - dataMin) ./ (dataMax - dataMin);
dataScaled(isnan(dataScaled)) = 0;

dataMonoScaled = (dataMonoScaled - dataMin) ./ (dataMax - dataMin);
dataMonoScaled(isnan(dataMonoScaled)) = 0;
% dataMonoScaled = normalize(dataMonoScaled,"range");

[coef, score, ~, ~, explained, ~] = pca(dataScaled, "Algorithm","eig");
fprintf("explained: %.2f\n", sum(explained(1:3)));

% scoreMono = dataMonoScaled * coef;

n1 = 22;
n2 = 28;
n3 = 31;
n4 = 31;
n5 = 11;
N = n1 + n2 + n3 + n4 + n5;

C1start = 1;
C2start = n1 + 1;
C3start = n1 + n2 + 1;
C4start = n1 + n2 + n3 + 1;
C5start = n1 + n2 + n3 + n4 + 1;

C1end = n1;
C2end = n1 + n2;
C3end = n1 + n2 + n3;
C4end = n1 + n2 + n3 + n4;
C5end = n1 + n2 + n3 + n4 + n5;

ind = [C1start, C1end; C2start, C2end; C3start, C3end; C4start, C4end; C5start, C5end];

z = 3; % 26 29 52
X = [score(:,1)];%; scoreMono(:,1)];
Y = [score(:,2)];%; scoreMono(:,2)];
Z = [score(:,z)];%; scoreMono(:,z)];

% X = score(:,1);
% Y = score(:,2);

figure("Name","Z = " + num2str(z) + ". súradnica po PCA")
hold on
scatter(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), 20, "red", "filled");
scatter(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), 20, "green", "filled");
scatter(X(ind(3,1):ind(3,2)), Y(ind(3,1):ind(3,2)), 20, "blue", "filled");
scatter(X(ind(4,1):ind(4,2)), Y(ind(4,1):ind(4,2)), 20, "magenta", "filled");
scatter(X(ind(5,1):ind(5,2)), Y(ind(5,1):ind(5,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
title("Z = " + num2str(z) + ". súradnica po PCA")


figure("Name","Z = " + num2str(z) + ". súradnica po PCA")
scatter3(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), Z(ind(1,1):ind(1,2)), 20, "red", "filled");
hold on
scatter3(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), Z(ind(2,1):ind(2,2)), 20, "green", "filled");
scatter3(X(ind(3,1):ind(3,2)), Y(ind(3,1):ind(3,2)), Z(ind(3,1):ind(3,2)), 20, "blue", "filled");
scatter3(X(ind(4,1):ind(4,2)), Y(ind(4,1):ind(4,2)), Z(ind(4,1):ind(4,2)), 20, "magenta", "filled");
scatter3(X(ind(5,1):ind(5,2)), Y(ind(5,1):ind(5,2)), Z(ind(5,1):ind(5,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
title("Z = " + num2str(z) + ". súradnica po PCA")

%%
% Parallel coordinates
% INFO: vypocet priemerov a std
means = cell(5,1);
sigmas = cell(5,2);

for i = 1:4
	means{i} = mean( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,1} = means{i} - std( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,2} = means{i} + std( score( ind(i,1):ind(i,2), : ) );
end

means{5} = mean( scoreMono( 1:11, : ) );
sigmas{5,1} = means{5} - std( scoreMono( 1:11, : ) );
sigmas{5,2} = means{5} + std( scoreMono( 1:11, : ) );

% INFO: vykreslenie
figure; 
hold on

features = size(score,2);

% INFO: toto je tu len kvoli legende
plot([0,0],[0,0],'-r');
plot([0,0],[0,0],'-g');
plot([0,0],[0,0],'-b');
plot([0,0],[0,0],'-m');
plot([0,0],[0,0],'-k');

% 91E0
plot(1:1:features,  means{1}(1:features),  '.-r', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{1,1}(1:features),':r',  "LineWidth", 2)
plot(1:1:features, sigmas{1,2}(1:features),':r',  "LineWidth", 2)

% 91F0
plot(1:1:features,  means{2}(1:features),  '.-g', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{2,1}(1:features),':g',  "LineWidth", 2)
plot(1:1:features, sigmas{2,2}(1:features),':g',  "LineWidth", 2)

% 91G0
plot(1:1:features,  means{3}(1:features),  '.-b', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{3,1}(1:features),':b',  "LineWidth", 2)
plot(1:1:features, sigmas{3,2}(1:features),':b',  "LineWidth", 2)

% 9110
plot(1:1:features,  means{4}(1:features),  '.-m', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{4,1}(1:features),':m',  "LineWidth", 2)
plot(1:1:features, sigmas{4,2}(1:features),':m',  "LineWidth", 2)

% monokultury
plot(1:1:features,  means{5}(1:features),  '.-k', "LineWidth", 3, "MarkerSize", 20)
plot(1:1:features, sigmas{5,1}(1:features),':k',  "LineWidth", 1.5)
plot(1:1:features, sigmas{5,2}(1:features),':k',  "LineWidth", 1.5)

legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
hold off


% figure; 
% hold on
% % features = length(selected);
% features = 72;
% 
% for i = C1start:C1end 
% 	plot(1:1:features, score(i,1:features), '.-r');
% end
% for i = C2start:C2end 
% 	plot(1:1:features, score(i,1:features), '.-g');
% end
% for i = C3start:C3end
% 	plot(1:1:features, score(i,1:features), '.-b');
% end
% for i = C4start:C4end 
% 	plot(1:1:features, score(i,1:features), '.-m');
% end
% for i = 1:11
% 	plot(1:1:features, score(i,1:features), '.-k'); 
% end
% % legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
% hold off






