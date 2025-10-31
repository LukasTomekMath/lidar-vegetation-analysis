% close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
CURVES_DIRECTORY = fullfile(DATA_DIRECTORY, 'curves');

SHAPEFILES_DIRECTORY = fullfile(DATA_DIRECTORY,'_shapeFiles');
SHP_ALL_LOTS = fullfile(SHAPEFILES_DIRECTORY,...
						'prehlad_lokalit_lls_1_cyklus/prehlad_lokalit_lls_1_cyklus.shp');

SHP_FOOTPRINTS = cell(42, 1);
MAT_FOOTPRINTS = cell(42, 1);
lot_i = "";
for i = 1:42
	if (i <= 9)
		lot_i = "LOT0" + num2str(i);
	else
		lot_i = "LOT" + num2str(i);
	end
	
	SHP_FOOTPRINTS{i} = fullfile(SHAPEFILES_DIRECTORY,...
							  "footprints_s-jtsk03_shp/" + lot_i + "/" + lot_i + "_las.shp");

	MAT_FOOTPRINTS{i} = fullfile(SHAPEFILES_DIRECTORY,"footprints_mat/" + num2str(i) + ".mat");

end

LOT_DIR = fullfile(MAIN_DIRECTORY, "data/lazFiles");

cd(MAIN_DIRECTORY);

addpath(MAIN_DIRECTORY);

lasColorMap = load('las_colormap.txt')/255; % color map for classified point cloud

COLOR_GROUND     = lasColorMap(2,:);
COLOR_VEG_LOW    = lasColorMap(3,:);
COLOR_VEG_MEDIUM = lasColorMap(4,:);
COLOR_VEG_HIGH   = lasColorMap(5,:);

%%
% LOAD .KML CURVES
[KMLFiles, KMLFilePath] = uigetfile('*.kml', 'Select KML file', CURVES_DIRECTORY,'MultiSelect','off');
KMLFiles = string(KMLFiles);
KMLFilePath = string(KMLFilePath);
tic
curves = cell(numel(KMLFiles), 1);

for i = 1:numel(curves)
	coords = readKML(KMLFilePath + KMLFiles(i));
	if isempty(coords)
		continue;
	end
	lon = coords(:,1)';
	lat = coords(:,2)';
	h = coords(:,3)';
	
	[x, y] = gps_to_JTSK03_transformation(lat, lon, h);
	curves{i}.poly = polyshape(x, y);
	curves{i}.name = KMLFiles(i);

	fprintf("Curve %d/%d loaded\n", i, numel(curves));
end
lazTime = toc;
clear h x y lat lon coords

fprintf("All KML curves loaded in %.2f s.\n", lazTime);

cd(MAIN_DIRECTORY);
h = 10; n = 1;

dataPP = dataPreprocessorCurve(curves{1}.poly, h, n);

% FIND NAMES OF .LAZ FILES FOR EACH LOADED CURVE
foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
foundLazFiles = cell(length(curves), 1);
LOTS_curves = loadLOTS(SHP_ALL_LOTS); % table of curves for
fprintf("LOTs curves loaded\n");

tic
for i = 1:length(curves) % iterate through all .KML curves
	foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
	% find indices of LOTs where the i-th .KML curve belongs to
	for j = 1:length(LOTS_curves) % iterate through all LOT curves
		% find intersection of i-th .KML curve with j-th LOT
		intersection = intersect(dataPP.Omega, LOTS_curves{j}.poly);

		if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
			continue;
		end

		foundLots = [foundLots, str2num(LOTS_curves{j}.lotNum)];

	end % end j
	
	% find names of .LAZ files to load for the i-th curve
	foundLazFiles{i} = {};
	lazFileFound = 0;
	FOOTPRINTS_curves = [];
	
	if isempty(foundLots)
		fprintf("Curve %s outside SVK\n", curves{i}.name);
		continue;
	end
	
	for j = 1:length(foundLots) % iterate over all found LOTs
% 		FOOTPRINTS_curves = loadFOOTPRINTS(SHP_FOOTPRINTS{foundLots(j)});
		load(MAT_FOOTPRINTS{foundLots(j)});
% 		fprintf("Footprints loaded\n");
		
		for k = 1:length(FOOTPRINTS_curves) % iterate over all footprints within found j-th LOT

			intersection = intersect(dataPP.Omega, FOOTPRINTS_curves{k}.poly);

			if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
				continue;
			end
		
			% save found .LAZ file name
			lazFileFound = lazFileFound + 1;
			foundLazFiles{i}{lazFileFound}.poly = FOOTPRINTS_curves{k}.poly;
			foundLazFiles{i}{lazFileFound}.lazName = FOOTPRINTS_curves{k}.fileName;
			foundLazFiles{i}{lazFileFound}.lazPath = strrep(fullfile(LOT_DIR, FOOTPRINTS_curves{k}.fileName),'\', '/');
		end
		
	end

	fprintf("Curve %d/%d done\n", i, length(curves));
% 	fprintf("Found .LAZ files:\n");
% 	for j = 1:length(foundLazFiles{i})
% 		fprintf("%s\n", foundLazFiles{i}{j});
% 	end
end % end i

lazTime = toc;
fprintf("LAZ files for all curves found, elapsed time: %.2f s\n", lazTime);

% PRINT FOUND .LAZ FILES
for i = 1:length(foundLazFiles)
	fprintf("\nFound .LAZ files for curve '%s':\n", curves{i}.name);
	for j = 1:length(foundLazFiles{i})
	fprintf("%s\n", foundLazFiles{i}{j}.lazName);
	end
end

%% PLOT KML CURVE AND FOUND FOOTPRINTS
% plot(curves{2},"LineWidth",5,"EdgeColor","red")
c = 1;
figure
axis equal
hold on
plot(dataPP.Omega,"LineWidth",4,"EdgeColor","magenta")
for i = 1:length(foundLazFiles{c})
	plot(foundLazFiles{c}{i}.poly,"FaceAlpha",0.05);
	pause(0.5)
end
plot(curves{c}.poly,"LineWidth",4,"EdgeColor","red","FaceAlpha",0)
hold off

%% Create out file for C++
outFile = fopen(curves{1}.name + "_data.txt", 'w');

fprintf(outFile, "%d;%d\n", dataPP.nx, dataPP.ny);
fprintf(outFile, "%d;%d\n", dataPP.nx/h, dataPP.ny/h);
fprintf(outFile, "%.10f;%.10f\n", dataPP.x1, dataPP.y2);
fprintf(outFile, "%d\n", length(foundLazFiles{c}));

for i = 1:length(foundLazFiles{c}) % nacitanie PC cez .MAT subory pre krivku "c"
	fprintf(outFile, "%s\n", foundLazFiles{c}{i}.lazPath);
end

fclose(outFile);

%% Visualize computed metrics
cd(DATA_DIRECTORY);

filename = '2024_02_25_18-14-35_91E0_bodikyprihradzi1_final.kml_data_h=10m.tif';
[metrics, rasterReference] = readgeoraster(filename, 'OutputType', 'double');
infoMetrics = georasterinfo(filename);

metrics = standardizeMissing(metrics, infoMetrics.MissingDataIndicator);
metrics = flipud(metrics);
rasterReference.ColumnsStartFrom = 'south';

[X, Y] = rasterReference.worldGrid();

band = 5;
figure('Name', 'Metrics Bodiky')
title 'Bodiky data'
alphaData = ones(size(X));
alphaData(isnan(metrics(:,:,band))) = 0;
imagesc(X(1,:), Y(:,1), metrics(:,:,band), 'AlphaData', alphaData)
set(gca,'YDir','normal') % spravne hodnoty na y osi
colormap bone
colorbar
axis equal
% plotCurvesPolyshape(curves)

%% Find all necessarry laz files for all curves
uniqueLaz = cell(1);
uniqueLaz{1} = '';
unique = 1;
for i = 1:length(foundLazFiles)
	for j = 1:length(foundLazFiles{i})
		if (~ismember(foundLazFiles{i}{j}.lazName, uniqueLaz))
			uniqueLaz{unique} = foundLazFiles{i}{j}.lazName;
			unique = unique + 1;
		end

	end
end

uniqueLaz = uniqueLaz';

%% Find missing laz files inside laz files folder
cd(LOT_DIR)

missingFiles = cell(1);
missingFiles{1} = '';
found = 1;
for i = 1:length(uniqueLaz)
	if (~isfile(uniqueLaz{i})) % if laz not found
		missingFiles{found} = uniqueLaz{i};
		found = found + 1;
	end
end

missingFiles = missingFiles';

%% Copy missing files - current folder path to hard-drive
for i = 1:length(missingFiles)
	fprintf("file %d/%d\n", i, length(missingFiles));
	if (isfile(missingFiles{i}))
		if (~isfile(fullfile(MAIN_DIRECTORY, "data/lazFiles", missingFiles{i})))
			copyfile(missingFiles{i}, fullfile(MAIN_DIRECTORY, "data/lazFiles", missingFiles{i}))
		end
	else
		fprintf("file %d not found\n", i);
	end
end

fprintf("Copy files done\n");

% cd(MAIN_DIRECTORY)

%% Vytvorenie stvorcov na paralelny vypocet
LOTS_curves = loadLOTS(SHP_ALL_LOTS);

%%
cd(MAIN_DIRECTORY)

lots_union_1 = LOTS_curves{1}.poly;
lots_union_2 = LOTS_curves{22}.poly;

for i = 2:21
	fprintf("lot %d/%d\n", i, length(LOTS_curves))
	lots_union_1 = union(lots_union_1, LOTS_curves{i}.poly);
end

for i = 23:42
	fprintf("lot %d/%d\n", i, length(LOTS_curves))
	lots_union_2 = union(lots_union_2, LOTS_curves{i}.poly);
end

%%
lots_union = union(lots_union_1, lots_union_2);

%%
figure
plot(lots_union,'FaceAlpha',0.5,'FaceColor','green')
hold on
axis equal


[BB_x, BB_y] = lots_union.boundingbox;

BBx = [BB_x(1), BB_x(2), BB_x(2), BB_x(1)];
BBy = [BB_y(2), BB_y(2), BB_y(1), BB_y(1)];

BB = polyshape(BBx, BBy);
% plot(BB, 'FaceAlpha',0.1,'FaceColor','black','EdgeColor','black')

h = 2000; n = 0;

dataFE = dataFeatureExtractorCurve(lots_union, h, n,{},{});

plot(dataFE.Omega,'FaceAlpha',0,'EdgeColor','magenta')

dataFE.plotMesh("PlotCenters",false)

hold off
% 588 obrazkov
%% stvorce pre cele SVK
index = 0;
x = [];
y = [];
squares = {};
Yc = flip(dataFE.yc);
Xc = dataFE.xc;

xc = Xc(1);
yc = Yc(1);

x = [xc - h/2 ,xc + h/2 ,xc + h/2 ,xc - h/2];
y = [yc - h/2, yc - h/2, yc + h/2, yc + h/2];

diff_x = abs(x(1)) - floor(abs(x(1)));
diff_y = abs(y(1)) - floor(abs(y(1)));

s = 1;
for i = 1:length(Yc)
	for j = 1:length(Xc)
		index = index + 1;
		xc = Xc(j);
		yc = Yc(i);

		x = [xc - h/2, xc + h/2, xc + h/2, xc - h/2];
		y = [yc + h/2, yc + h/2, yc - h/2, yc - h/2];
		
		x = x + diff_x;
		y = y + diff_y;

		temp = polyshape(x,y);

		squares{s}.name = num2str(i-1) + "-" + num2str(j-1) + "_area";
		squares{s}.poly = temp;
		squares{s}.index = index;

		s = s + 1;
	end
end

squares = squares';

%% stvorce pre loty 2 a 3 do diplomovky
% index = 0;
% x = [];
% y = [];
% squares = {};
% Yc = flip(dataFE.yc);
% Xc = dataFE.xc;
% s = 1;
% for i = 1:length(Yc)
% 	for j = 1:length(Xc)
% 		index = index + 1;
% 		xc = Xc(j);
% 		yc = Yc(i);
% 
% 		x = [xc - h/2 ,xc + h/2 ,xc + h/2 ,xc - h/2];
% 		y = [yc - h/2, yc - h/2, yc + h/2, yc + h/2];
% 		temp = polyshape(x,y);
% 
% 		intersection = intersect(temp, lots_2_3);
% 
% 		if isempty(intersection.Vertices)
% 			continue;
% 		end
% 
% 		intersection = intersect(temp, LOTS_curves{1}.poly);
% 
% 		if ~isempty(intersection.Vertices)
% 			continue;
% 		end
% 
% 		squares{s}.name = num2str(i-1) + "-" + num2str(j-1) + "_area";
% 		squares{s}.poly = temp;
% 		squares{s}.index = index;
% 
% 		s = s + 1;
% 	end
% end
% 
% squares = squares';

%%
figure
plot(lots_union,'FaceAlpha',0.5,'FaceColor','green')
axis equal
hold on
for i = length(squares):length(squares)%length(squares)
	plot(squares{i}.poly, 'EdgeColor','black','FaceAlpha',0, 'LineWidth',1)
% 	[xCent, yCent] = squares{i}.poly.centroid;
% 	text(xCent-500,yCent,num2str(squares{i}.index),"FontSize",15);
% 	pause(0.001)
end
xline(-164375)
xline(-166375)
xline(-592375)
xline(-590375)
yline(-1131730)
yline(-1133730)
yline(-1333730)
yline(-1335730)

hold off

%% export squares for whole SVK
out = fopen("_grid_squares.txt","w");

for i = 1:length(squares)
	if i == length(squares)
		fprintf(out, "%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f",...
		squares{i}.poly.Vertices(1,1),squares{i}.poly.Vertices(1,2),...
		squares{i}.poly.Vertices(2,1),squares{i}.poly.Vertices(2,2),...
		squares{i}.poly.Vertices(3,1),squares{i}.poly.Vertices(3,2),...
		squares{i}.poly.Vertices(4,1),squares{i}.poly.Vertices(4,2));
	else
		fprintf(out, "%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\n",...
			squares{i}.poly.Vertices(1,1),squares{i}.poly.Vertices(1,2),...
			squares{i}.poly.Vertices(2,1),squares{i}.poly.Vertices(2,2),...
			squares{i}.poly.Vertices(3,1),squares{i}.poly.Vertices(3,2),...
			squares{i}.poly.Vertices(4,1),squares{i}.poly.Vertices(4,2));
	end
end

fclose(out);

%% Find LAZ files for all squares
foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
foundLazFiles = cell(length(squares), 1);

tic
for i = 1:length(squares) % iterate through all .KML curves
	foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
	% find indices of LOTs where the i-th .KML curve belongs to
	for j = 1:length(LOTS_curves) % iterate through all LOT curves
		% find intersection of i-th .KML curve with j-th LOT
		intersection = intersect(squares{i}.poly, LOTS_curves{j}.poly);

		if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
			continue;
		end

		foundLots = [foundLots, str2num(LOTS_curves{j}.lotNum)];

	end % end j
	
	% find names of .LAZ files to load for the i-th curve
	foundLazFiles{i} = {};
	lazFileFound = 0;
	FOOTPRINTS_curves = [];
	
	if isempty(foundLots)
		fprintf("Square %d outside SVK\n", i);
		continue;
	end
	
	for j = 1:length(foundLots) % iterate over all found LOTs
% 		FOOTPRINTS_curves = loadFOOTPRINTS(SHP_FOOTPRINTS{foundLots(j)});
		load(MAT_FOOTPRINTS{foundLots(j)});
% 		fprintf("Footprints loaded\n");
		
		for k = 1:length(FOOTPRINTS_curves) % iterate over all footprints within found j-th LOT

			intersection = intersect(squares{i}.poly, FOOTPRINTS_curves{k}.poly);

			if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
				continue;
			end
		
			% save found .LAZ file name
			lazFileFound = lazFileFound + 1;
			foundLazFiles{i}{lazFileFound}.poly = FOOTPRINTS_curves{k}.poly;
			foundLazFiles{i}{lazFileFound}.lazName = FOOTPRINTS_curves{k}.fileName;
			foundLazFiles{i}{lazFileFound}.lazPath = strrep(fullfile(LOT_DIR, FOOTPRINTS_curves{k}.fileName),'\', '/');
		end
		
	end

	fprintf("Square %d/%d done\n", i, length(squares));

end % end i

lazTime = toc;
fprintf("LAZ files for all curves found, elapsed time: %.2f s\n", lazTime);

%%
c = 2;
figure
axis equal
hold on
plot(squares{c}.poly,"LineWidth",5,"EdgeColor","red")
for i = 1:length(foundLazFiles{c})
	plot(foundLazFiles{c}{i}.poly,"FaceAlpha",0.05);
	pause(0.5)
end
plot(squares{c}.poly,"LineWidth",5,"EdgeColor","red","FaceAlpha",0)
hold off

%%
h = 10; n = 0;
dataPP = dataPreprocessorCurve(squares{373}.poly, h, n);

%% Create out files for C++
cd(MAIN_DIRECTORY)

fprintf("Creating outFiles... ");
outFilesAll = fopen("outFiles/_allOutFiles.txt", 'w');
fprintf(outFilesAll, "%d\n", length(squares));
for c = 1:length(squares)
	outFile = fopen("outFiles/" + num2str(squares{c}.index) + "_" + squares{c}.name + "_data.txt", 'w');
	fprintf(outFilesAll, num2str(squares{c}.index) + "_" + squares{c}.name + "_data.txt\n");

	h = 10; n = 0;
	dataPP = dataPreprocessorCurve(squares{c}.poly, h, n);
	
	fprintf(outFile, "%d;%d\n", dataPP.nx, dataPP.ny);
	fprintf(outFile, "%d;%d\n", dataPP.nx/h, dataPP.ny/h);
	fprintf(outFile, "%.10f;%.10f\n", dataPP.x1, dataPP.y2);
	fprintf(outFile, "%d\n", length(foundLazFiles{c}));

	for i = 1:length(foundLazFiles{c}) % nacitanie PC cez .MAT subory pre krivku "c"
		fprintf(outFile, "%s\n", foundLazFiles{c}{i}.lazPath);
	end

	fclose(outFile);
	fprintf("OutFile %d/%d done\n", c, length(foundLazFiles));
end

fclose(outFilesAll);
fprintf("done\n");

%%
% LAZ_count; nPts; read_time; norm_time; redist_time; comp_metrics_time;
data = readmatrix('stats.txt');
data(:,3) = data(:,3) / 60; % convert seconds to minutes

means = mean(data);

%%
figure % LAZ counts
boxplot(data(:, 1))

figure % Pts count
boxplot(data(:,2));

%%
types = ["Normalization","Redistribution","Metrics computation"];
X = categorical(types);
X = reordercats(X, types);

figure
bar(X, means(4:6))

figure
boxplot(data(:,4:6),'Labels', {'Normalization','Redistribution','Metrics computation'})
hold on
scatter([1 2 3], means(4:6),50, [0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5],...
					"filled","o","MarkerEdgeColor","k","LineWidth",2)






