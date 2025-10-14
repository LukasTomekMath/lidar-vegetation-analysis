#include "PointCloudMetrics.h"

DataHandler::DataHandler()
{
}

// TODO: pridat metodu na initalizovanie `m_areaInfo` premennej

DataHandler::~DataHandler()
{
	std::cout << "Destructor\n";

	for (size_t i = 0; i < m_metrics.size(); i++)
	{
		delete m_metrics[i];
	}
}

void DataHandler::setupAreaInfo(const std::string& lazFileName, const double upperLeftX, const double upperLeftY)
{
	m_areaInfo.width_n = 2000;
	m_areaInfo.height_n = 2000;
	m_areaInfo.width = 1000;
	m_areaInfo.height = 1000;
	m_areaInfo.desiredPixelSize = m_areaInfo.width_n / m_areaInfo.width;;
	m_areaInfo.xLeft = upperLeftX;
	m_areaInfo.yLeft = upperLeftY;
	m_areaInfo.lazFileName = lazFileName;
}

bool DataHandler::performCalculation()
{
	// std::cout << "width_n: " << m_areaInfo.width_n << "\nheight_n: " << m_areaInfo.height_n << "\n";
	// std::cout << "width: " << m_areaInfo.width << "\nheight: " << m_areaInfo.height << "\n";
	// std::cout << "(x,y): [" << m_areaInfo.xLeft << ", " << m_areaInfo.yLeft << "]\n";

	auto reading_start = std::chrono::high_resolution_clock::now();

	if (!readLazFile())
	{
		std::cout << "failed to read '" << m_areaInfo.lazFileName << "' file, computation will not proceed further" << std::endl;
		return false;
	}

	auto reading_end = std::chrono::high_resolution_clock::now();
	m_times.readTime = std::chrono::duration_cast<std::chrono::milliseconds>(reading_end - reading_start).count() / 1000.0;

	auto normalization_start = std::chrono::high_resolution_clock::now();
	normalizePoints();
	auto normalization_end = std::chrono::high_resolution_clock::now();
	m_times.normalizationTime = std::chrono::duration_cast<std::chrono::milliseconds>(normalization_end - normalization_start).count() / 1000.0;

	auto redistribution_start = std::chrono::high_resolution_clock::now();
	redistributePoints();
	auto redistribution_end = std::chrono::high_resolution_clock::now();
	m_times.redistributionTime = std::chrono::duration_cast<std::chrono::milliseconds>(redistribution_end - redistribution_start).count() / 1000.0;

	//computation:
	auto metrics_start = std::chrono::high_resolution_clock::now();
	computeMetrics();
	auto metrics_end = std::chrono::high_resolution_clock::now();
	m_times.metricsTime = std::chrono::duration_cast<std::chrono::milliseconds>(metrics_end - metrics_start).count() / 1000.0;

	auto export_start = std::chrono::high_resolution_clock::now();
	exportMetrics(m_areaName);
	auto export_end = std::chrono::high_resolution_clock::now();
	m_times.exportTime = std::chrono::duration_cast<std::chrono::milliseconds>(export_end - export_start).count() / 1000.0;


	return true;
}

bool DataHandler::readLazFile()
{
	LASreadOpener lasReadOpener;
	LASreader* reader = nullptr;
	std::string fileName = "";
	int nPointsInMesh = 0;

	int nPixels = m_areaInfo.width_n * m_areaInfo.height_n;
	m_meshPixels.resize(nPixels, PixelPoints());
	//m_DTM = new double[nPixels] { 0.0 };

	auto start = std::chrono::high_resolution_clock::now();

	// TODO: tento for cyklus sa moze dat prec a nacita sa iba jeden LAZ subor, ktoreho nazov by mohol byt ulozeny v `m_areaInfo`, inak kod ohladom citania a ukladania hodnot netreba menit
	lasReadOpener.add_file_name(m_areaInfo.lazFileName.c_str());

	if (!lasReadOpener.active())
	{
		std::cout << "Laz file " << m_areaInfo.lazFileName << "not active\n";
		return false;
	}

	fileName = lasReadOpener.get_file_name();
	std::cout << "trying to read file " << m_areaInfo.lazFileName << std::endl;

	reader = lasReadOpener.open();
	if (reader == nullptr)
	{
		std::cout << "Failed to open laz file" << m_areaInfo.lazFileName << std::endl;
		return false;
	}

	//printf("xBounds: [%.2lf, %.2lf]\n", reader->get_min_x(), reader->get_max_x());
	//printf("yBounds: [%.2lf, %.2lf]\n", reader->get_min_y(), reader->get_max_y());
	//printf("zBounds: [%.2lf, %.2lf]\n", reader->get_min_z(), reader->get_max_z());
	//std::cout << "nPoints: " << reader->npoints << "\n";

	double xLeft = m_areaInfo.xLeft;
	double yLeft = m_areaInfo.yLeft;
	double ptX = 0, ptY = 0, ptZ = 0;
	int ptClass = 0;
	double deltaX = 0, deltaY = 0;
	int i = -1; // row index
	int j = -1; // column index
	int pixelIdx = -1;

	//printf("Scale factors: [%.6lf, %.6lf, %.6lf]\n", reader->header.x_scale_factor, reader->header.y_scale_factor, reader->header.z_scale_factor);

	while (reader->read_point())
	{
		ptX = reader->point.get_x();
		ptY = reader->point.get_y();
		ptZ = reader->point.get_z();
		ptClass = reader->point.get_classification();

		// check if point is from class 2, 3, 4 or 5 (ground or vegetation)
		if (ptClass < GROUND || ptClass > VEG_HIGH)
			continue;

		deltaX = ptX - xLeft;
		deltaY = ptY - yLeft;

		i = static_cast<int>(std::floor(std::abs(deltaY)));
		j = static_cast<int>(std::floor(std::abs(deltaX)));

		// check if point is inside mesh
		if (deltaY > 0.0 || i > m_areaInfo.height_n - 1)
			continue;
		
		if (deltaX < 0.0 || j > m_areaInfo.width_n - 1)
			continue;

		// save point to the corresponding mesh pixel
		pixelIdx = i * m_areaInfo.width_n + j;

		if (ptClass == GROUND)
		{
			m_meshPixels[pixelIdx].groundPts.xCoords.push_back(ptX);
			m_meshPixels[pixelIdx].groundPts.yCoords.push_back(ptY);
			m_meshPixels[pixelIdx].groundPts.zCoords.push_back(ptZ);

			nPointsInMesh++;

			continue;
		}

		m_meshPixels[pixelIdx].vegetationPts.xCoords.push_back(ptX);
		m_meshPixels[pixelIdx].vegetationPts.yCoords.push_back(ptY);
		m_meshPixels[pixelIdx].vegetationPts.zCoords.push_back(ptZ);
		m_meshPixels[pixelIdx].vegetationPts.ptClass.push_back(ptClass);

		nPointsInMesh++;
	}

	reader->close();
	delete reader;

	std::cout << "nPointsInMesh: " << nPointsInMesh << "\n";
	
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	printf("\nReading LAZ files: %.4lf s\n", (double)duration.count() / 1000.0);

	return true;
}

// ? TODO: toto by sa mozno mohlo pouzit namiesto upravovania `readLazFiles` metody, kedze by to mal but rovnaky kod
// bool DataHandler::readLasFile(std::string lasFileName)
// {
// 	LASreadOpener lasReadOpener;
// 	LASreader* reader = nullptr;
// 	int nPointsInMesh = 0;

// 	int nPixels = m_areaInfo.width_n * m_areaInfo.height_n;
// 	m_meshPixelsRedistributed.resize(nPixels, PixelPoints());
// 	int desiredPixelSize = m_areaInfo.desiredPixelSize;

// 	auto start = std::chrono::high_resolution_clock::now();

// 	lasReadOpener.add_file_name(lasFileName.c_str());

// 	if (!lasReadOpener.active())
// 	{
// 		std::cout << "Las file " << lasFileName << "not active\n";
// 		return false;
// 	}

// 	std::cout << "reading file " << lasFileName << "\n";

// 	reader = lasReadOpener.open();
// 	if (reader == nullptr)
// 	{
// 		std::cout << "Failed to open las file\n";
// 		return false;
// 	}

// 	//printf("xBounds: [%.2lf, %.2lf]\n", reader->get_min_x(), reader->get_max_x());
// 	//printf("yBounds: [%.2lf, %.2lf]\n", reader->get_min_y(), reader->get_max_y());
// 	//printf("zBounds: [%.2lf, %.2lf]\n", reader->get_min_z(), reader->get_max_z());
// 	//std::cout << "nPoints: " << reader->npoints << "\n";

// 	double xLeft = m_areaInfo.xLeft;
// 	double yLeft = m_areaInfo.yLeft;
// 	double ptX = 0, ptY = 0, ptZ = 0;
// 	int ptClass = 0;
// 	double deltaX = 0, deltaY = 0;
// 	int i = -1; // row index
// 	int j = -1; // column index
// 	int pixelIdx = -1;
// 	int nPoints = 0;

// 	//printf("Scale factors: [%.6lf, %.6lf, %.6lf]\n", reader->header.x_scale_factor, reader->header.y_scale_factor, reader->header.z_scale_factor);
// 	while (reader->read_point())
// 	{
// 		nPoints++;
// 		ptX = reader->point.get_x();
// 		ptY = reader->point.get_y();
// 		ptZ = reader->point.get_z();
// 		ptClass = reader->point.get_classification();

// 		// check if point is from class 2, 3, 4 or 5 (ground or vegetation)
// 		if (ptClass < GROUND || ptClass > VEG_HIGH)
// 			continue;

// 		// std::cout << "point class: " << ptClass << std::endl;

// 		deltaX = ptX - xLeft;
// 		deltaY = ptY - yLeft;

// 		i = static_cast<int>(std::floor(std::abs(deltaY)));
// 		j = static_cast<int>(std::floor(std::abs(deltaX)));

// 		// check if point is inside mesh
// 		if (deltaY > 0.0 || i >= m_areaInfo.height_n - 1)
// 			continue;

// 		if (deltaX < 0.0 || j >= m_areaInfo.width_n - 1)
// 			continue;

// 		// save point to the corresponding mesh pixel
// 		pixelIdx = i * m_areaInfo.width + j;
// 		//std::cout << "idx: " << pixelIdx << "\n";

// 		if (ptClass == GROUND)
// 		{
// 			m_meshPixelsRedistributed[pixelIdx].groundPts.xCoords.emplace_back(ptX);
// 			m_meshPixelsRedistributed[pixelIdx].groundPts.yCoords.emplace_back(ptY);
// 			m_meshPixelsRedistributed[pixelIdx].groundPts.zCoords.emplace_back(ptZ);

// 			nPointsInMesh++;

// 			continue;
// 		}

// 		m_meshPixelsRedistributed[pixelIdx].vegetationPts.xCoords.emplace_back(ptX);
// 		m_meshPixelsRedistributed[pixelIdx].vegetationPts.yCoords.emplace_back(ptY);
// 		m_meshPixelsRedistributed[pixelIdx].vegetationPts.zCoords.emplace_back(ptZ);
// 		m_meshPixelsRedistributed[pixelIdx].vegetationPts.ptClass.emplace_back(ptClass);

// 		nPointsInMesh++;
// 	}

// 	reader->close();
// 	delete reader;

// 	std::cout << "nPointsInMesh: " << nPointsInMesh << ", total points read: " << nPoints << std::endl;

// 	auto end = std::chrono::high_resolution_clock::now();
// 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

// 	printf("\nReading LAS file: %.4lf s\n", (double)duration.count() / 1000.0);

// 	return true;
// }

void DataHandler::normalizePoints()
{
	double minZ = -1.0;
	double minV = -1.0;
	double minG = -1.0;

	bool noGroundPts = true;
	bool noVegetationPts = true;

	auto start = std::chrono::high_resolution_clock::now();
	int i = 0;
#pragma omp parallel for private(minZ, minV, minG, noGroundPts, noVegetationPts)
	for (i = 0; i < m_meshPixels.size(); i++)
	{
		//m_DTM[i] = NAN;

		// Check if there are points in pixel
		noGroundPts = m_meshPixels[i].groundPts.zCoords.empty();
		noVegetationPts = m_meshPixels[i].vegetationPts.zCoords.empty();
		
		if (noGroundPts && noVegetationPts)
			continue;

		if (noGroundPts)
			minG = DBL_MAX;
		else
			minG = *std::min_element(m_meshPixels[i].groundPts.zCoords.begin(), m_meshPixels[i].groundPts.zCoords.end());

		if (noVegetationPts)
			minV = DBL_MAX;
		else
			minV = *std::min_element(m_meshPixels[i].vegetationPts.zCoords.begin(), m_meshPixels[i].vegetationPts.zCoords.end());

		minZ = (minG < minV) ? minG : minV;

		//m_DTM[i] = StatFunctions::mean(m_meshPixels[i].groundPts.zCoords);

		for (auto& z : m_meshPixels[i].vegetationPts.zCoords)
		{
			z -= minZ;
		}

		for (auto& z : m_meshPixels[i].groundPts.zCoords)
		{
			z -= minZ;
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	printf("\nNormalization: %.4lf s\n", (double)duration.count() / 1000.0);

}

void DataHandler::redistributePoints()
{
	int nPixels = m_areaInfo.width * m_areaInfo.height;
	m_meshPixelsRedistributed.resize(nPixels, PixelPoints());

	int desiredPixelSize = m_areaInfo.desiredPixelSize;
	int I = -1, J = -1;
	int index = -1;
	int INDEX = -1;
	int j = -1, k = -1, l = -1;

	std::vector<double>* v1 = nullptr;
	std::vector<double>* v2 = nullptr;

	std::vector<int>* v1i = nullptr;
	std::vector<int>* v2i = nullptr;

	auto start = std::chrono::high_resolution_clock::now();

	// iterate over computation mesh with desired pixel size -> (i,j)
#pragma omp parallel for private(j,k,l,I,J,v1,v2,v1i,v2i,index,INDEX)
	for (int i = 0; i < m_areaInfo.height; i++)
	{
		for (j = 0; j < m_areaInfo.width; j++)
		{
			// iterate over normalization mesh -> (I,J): normalization mesh, (k,l): mask
			for (k = 0; k < desiredPixelSize; k++)
			{
				for (l = 0; l < desiredPixelSize; l++)
				{
					// compute pixel indices in normalization mesh
					I = i * desiredPixelSize + k;
					J = j * desiredPixelSize + l;

					if ((I < m_areaInfo.height_n) && (J < m_areaInfo.width_n)) // check if inside normalization mesh
					{
						INDEX = I * m_areaInfo.width_n + J; // compute correct INDEX in normalization mesh
						index = i * m_areaInfo.width   + j; // compute correct index in desired computational mesh
						
						// move vegetation points
						v1 = &m_meshPixelsRedistributed[index].vegetationPts.xCoords;
						v2 = &m_meshPixels[INDEX].vegetationPts.xCoords;
						v1->insert(v1->end(), v2->begin(), v2->end()); // insert v2 to v1
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

						v1 = &m_meshPixelsRedistributed[index].vegetationPts.yCoords;
						v2 = &m_meshPixels[INDEX].vegetationPts.yCoords;
						v1->insert(v1->end(), v2->begin(), v2->end()); // insert v2 to v1
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

						v1 = &m_meshPixelsRedistributed[index].vegetationPts.zCoords;
						v2 = &m_meshPixels[INDEX].vegetationPts.zCoords;
						v1->insert(v1->end(), v2->begin(), v2->end()); // insert v2 to v1
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

						v1i = &m_meshPixelsRedistributed[index].vegetationPts.ptClass;
						v2i = &m_meshPixels[INDEX].vegetationPts.ptClass;
						v1i->insert(v1i->end(), v2i->begin(), v2i->end()); // insert v2 to v1
						v2i->clear(); v2i->shrink_to_fit(); // clear v2 contents

						// move ground points
						v1 = &m_meshPixelsRedistributed[index].groundPts.xCoords;
						v2 = &m_meshPixels[INDEX].groundPts.xCoords;
						v1->insert(v1->end(), v2->begin(), v2->end()); // insert v2 to v1
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

						v1 = &m_meshPixelsRedistributed[index].groundPts.yCoords;
						v2 = &m_meshPixels[INDEX].groundPts.yCoords; // insert v2 to v1
						v1->insert(v1->end(), v2->begin(), v2->end()); // clear v2 contents
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

						v1 = &m_meshPixelsRedistributed[index].groundPts.zCoords;
						v2 = &m_meshPixels[INDEX].groundPts.zCoords; // insert v2 to v1
						v1->insert(v1->end(), v2->begin(), v2->end()); // clear v2 contents
						v2->clear(); v2->shrink_to_fit(); // clear v2 contents

					}
				}
			} // END of mask iteration

		} // proceed to new pixel
	}

	m_meshPixels.clear();
	m_meshPixels.shrink_to_fit();

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	printf("\nPoints redistribution: %.4lf s\n", (double)duration.count() / 1000.0);
}

void DataHandler::computeMetrics()
{
	int nPixels = m_areaInfo.width * m_areaInfo.height;

	std::cout << "Compute metrics started...\n";
	auto start = std::chrono::high_resolution_clock::now();

	m_metrics.resize(nMetrics);
	for (int i = 0; i < nMetrics; i++)
	{
		m_metrics[i] = new double[nPixels] {};
	}

	double meanZ = -1.0, varZ = -1.0, stdZ = -1.0;
	int nAllPts = -1;
	int nGroundPts = -1;
	int nVegetationPts = -1;
	double temp = -1.0;
	std::vector<double>* zValues = nullptr;
	int m = 0;

#pragma omp parallel for private(zValues, nGroundPts, nVegetationPts, nAllPts, meanZ, varZ, stdZ, temp, m)
	for (int i = 0; i < nPixels; i++)
	{
		std::vector<double> emptyZ = { 0.0 };

		// set default metrics value to NAN
		for (m = 0; m < nMetrics; m++)
			m_metrics[m][i] = NO_DATA_VALUE;
		//std::cout << "nans assigned\n";

		zValues = &m_meshPixelsRedistributed[i].vegetationPts.zCoords;

		nGroundPts = m_meshPixelsRedistributed[i].groundPts.zCoords.size();
		nVegetationPts = zValues->size();
		nAllPts = nGroundPts + nVegetationPts;

		//std::cout << "nPoints found\n";

		if (nAllPts == 0)
			continue;

		if (nVegetationPts == 0)
			zValues = &emptyZ;

		// ECOSYSTEM HEIGHT METRICS
		m_metrics[Hmax][i] = *std::max_element(zValues->begin(), zValues->end());

		meanZ = StatFunctions::mean(*zValues);
		m_metrics[Hmean][i] = meanZ;

		m_metrics[Hmedian][i] = StatFunctions::percentile(*zValues, 50.0);

		m_metrics[Hp25][i] = StatFunctions::percentile(*zValues, 25.0);

		m_metrics[Hp75][i] = StatFunctions::percentile(*zValues, 75.0);

		m_metrics[Hp95][i] = StatFunctions::percentile(*zValues, 95.0);

		m_metrics[PPR][i] = static_cast<double>(nGroundPts) / nAllPts;

		temp = StatFunctions::pointsAbove(*zValues, meanZ);
		m_metrics[DAM_z][i] = (nVegetationPts == 0) ? 0.0 : (temp / nVegetationPts);

		//std::cout << "Metrics 1/3 done\n";

		// ECOSYSTEM COVER METRICS
		temp = StatFunctions::pointsBelow(*zValues, 1.0);
		m_metrics[BR_below_1][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;

		temp = StatFunctions::pointsBetween(*zValues, 1.0, 2.0);
		m_metrics[BR_1_2][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsBetween(*zValues, 2.0, 3.0);
		m_metrics[BR_2_3][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsAbove(*zValues, 3.0);
		m_metrics[BR_above_3][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsBetween(*zValues, 3.0, 4.0);
		m_metrics[BR_3_4][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsBetween(*zValues, 4.0, 5.0);
		m_metrics[BR_4_5][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsBelow(*zValues, 5.0);
		m_metrics[BR_below_5][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsBetween(*zValues, 5.0, 20.0);
		m_metrics[BR_5_20][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		temp = StatFunctions::pointsAbove(*zValues, 20.0);
		m_metrics[BR_above_20][i] = (nVegetationPts == 0) ? 0.0 : temp / nVegetationPts;
		
		//std::cout << "Metrics 2/3 done\n";

		// ECOSYSTEM STRUCTURAL COMPLEXITY METRICS
		m_metrics[Hkurt][i] = StatFunctions::kurtosis(*zValues, meanZ);

		m_metrics[Hskew][i] = StatFunctions::skewness(*zValues, meanZ);

		varZ = StatFunctions::variance(*zValues, meanZ);
		m_metrics[Hvar][i] = varZ;

		stdZ = std::sqrt(varZ);
		m_metrics[Hstd][i] = stdZ;

		m_metrics[Coeff_var_z][i] = stdZ / meanZ;

		m_metrics[Shannon][i] = StatFunctions::shannonIndex(*zValues);

		//std::cout << "Metrics 3/3 done\n";

	}

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	printf("done: %.4lf s\n", (double)duration.count() / 1000.0);
}

void DataHandler::exportMetrics(std::string fileName)
{
	//time_t now = time(0);
	//tm* ltm = localtime(&now);
	
	//std::string timeStamp = QString("../../%1_%2_%3_%4-%5-%6_").arg(ltm->tm_year + 1900,4).arg(ltm->tm_mon,2, 10, QLatin1Char('0')).arg(ltm->tm_mday, 2, 10, QLatin1Char('0')).arg(ltm->tm_hour, 2, 10, QLatin1Char('0')).arg(ltm->tm_min, 2, 10, QLatin1Char('0')).arg(ltm->tm_sec, 2, 10, QLatin1Char('0')).toStdString();
	//timeStamp.append(fileName);
	//fileName = timeStamp;

	// fileName = std::string("../../").append(fileName);
	//fileName = std::string("../../_exportedTIFs/").append(fileName);
	// fileName = std::string(".\\") + fileName.append(".tif");
	// fileName += std::string("_h=") + std::to_string(m_areaInfo.desiredPixelSize) + std::string("m");
	fileName = "../../../raster_metrics/"+ fileName + ".tif";

	std::cout << "Exporting metrics..." << std::endl;
	// load drivers
	GDALAllRegister();
	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");

	std::cout << "Driver registered" << std::endl;

	GDALDataset* poDstDS = nullptr;
	char** papszOptions = NULL;
	poDstDS = poDriver->Create(fileName.c_str(), m_areaInfo.width, m_areaInfo.height, m_metrics.size(), GDT_Float64, papszOptions);

	std::cout << "Dataset created" << std::endl;

	// geotransform -> x,y lavy horny vrchol, dlzky v metroch pre pixely, rotacie (tie su 0)
	double adfGeoTransform[6] = { m_areaInfo.xLeft,
								  static_cast<double>(m_areaInfo.desiredPixelSize),
		                          0.0,
		                          m_areaInfo.yLeft,
		                          0.0,
		                         static_cast<double>(- m_areaInfo.desiredPixelSize)};
	poDstDS->SetGeoTransform(adfGeoTransform);

	std::cout << "Set geotransform" << std::endl;

	// Spatial reference system for tif
	OGRSpatialReference oSRS;
	char* pszSRS_WKT = NULL;
	oSRS.importFromProj4("+proj=krovak +lat_0=49.5 +lon_0=24.8333333333333 +alpha=30.2881397527778 +k=0.9999 +x_0=0 +y_0=0 +ellps=bessel +units=m +no_defs +type=crs");
	oSRS.exportToWkt(&pszSRS_WKT);
	poDstDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);

	std::cout << "Set spatial reference" << std::endl;

	// Export bands
	GDALRasterBand* poBand;

	for (int i = 0; i < m_metrics.size(); i++)
	{
		poBand = poDstDS->GetRasterBand(i + 1);
		poBand->SetDescription(m_bandNames[i].c_str());
		poBand->SetColorInterpretation(GCI_GrayIndex);
		poBand->SetNoDataValue(NO_DATA_VALUE);
		poBand->RasterIO(GF_Write, 0, 0, m_areaInfo.width, m_areaInfo.height,
			m_metrics[i], m_areaInfo.width, m_areaInfo.height, GDT_Float64, 0, 0);
	}

	std::cout << "Data written to file" << std::endl;

	// Close file
	GDALClose((GDALDatasetH)poDstDS);

	std::cout << "Export done" << std::endl;
}

bool DataHandler::exportLAS(std::string fileName)
{
	if (fileName.substr(fileName.length() - 4, 4) != ".las")
	{
		fileName.append(".las");
	}

	fileName = std::string("../../").append(fileName);

	std::cout << "Exporting PC to LAS file: '" << fileName << "'\n";

	LASwriteOpener lasWriteOpener;
	lasWriteOpener.set_file_name(fileName.c_str());

	if (!lasWriteOpener.active())
	{
		std::cout << "Las file " << lasWriteOpener.get_file_name() << "not active\n";
		return false;
	}

	LASheader lasHeader;

	lasHeader.x_scale_factor = 0.01;
	lasHeader.y_scale_factor = 0.01;
	lasHeader.z_scale_factor = 0.01;
	lasHeader.x_offset = 0.0;
	lasHeader.y_offset = 0.0;
	lasHeader.z_offset = 0.0;
	lasHeader.point_data_format = 1;
	lasHeader.point_data_record_length = 28;

	LASpoint lasPoint;
	lasPoint.init(&lasHeader, lasHeader.point_data_format, lasHeader.point_data_record_length, 0);

	// open laswriter

	LASwriter* lasWriter = lasWriteOpener.open(&lasHeader);
	if (lasWriter == 0)
	{
		std::cout << "lasWriter NOT opened\n";
		return false;
	}

	VegetationPoints* vPts = nullptr;
	GroundPoints* gPts = nullptr;

	// write points
	for (int i = 0; i < m_meshPixelsRedistributed.size(); i++)
	{
		// write vegetation points
		vPts = &m_meshPixelsRedistributed[i].vegetationPts;

		for (int j = 0; j < vPts->zCoords.size(); j++)
		{
			// populate the point
			lasPoint.set_x(vPts->xCoords[j]);
			lasPoint.set_y(vPts->yCoords[j]);
			lasPoint.set_z(vPts->zCoords[j]);
			lasPoint.set_classification(vPts->ptClass[j]);
			
			// write the point
			lasWriter->write_point(&lasPoint);

			// add it to the inventory
			lasWriter->update_inventory(&lasPoint);
		}

		// write ground points
		gPts = &m_meshPixelsRedistributed[i].groundPts;
		
		for (int j = 0; j < gPts->zCoords.size(); j++)
		{
			// populate the point
			lasPoint.set_x(gPts->xCoords[j]);
			lasPoint.set_y(gPts->yCoords[j]);
			lasPoint.set_z(gPts->zCoords[j]);
			lasPoint.set_classification((U8)GROUND);

			// write the point
			lasWriter->write_point(&lasPoint);

			// add it to the inventory
			lasWriter->update_inventory(&lasPoint);
		}
	}

	// update the header
	lasWriter->update_header(&lasHeader, TRUE);

	// close the writer
	I64 total_bytes = lasWriter->close();

	std::cout << "Total bytes written: " << total_bytes << "\n";

	delete lasWriter;

	return true;
}

bool DataHandler::exportLAZ(std::string fileName)
{
	if (fileName.substr(fileName.length() - 4, 4) != ".laz")
	{
		fileName.append(".laz");
	}

	fileName = std::string("../../").append(fileName);

	std::cout << "Exporting PC to LAZ file: '" << fileName << "'\n";

	// TODO: dorobit podla https://github.com/LAStools/LAStools/blob/master/LASzip/example/laszipdllexample.cpp

	return true;
}

// void DataHandler::exportDTM(std::string fileName)
// {
// 	time_t now = time(0);
// 	tm* ltm = localtime(&now);

// 	std::string timeStamp = QString("../../DTM_%1_%2_%3_%4-%5-%6_").arg(ltm->tm_year + 1900, 4).arg(ltm->tm_mon, 2, 10, QLatin1Char('0')).arg(ltm->tm_mday, 2, 10, QLatin1Char('0')).arg(ltm->tm_hour, 2, 10, QLatin1Char('0')).arg(ltm->tm_min, 2, 10, QLatin1Char('0')).arg(ltm->tm_sec, 2, 10, QLatin1Char('0')).toStdString();
// 	timeStamp.append(fileName);
// 	fileName = timeStamp;

// 	if (fileName.substr(fileName.length() - 4, 4) != ".tif")
// 	{
// 		fileName.append(".tif");
// 	}

// 	std::cout << "Exporting DTM...";
// 	// load drivers
// 	GDALAllRegister();
// 	GDALDriver* poDriver;
// 	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");

// 	GDALDataset* poDstDS = nullptr;
// 	char** papszOptions = NULL;
// 	poDstDS = poDriver->Create(fileName.c_str(), m_areaInfo.width_n, m_areaInfo.height_n, 1, GDT_Float64, papszOptions);

// 	// geotransform -> x,y lavy horny vrchol, dlzky v metroch pre pixely, rotacie (tie su 0)
// 	double adfGeoTransform[6] = { m_areaInfo.xLeft,
// 								  1.0,
// 								  0.0,
// 								  m_areaInfo.yLeft,
// 								  0.0,
// 								  -1.0 };
// 	poDstDS->SetGeoTransform(adfGeoTransform);

// 	// Spatial reference system for tif
// 	OGRSpatialReference oSRS;
// 	char* pszSRS_WKT = NULL;
// 	oSRS.importFromEPSG(8353);
// 	oSRS.exportToWkt(&pszSRS_WKT);
// 	poDstDS->SetProjection(pszSRS_WKT);
// 	CPLFree(pszSRS_WKT);

// 	// Export bands
// 	GDALRasterBand* poBand;

// 	poBand = poDstDS->GetRasterBand(1);
// 	poBand->SetDescription("Digital Terrain Model");
// 	poBand->RasterIO(GF_Write, 0, 0, m_areaInfo.width_n, m_areaInfo.height_n,
// 		m_DTM, m_areaInfo.width_n, m_areaInfo.height_n, GDT_Float64, 0, 0);
	
// 	//poDstDS->SetMetadata();
// 	// Close file
// 	GDALClose((GDALDatasetH)poDstDS);

// 	std::cout << " done\n";
// }

void DataHandler::reset()
{
	int i = 0;
	int j = 0;
//#pragma omp parallel for private(j)
	for (i = 0; i < m_meshPixelsRedistributed.size(); i++)
	{
		m_meshPixelsRedistributed[i].vegetationPts.xCoords.clear();
		m_meshPixelsRedistributed[i].vegetationPts.xCoords.shrink_to_fit();

		m_meshPixelsRedistributed[i].vegetationPts.yCoords.clear();
		m_meshPixelsRedistributed[i].vegetationPts.yCoords.shrink_to_fit();

		m_meshPixelsRedistributed[i].vegetationPts.zCoords.clear();
		m_meshPixelsRedistributed[i].vegetationPts.zCoords.shrink_to_fit();

		m_meshPixelsRedistributed[i].vegetationPts.ptClass.clear();
		m_meshPixelsRedistributed[i].vegetationPts.ptClass.shrink_to_fit();

		m_meshPixelsRedistributed[i].groundPts.xCoords.clear();
		m_meshPixelsRedistributed[i].groundPts.xCoords.shrink_to_fit();

		m_meshPixelsRedistributed[i].groundPts.yCoords.clear();
		m_meshPixelsRedistributed[i].groundPts.yCoords.shrink_to_fit();

		m_meshPixelsRedistributed[i].groundPts.zCoords.clear();
		m_meshPixelsRedistributed[i].groundPts.zCoords.shrink_to_fit();
	}
	std::cout << "vectors cleared\n";

	//m_meshPixelsRedistributed.clear();
	//m_meshPixelsRedistributed.shrink_to_fit();
	std::cout << "meshPixelsRedistributed cleared\n";

	for (int i = 0; i < m_metrics.size(); i++)
	{
		delete[] m_metrics[i];
	}

	m_metrics.clear();
	std::cout << "metrics cleared\n";
}

