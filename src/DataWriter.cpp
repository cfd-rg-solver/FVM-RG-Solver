#include "datawriter.h"
DataWriter::DataWriter(std::string pathName_)
{

	directory = pathName_;
	fs::path dataDir = directory / "data";
	fs::remove_all(dataDir);
	fs::create_directories(dataDir);
}

fs::path DataWriter::createTimeDirectory(double time)
{
	fs::path localDir = directory / "data" / std::to_string(time);
	fs::create_directories(localDir);
	return localDir;
}

void DataWriter::writeData(vector<macroParam> data, double time)
{
	fs::path localDir = createTimeDirectory(time);

	ofstream pressure(localDir / "pressure.txt", std::ios::out);
	ofstream velocity(localDir / "velocity.txt", std::ios::out);
	ofstream velocity_tau(localDir / "velocity_tau.txt", std::ios::out);
	ofstream velocity_normal(localDir / "velocity_normal.txt", std::ios::out);
	ofstream temp(localDir / "temp.txt", std::ios::out);
	ofstream temp_vibr(localDir / "temp_vibr.txt", std::ios::out);
	ofstream density(localDir / "density.txt", std::ios::out);
	ofstream gamma(localDir / "gamma.txt", std::ios::out);
	//    ofstream densityO2(localDir/"densityO2.txt",std::ios::out);
	//    ofstream densityO(localDir/"densityO.txt",std::ios::out);

	//    ofstream e(localDir/"e.txt",std::ios::out);

	pressure << "y" << " " << "p" << endl;
	velocity << "y" << " " << "v" << endl;
	velocity_tau << "y" << " " << "v_t" << endl;
	velocity_normal << "y" << " " << "v_n" << endl;
	temp << "y" << " " << "T" << endl;
	temp_vibr << "y" << " " << "Tvibr" << endl;
	density << "y" << " " << "rho" << endl;
	gamma << "y" << " " << "gamma" << endl;
	//    densityO2<<"y"<<" "<<"rho_O2"<<endl;
	//    densityO<<"y"<<" "<<"rho_O"<<endl;

	for (size_t i = 0; i < data.size(); i++)
	{
		pressure << dh * i << " " << data[i].pressure << endl;
		velocity << dh * i << " " << data[i].velocity << endl;
		velocity_tau << dh * i << " " << data[i].velocity_tau << endl;
		velocity_normal << dh * i << " " << data[i].velocity_normal << endl;
		temp << dh * i << " " << data[i].temp << endl;
		temp_vibr << dh * i << " " << data[i].temp_vibr << endl;
		density << dh * i << " " << data[i].density << endl;
		gamma << dh * i << " " << data[i].gamma << endl;
		//        if(data[i].mixture.NumberOfComponents == 2)
		//        {
		//            densityO2<<dh*i<<" "<<data[i].densityArray[0]<<endl;
		//            densityO<<dh*i<<" "<<data[i].densityArray[1]<<endl;
		//        }
		//        double gamma = 1.4;
		//        double e_ = data[i].pressure/((gamma - 1) * data[i].density);
		//        e<<dh*i<<" "<<e_<<endl;
	}
	pressure.close();
	velocity.close();
	temp.close();
	temp_vibr.close();
	velocity_normal.close();
	velocity_tau.close();
	density.close();
	gamma.close();
}

void DataWriter::setDelta_h(double dh_)
{
	dh = dh_;
}


bool DataReader::read()
{
	if (!fillDataVector(pres, "/pressure.txt"))
		return 0;
	if (!fillDataVector(vel, "/velocity.txt"))
		return 0;
	if (!fillDataVector(temp, "/temp.txt"))
		return 0;
	if (!fillDataVector(temp_vibr, "/temp_vibr.txt"))
		return 0;
	if (!fillDataVector(density, "/density.txt"))
		return 0;
	if (!fillDataVector(vel_tau, "/velocity_tau.txt"))
		return 0;
	if (!fillDataVector(vel_normal, "/velocity_normal.txt"))
		return 0;
	if (!fillDataVector(gamma, "/gamma.txt"))
		return 0;
	return 1;
}

void DataReader::getPoints(vector<macroParam>& points)
{
	points.clear();
	for (size_t i = 0; i < pres.size(); i++)
	{
		macroParam tmp;
		tmp.densityArray = { density[i] };
		tmp.fractionArray = { 1 };
		tmp.density = density[i];
		tmp.pressure = pres[i];
		tmp.velocity_tau = vel_tau[i];
		tmp.velocity_normal = vel_normal[i];
		tmp.velocity = vel[i];
		tmp.temp = temp[i];
		tmp.temp_vibr = temp_vibr[i];
		tmp.gamma = gamma[i];
		tmp.gas = "Ar";
		points.push_back(tmp);
	}
	return;
}

bool DataReader::fillDataVector(vector<double>& data, string dataFileName)
{
	std::ifstream file(pathName + dataFileName); // окрываем файл для чтения
	if (!file.is_open())
	{
		cout << "WARNING: no" << dataFileName << " file to read" << endl;
		return 0;
	}
	double h, value;
	while (file >> h >> value)
	{
		data.push_back(value);
	}
	file.close();     // закрываем файл
	return 1;
}
