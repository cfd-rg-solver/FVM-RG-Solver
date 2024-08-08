#include <iostream>
#include <omp.h>
#include <functional>
// #include "nn.h"
// #include <time.h>
#include "systemofequation.h"
#include "numeric.h"

void SystemOfEquation::setNumberOfCells(size_t cells)
{
	numberOfCells = cells;
}

void SystemOfEquation::setMixture(Mixture mixture_)
{
	mixture = mixture_;
	numberOfComponents = mixture.NumberOfComponents;
	prepareIndex();
}


void SystemOfEquation::setBorderCondition(BorderCondition* border_)
{
	border = border_;
}

void SystemOfEquation::setCoeffSolver(CoeffSolver* coeffSolver_)
{
	coeffSolver = coeffSolver_;
}

void SystemOfEquation::setSolverParams(solverParams solParam_)
{
	solParam = solParam_;
}

void SystemOfEquation::setEqSolver(NonLinearEqSolver* eqSolver_)
{
	eqSolver = eqSolver_;
}

double SystemOfEquation::getDensity(size_t i, size_t component)
{
	if (mixture.NumberOfComponents == 1)
		return getDensity(i);
	else return 0;
}

double SystemOfEquation::getMaxVelocity()
{
	double res = 0;
	for (size_t i = 0; i < numberOfCells; i++)
	{
		double tmp = fabs(getVelocity(i));
		if (res < tmp)
		{
			res = tmp;
		}
	}
	return res;
}

void SystemOfEquation::prepareVectorSizes()
{
	U.resize(systemOrder); // т.е. 1 уравнения неразрывности + 1 уравнения движения + 1 уравнение энергии
	for (size_t i = 0; i < U.size(); i++)
		U[i].resize(numberOfCells);


	F.resize(systemOrder);
	for (size_t i = 0; i < F.size(); i++)
		F[i].resize(numberOfCells);

	Flux.resize(systemOrder);
	for (size_t i = 0; i < Flux.size(); i++)
		Flux[i].resize(numberOfCells - 1); // на одну меньше, т.к. через грани

	R.resize(systemOrder);
	for (size_t i = 0; i < R.size(); i++)
		R[i].resize(numberOfCells);
}


void Couette2::prepareSolving(vector<macroParam>& points)
{

#pragma omp parallel for schedule(static)
	for (auto i = 0; i < numberOfCells; i++)
	{
		U[0][i] = points[i].density;
		for (size_t j = 1; j < numberOfComponents; j++)
			U[j][i] = points[i].densityArray[j];
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[v_normal][i] = points[i].density * points[i].velocity_normal;


		U[energy][i] = points[i].pressure / (2. / 3.) + 0.5 * pow(points[i].velocity, 2) * points[i].density;
		//U[energy][i] = (3.*points[i].pressure)/2. + 0.5*pow(points[i].velocity,2)*points[i].density;
	}
}

void Couette2::prepareIndex()
{
	systemOrder = numberOfComponents + 3;

	v_tau = numberOfComponents;
	v_normal = numberOfComponents + 1;
	energy = numberOfComponents + 2;
}

double Couette2::getPressure(size_t i)
{
	double rho = getDensity(i);
	return (getEnergy(i) - 0.5 * pow(getVelocity(i), 2)) * 2. / 3. * rho;
	//return rho * (UniversalGasConstant/mixture.molarMass()) * getTemp(i);
	//return (solParam.Gamma - 1.)*(U[energy][i]/rho - pow(getVelocity(i),2)/2.);
}

double Couette2::getDensity(size_t i)
{
	return U[0][i];
}

double Couette2::getVelocity(size_t i)
{
	double v_t = U[v_tau][i] / getDensity(i);
	double v_n = U[v_normal][i] / getDensity(i);
	double v = sqrt(pow(v_t, 2) + pow(v_n, 2));
	return v;
}

double Couette2::getVelocityTau(size_t i)
{
	return U[v_tau][i] / getDensity(i);
}

double Couette2::getVelocityNormal(size_t i)
{
	return U[v_normal][i] / getDensity(i);
}

double Couette2::getEnergy(size_t i)
{
	return U[energy][i] / getDensity(i);
}

double Couette2::getTemp(size_t i)
{
	double U = getEnergy(i) - pow(getVelocity(i), 2) / 2.;
	double n_kB = UniversalGasConstant / mixture.molarMass() * getDensity(i);
	return U * 2. / 3. / (n_kB)*getDensity(i);
	//return 0.67*mixture.mass(0) / kB * U;
}

double Couette2::getGamma(size_t i) {
	return solParam.Gamma;
}

void Couette2::updateU(double dh, double dt)
{
#pragma omp parallel for schedule(static)
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (/*R[j][i]*/0 - (Flux[j][i] - Flux[j][i - 1]) / dh) * dt;
		}
	}
}

void Couette2::updateBorderU(vector<macroParam>& points) // ! change calculation of conservative varibles vector in the fictious cells on the basis on BC type
{
	for (int i : {0, (int)(numberOfCells - 1)})
	{
		U[0][i] = points[i].density;
		for (size_t j = 1; j < numberOfComponents; j++)
			U[j][i] = points[i].densityArray[j];
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[v_normal][i] = points[i].density * points[i].velocity_normal;
		U[energy][i] = points[i].pressure / (2. / 3.) + 0.5 * pow(points[i].velocity, 2) * points[i].density;
		//U[energy][i] = (3*points[i].pressure)/2 + 0.5*pow(points[i].velocity,2)*points[i].density;
	}
	return;
}
//void Couette2::computeF(vector<macroParam> &points, double dh)
//{
//    Mixture mixture = points[0].mixture;
//    for(size_t i = 0 ; i < numberOfCells; i++)
//    {
//        // Рассчитываем производные в точке i
//        double dv_tau_dy;
//        double dv_normal_dy;
//        double dT_dy;
//        if(i!=numberOfCells-1)
//        {
//            dT_dy = (points[i+1].temp - points[i].temp) / (dh);
//            dv_tau_dy = (points[i+1].velocity_tau - points[i].velocity_tau) / (dh);
//            dv_normal_dy = (points[i+1].velocity_normal - points[i].velocity_normal) / (dh);
//        }
//        else
//        {
//            dT_dy = -(points[i].temp - points[i-1].temp) / (dh);
//            dv_tau_dy = -(points[i].velocity_tau - points[i-1].velocity_tau) / (dh);
//            dv_normal_dy = -(points[i].velocity_normal - points[i-1].velocity_normal) / (dh);
//        }

//        vector<double> dy_dy(numberOfComponents);

//        //учёт граничных условий
//        if(i == 0 || i == numberOfCells-1)
//            fill(dy_dy.begin(), dy_dy.end(),border->get_dyc_dy());
//        else
//        {
//            for(size_t j = 0 ; j <numberOfComponents; j++)
//            {
//                dy_dy[j] = (points[i+1].fractionArray[j] - points[i].fractionArray[j])/ (dh);
//            }
//        }
//        // Расчет поточных членов
//        // .....
//        // сейчас так:
//        double etta = coeffSolver->shareViscositySimple(points[i]);
//        double lambda = coeffSolver->lambda(points[i]);
//        double bulk = coeffSolver->bulcViscositySimple(mixture,points[i].temp, points[i].density, points[i].pressure);

//        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
//        {
//            if(j!=0)
//                F[j][i] = -points[i].density * mixture.getEffDiff(j) * dy_dy[j];
//            else
//                F[j][i] = points[i].density * points[i].velocity_normal;
//        }
//        F[v_tau][i] = points[i].density * points[i].velocity_tau * points[i].velocity_normal  -etta * dv_tau_dy;
//        F[v_normal][i] = points[i].density *pow(points[i].velocity_normal,2) + points[i].pressure - (bulk + 4/3*etta)* dv_normal_dy;
//        F[energy][i] = 0;
//        for(size_t j = 0 ; j <numberOfComponents; j++)
//        {
//            F[energy][i]+= - points[i].density * mixture.getEffDiff(j)*dy_dy[j] * mixture.getEntalp(i);
//        }
//        F[energy][i] += -lambda*dT_dy - etta*points[i].velocity_tau*dv_tau_dy + (points[i].pressure - (bulk + 4/3*etta)* dv_normal_dy) * points[i].velocity_normal;
//    }
//}


void Couette2::computeF(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные в точке i
		double dv_tau_dy;
		double dv_normal_dy;
		double dT_dy;
		dT_dy = (p2.temp - p0.temp) / denominator;
		dv_tau_dy = (p2.velocity_tau - p0.velocity_tau) / denominator;
		dv_normal_dy = (p2.velocity_normal - p0.velocity_normal) / denominator;

		vector<double> dy_dy(numberOfComponents);

		//учёт граничных условий
		if (i == 0 || i == numberOfCells - 1)
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		else
		{
			for (size_t j = 0; j < numberOfComponents; j++)
			{
				dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j]) / denominator;
			}
		}
		// Расчет поточных членов
		// .....
		// сейчас так:
		double etta = coeffSolver->shareViscositySimple(p1);
		double lambda = coeffSolver->lambda(p1);
		double bulk = coeffSolver->bulcViscositySimple(p1);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			if (j != 0)
				F[j][i] = -p1.density * coeffSolver->effDiffusion(p1, j) * dy_dy[j];
			else
				F[j][i] = p1.density * p1.velocity_normal;
		}
		F[v_tau][i] = p1.density * p1.velocity_tau * p1.velocity_normal - etta * dv_tau_dy;
		F[v_normal][i] = p1.density * pow(p1.velocity_normal, 2) + p1.pressure - (bulk + 4. / 3. * etta) * dv_normal_dy;
		F[energy][i] = 0;
		for (size_t j = 0; j < numberOfComponents; j++)
		{
			F[energy][i] += -p1.density * coeffSolver->effDiffusion(p1, j) * dy_dy[j] * (solParam.Gamma * kB * p1.temp / p1.mixture.components[j].mass)/*mixture.getEntalp(i)*/;
		}
		F[energy][i] += -lambda * dT_dy - etta * p1.velocity_tau * dv_tau_dy + (p1.pressure - (bulk + 4. / 3. * etta) * dv_normal_dy) * p1.velocity_normal;
	}
}

void Couette2Alt::prepareVectorSizes()
{
	SystemOfEquation::prepareVectorSizes();

	Fv.resize(systemOrder);
	for (size_t i = 0; i < Fv.size(); i++)
		Fv[i].resize(numberOfCells);
}

void Couette2Alt::updateU(double dh, double dt)
{
#pragma omp parallel for schedule(static)
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (/*R[j][i]*/0 - (Flux[j][i] - Flux[j][i - 1]) / dh - (Fv[j][i] - Fv[j][i - 1]) / (/*2.**/dh)) * dt;
		}
	}
}

void Couette2Alt::computeF(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные в точке i


		vector<double> dy_dy(numberOfComponents);

		//учёт граничных условий
		if (i == 0 || i == numberOfCells - 1)
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		else
		{
			for (size_t j = 0; j < numberOfComponents; j++)
			{
				dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j]) / denominator;
			}
		}

		double dT_dy;
		dT_dy = (p2.temp - p0.temp) / denominator;
		//        double lambda = coeffSolver->lambda(p1);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			if (j != 0)
				F[j][i] = 0;
			else
				F[j][i] = p1.density * p1.velocity_normal;
		}
		F[v_tau][i] = p1.density * p1.velocity_tau * p1.velocity_normal;
		F[v_normal][i] = p1.density * pow(p1.velocity_normal, 2) + p1.pressure;
		F[energy][i] =  /*-lambda*dT_dy*/ +p1.pressure * p1.velocity_normal + p1.density * p1.velocity_normal * getEnergy(i);
	}
}

void Couette2Alt::computeFv(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные в точке i
		double dv_tau_dy;
		double dv_normal_dy;
		double dT_dy;
		dT_dy = (p2.temp - p0.temp) / denominator;
		dv_tau_dy = (p2.velocity_tau - p0.velocity_tau) / denominator;
		dv_normal_dy = (p2.velocity_normal - p0.velocity_normal) / denominator;

		vector<double> dy_dy(numberOfComponents);

		//учёт граничных условий
		if (i == 0 || i == numberOfCells - 1)
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		else
		{
			for (size_t j = 0; j < numberOfComponents; j++)
			{
				dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j]) / denominator;
			}
		}
		// Расчет поточных членов
		// .....
		// сейчас так:
		double lambda = coeffSolver->lambda(p1);
		double etta = coeffSolver->shareViscositySimple(p1);
		double bulk = coeffSolver->bulcViscositySimple(p1);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			if (j != 0)
				Fv[j][i] = -p1.density * coeffSolver->effDiffusion(p1, j) * dy_dy[j];
			else
				Fv[j][i] = 0;
		}
		Fv[v_tau][i] = -etta * dv_tau_dy;
		Fv[v_normal][i] = -(bulk + 4. / 3. * etta) * dv_normal_dy;
		Fv[energy][i] = 0;
		for (size_t j = 0; j < numberOfComponents; j++)
		{
			double effDiffCoeff = coeffSolver->effDiffusion(p1, j);
			Fv[energy][i] += -p1.density * effDiffCoeff * dy_dy[j] * (solParam.Gamma * kB * p1.temp / p1.mixture.components[j].mass)/*mixture.getEntalp(i)*/;
		}
		Fv[energy][i] = -lambda * dT_dy - etta * p1.velocity_tau * dv_tau_dy - (bulk + 4. / 3. * etta) * dv_normal_dy * p1.velocity_normal;
	}
}

void Couette2AltBinary::prepareSolving(vector<macroParam>& points)
{
	temperature.resize(numberOfCells);
#pragma omp parallel for schedule(static)
	for (auto i = 0; i < numberOfCells; i++)
	{
		U[0][i] = points[i].density;
		for (size_t j = 1; j < numberOfComponents; j++)
			U[j][i] = points[i].densityArray[j];
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[v_normal][i] = points[i].density * points[i].velocity_normal;

		U[energy][i] = energyCalculator->calcEnergy(points[i]);
		temperature[i] = points[i].temp;
	}
}

double Couette2AltBinary::getPressure(size_t i)
{
	double sum = 0;
	for (size_t j = 0; j < numberOfComponents; j++)
	{
		double y_c = getDensity(i, j) / getDensity(i);
		double M_c = mixture.molarMass(j);
		sum += y_c / M_c;
	}
	double M = 1 / sum;
	double pres = getDensity(i) * UniversalGasConstant * getTemp(i) / M;
	return pres;
}

double Couette2AltBinary::getDensity(size_t i)
{
	return U[0][i];
}

double Couette2AltBinary::getDensity(size_t i, size_t component)
{
	if (component == 0)
	{
		double sum = 0;
		for (size_t j = 1; j < numberOfComponents; j++)
		{
			sum += U[j][i];
		}
		return U[component][i] - sum;
	}
	else
		return U[component][i];
}

double Couette2AltBinary::getTemp(size_t i)
{
	return temperature[i];
}

double Couette2AltBinary::getGamma(size_t i) {
	return solParam.Gamma;
}

void Couette2AltBinary::updateU(double dh, double dt)
{
#pragma omp parallel for schedule(static)
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (/*R[j][i]*/0 - (Flux[j][i] - Flux[j][i - 1]) / dh - (Fv[j][i] - Fv[j][i - 1]) / (/*2.**/dh)) * dt;
		}
	}
	calcAndRememberTemp();
}

void Couette2AltBinary::updateBorderU(vector<macroParam>& points) // ! change calculation of conservative varibles vector in the fictious cells on the basis on BC type
{
	for (int i : {0, (int)(numberOfCells - 1)})
	{
		U[0][i] = points[i].density;
		for (size_t j = 1; j < numberOfComponents; j++)
			U[j][i] = points[i].densityArray[j];
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[v_normal][i] = points[i].density * points[i].velocity_normal;

		U[energy][i] = energyCalculator->calcEnergy(points[i]);
	}
	return;
}

void Couette2AltBinary::computeF(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные в точке i


		vector<double> dy_dy(numberOfComponents);

		//учёт граничных условий
		if (i == 0 || i == numberOfCells - 1)
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		else
		{
			for (size_t j = 0; j < numberOfComponents; j++)
			{
				dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j]) / denominator;
			}
		}

		double dT_dy;
		dT_dy = (p2.temp - p0.temp) / denominator;
		//        double lambda = coeffSolver->lambda(p1);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			if (j != 0)
				F[j][i] = p1.density * p1.fractionArray[j] * p1.velocity_normal;
			else
				F[j][i] = p1.density * p1.velocity_normal;
		}
		F[v_tau][i] = p1.density * p1.velocity_tau * p1.velocity_normal;
		F[v_normal][i] = p1.density * pow(p1.velocity_normal, 2) + p1.pressure;
		F[energy][i] = p1.pressure * p1.velocity_normal + p1.density * p1.velocity_normal * getEnergy(i);
	}
}

void Couette2AltBinary::computeFv(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные в точке i
		double dv_tau_dy;
		double dv_normal_dy;
		double dT_dy;
		dT_dy = (p2.temp - p0.temp) / denominator;
		dv_tau_dy = (p2.velocity_tau - p0.velocity_tau) / denominator;
		dv_normal_dy = (p2.velocity_normal - p0.velocity_normal) / denominator;

		vector<double> dy_dy(numberOfComponents);

		//учёт граничных условий
		if (i == 0 || i == numberOfCells - 1)
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		else
		{
			for (size_t j = 0; j < numberOfComponents; j++)
			{
				dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j]) / denominator;
			}
		}
		// Расчет поточных членов
		// .....
		// сейчас так:
		double lambda = coeffSolver->lambda(p1);
		double etta = coeffSolver->shareViscositySimple(p1);
		double bulk = coeffSolver->bulcViscositySimple(p1);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			if (j != 0)
				Fv[j][i] = -p1.density * coeffSolver->effDiffusion(p1, j) * dy_dy[j];
			else
				Fv[j][i] = 0;
		}
		Fv[v_tau][i] = -etta * dv_tau_dy;
		Fv[v_normal][i] = -(bulk + 4. / 3. * etta) * dv_normal_dy;
		Fv[energy][i] = 0;
		for (size_t j = 0; j < numberOfComponents; j++)
		{
			double effDiffCoeff = coeffSolver->effDiffusion(p1, j);
			Fv[energy][i] += -p1.density * effDiffCoeff * dy_dy[j] * energyCalculator->getEntalp(points[i], j);
		}
		Fv[energy][i] = -lambda * dT_dy - etta * p1.velocity_tau * dv_tau_dy - (bulk + 4. / 3. * etta) * dv_normal_dy * p1.velocity_normal;
	}
}
void Couette2AltBinary::calcAndRememberTemp()
{
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0(mixture);
		p0.temp = temperature[i];
		p0.densityArray.resize(numberOfComponents);
		p0.fractionArray.resize(numberOfComponents);
		for (size_t j = 0; j < numberOfComponents; j++)
		{

			p0.densityArray[j] = getDensity(i, j);
			p0.fractionArray[j] = p0.densityArray[j] / getDensity(i);
			p0.velocity = getVelocity(i);
		}
		p0.density = getDensity(i);

		temperature[i] = eqSolver->solveEq(energyCalculator, p0, U[energy][i]);
	}
	return;
}

void Soda::prepareSolving(vector<macroParam>& points)
{
	for (auto i = 0; i < numberOfCells; i++)
	{
		U[0][i] = points[i].density;
		U[v_tau][i] = points[i].density * points[i].velocity;
		double e = points[i].pressure / ((gamma - 1) * points[i].density);
		U[energy][i] = points[i].density * 0.5 * pow(points[i].velocity, 2) + points[i].density * e;
	}
}

void Soda::prepareIndex()
{
	systemOrder = numberOfComponents + 2;

	v_tau = numberOfComponents;
	v_normal = numberOfComponents;
	energy = numberOfComponents + 1;
}

double Soda::getPressure(size_t i)
{
	double rho = getDensity(i);
	return ((getEnergy(i) * rho - 0.5 * rho * pow(getVelocity(i), 2)) * (gamma - 1));

}

double Soda::getDensity(size_t i)
{
	return U[0][i];
}

double Soda::getVelocity(size_t i)
{
	return U[v_tau][i] / getDensity(i);
}

double Soda::getVelocityTau(size_t i)
{
	return U[v_tau][i] / getDensity(i);
}

double Soda::getVelocityNormal(size_t i)
{
	return 0;
}

double Soda::getSoundSpeed(size_t i)
{
	return sqrt(gamma * getPressure(i) / getDensity(i));
}

double Soda::getGamma(size_t i) {
	return solParam.Gamma;
}

double Soda::getEnergy(size_t i)
{
	return U[energy][i] / getDensity(i); // important
}

void Soda::updateU(double dh, double dt)
{
	int last = numberOfCells - 1;
	for (int j = 0; j < systemOrder; j++)
	{
		U[j][0] += 0;
		U[j][last] += 0;
	}
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (/*R[j][i]*/0 - (Flux[j][i] - Flux[j][i - 1]) / dh) * dt;
		}
	}

}

void Soda::computeF(vector<macroParam>& points, double dh)
{
	for (size_t i = 0; i < numberOfCells; i++)
	{
		F[0][i] = points[i].density * points[i].velocity;
		F[v_tau][i] = points[i].density * pow(points[i].velocity, 2) + points[i].pressure;
		F[energy][i] = points[i].velocity * (U[energy][i] + points[i].pressure);
	}
}

//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

void Shockwave1::prepareSolving(vector<macroParam>& points)
{
#pragma omp parallel for schedule(static)
	for (auto i = 0; i < numberOfCells; i++)
	{
		U[0][i] = points[i].density;
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[energy][i] = points[i].pressure / (2. / 3.) + 0.5 * pow(points[i].velocity_tau, 2) * points[i].density;
	}
}

void Shockwave1::prepareVectorSizes()
{
	SystemOfEquation::prepareVectorSizes();

	Fv.resize(systemOrder);
	for (size_t i = 0; i < Fv.size(); i++)
		Fv[i].resize(numberOfCells);
}

void Shockwave1::prepareIndex()
{
	systemOrder = numberOfComponents + 2; // однокомпонентная постановка
	v_tau = numberOfComponents;
	energy = numberOfComponents + 1;
}

double Shockwave1::getDensity(size_t i)
{
	return U[0][i];
}

double Shockwave1::getVelocity(size_t i)
{
	double v = U[v_tau][i] / getDensity(i);
	return v;
}

double Shockwave1::getVelocityTau(size_t i)
{
	double v = U[v_tau][i] / getDensity(i);
	return v;
}


double Shockwave1::getTemp(size_t i)
{
	// double E_energy = U[energy][i] / getDensity(i); // полная энергия E
	// double U_energy = E_energy - 0.5 * pow(getVelocity(i), 2); // внутренняя энергия U

	// однокомпонентная - U_energy = 3*n*k*T/(2*density) + k*T/mass + <e_i>_vibr/mass + e_c/mass
	// многокомпонентная - U_energy = 3*n*k*T/(2*density) +
	// + sum([k*T/mass[i])*fractionArray[i] for i in range(numberOfComponents)]) +
	// + sum([fractionArray[i]*<e_i>_vibr/mass[i] for i in range(numberOfComponents)]) +
	// + sum([fractionArray[i]*e_c/mass[i] for i in range(numberOfComponents)])

	// double n_kB = UniversalGasConstant / mixture.molarMass() * getDensity(i);
	// double T = U_energy * 2./3. / (n_kB) * getDensity(i);

	double U = getEnergy(i) - pow(getVelocity(i), 2) / 2.;
	double n_kB = UniversalGasConstant / mixture.molarMass() * getDensity(i);
	return U * 2. / 3. / (n_kB)*getDensity(i);

	// return T;
}

double Shockwave1::getEnergy(size_t i) {
	return U[energy][i] / getDensity(i); // из однотемпературной модели
}

double Shockwave1::getPressure(size_t i) {
	// return UniversalGasConstant * getTemp(i) * getDensity(i) / mixture.molarMass();
	double rho = getDensity(i);
	return (getEnergy(i) - 0.5 * pow(getVelocity(i), 2)) * 2. / 3. * rho;
}

void Shockwave1::updateU(double dh, double dt)
{
#pragma omp parallel for schedule(static)
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (0 - (Flux[j][i] - Flux[j][i - 1]) / dh - (Fv[j][i + 1] - Fv[j][i - 1]) / (2 * dh)) * dt;
		}
	}
}

void Shockwave1::updateBorderU(vector<macroParam>& points) {
	for (int i : {0, (int)(numberOfCells - 1)})
	{
		U[0][i] = points[i].density;
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[energy][i] = points[i].pressure / (2. / 3.) + 0.5 * pow(points[i].velocity_tau, 2) * points[i].density;
		//or energyCalculator->calcEnergy(points[i]);
	}
	return;
}

double Shockwave1::getGamma(size_t i) {
	return solParam.Gamma;
}

void Shockwave1::computeF(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
	macroParam p1;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		// Переобозначаем величины в ячейках (не в фиктивных):
		p1 = points[i];
		// 1-е уравнение (однокомпонентная постановка) в векторе F с консервативными составляющими:
		F[0][i] = p1.density * p1.velocity_tau;
		// Последние 2 уравнения в векторе F с консервативными составляющими:
		F[v_tau][i] = p1.density * pow(p1.velocity_tau, 2) + p1.pressure;
		F[energy][i] = p1.density * p1.velocity_tau * getEnergy(i) + p1.pressure * p1.velocity_tau;
	}
}

void Shockwave1::computeFv(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		// Переобозначаем величины в ячейках (не в фиктивных):
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}

		// Рассчитываем производные:
		double dv_tau_dy = (p2.velocity_tau - p0.velocity_tau) / denominator;
		double dT_dy = (p2.temp - p0.temp) / denominator;
		vector<double> dy_dy(1);
		// Учёт граничных условий:
		if (i == 0 || i == numberOfCells - 1) {
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		}
		else {
			dy_dy[0] = (p2.fractionArray[0] - p0.fractionArray[0]) / denominator;
		}

		// Расчет потоковых членов:
		double lambda = coeffSolver->lambda(p1);
		double etta = coeffSolver->shareViscositySimple(p1);
		double bulk = coeffSolver->bulcViscositySimple(p1);

		// 1-е уравнение (однокомпонентная постановка) в векторе F с вязкими составляющими:
		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			Fv[j][i] = 0;
		}
		// Последние 2 уравнения в векторе F с вязкими составляющими:
		Fv[v_tau][i] = -(bulk + 4. / 3. * etta) * dv_tau_dy;
		Fv[energy][i] = -lambda * dT_dy - (bulk + 4. / 3. * etta) * dv_tau_dy * p1.velocity_tau;
	}
}

//////////////////////////////////////////////////////////////

void Shockwave2::prepareIndex()
{
	systemOrder = numberOfComponents + 2; // one-component gas, 3 equations
	v_tau = numberOfComponents;
	energy = numberOfComponents + 1;
}

void Shockwave2::prepareSolving(vector<macroParam>& points)
{
	temperature.resize(numberOfCells);
#pragma omp parallel for schedule(static)
	for (auto i = 0; i < numberOfCells; i++)
	{
		U[0][i] = points[i].density;
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[energy][i] = energyCalculator->calcEnergy(points[i]);
		temperature[i] = points[i].temp;
	}
}

void Shockwave2::prepareVectorSizes()
{
	SystemOfEquation::prepareVectorSizes();

	Fv.resize(systemOrder);
	for (size_t i = 0; i < Fv.size(); i++)
		Fv[i].resize(numberOfCells);
}

double Shockwave2::getDensity(size_t i)
{
	return U[0][i];
}

double Shockwave2::getVelocity(size_t i)
{
	double v = U[v_tau][i] / getDensity(i);
	return v;
}

double Shockwave2::getVelocityTau(size_t i)
{
	double v = U[v_tau][i] / getDensity(i);
	return v;
}

double Shockwave2::getEnergy(size_t i) {
	return U[energy][i] / getDensity(i); // because one-temp model
}

double Shockwave2::getTemp(size_t i)
{
	return temperature[i];
}

double Shockwave2::getGamma(size_t i)
{
	macroParam p0(mixture);
	p0.temp = temperature[i];
	p0.velocity = getVelocity(i);
	p0.density = getDensity(i);
	double gamma = energyCalculator->getGamma(p0);
	return gamma;
}

double Shockwave2::getPressure(size_t i)
{
	return getDensity(i) * UniversalGasConstant * getTemp(i) / mixture.molarMass();
}

double Shockwave2::getNormalStress(size_t i)
{ 
	// Normal stress - P = p - 4/3*eta*div(v) - zeta*div(v)
    return Fv[v_tau][i] + getPressure(i);
}

double Shockwave2::getHeatFlux(size_t i)
{
	// Heat flux - q = -lambda*grad(T)
    return Fv[energy][i] - Fv[v_tau][i] * getVelocity(i);
}

void Shockwave2::updateU(double dh, double dt)
{
#pragma omp parallel for schedule(static)
	for (auto i = 1; i < numberOfCells - 1; i++)
	{
		for (int j = 0; j < systemOrder; j++)
		{
			U[j][i] += (0 - (Flux[j][i] - Flux[j][i - 1]) / dh - (Fv[j][i + 1] - Fv[j][i - 1]) / (2 * dh)) * dt;
		}
	}
	calcAndRememberTemp();
}

void Shockwave2::updateBorderU(vector<macroParam>& points) {
	for (int i : {0, (int)(numberOfCells - 1)})
	{
		U[0][i] = points[i].density;
		U[v_tau][i] = points[i].density * points[i].velocity_tau;
		U[energy][i] = energyCalculator->calcEnergy(points[i]); // multi-modes energy calc
	}
}

void Shockwave2::computeF(vector<macroParam>& points, double dh)
{
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p1;
		p1 = points[i];
		F[0][i] = p1.density * p1.velocity_tau;
		F[v_tau][i] = p1.density * pow(p1.velocity_tau, 2) + p1.pressure;
		F[energy][i] = p1.density * p1.velocity_tau * getEnergy(i) + p1.pressure * p1.velocity_tau;
	}
}

void Shockwave2::computeFv(vector<macroParam>& points, double dh)
{
	// clock_t start = clock();
	Mixture mixture = points[0].mixture;
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		// Переобозначаем величины в ячейках (не в фиктивных):
		macroParam p0, p1, p2;
		double denominator = 1.;
		if (i != 0 && i != numberOfCells - 1)
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = 2. * dh;
		}
		else if (i == 0)
		{
			p0 = points[i];
			p1 = points[i];
			p2 = points[i + 1];
			denominator = dh;
		}
		else if (i == (numberOfCells - 1))
		{
			p0 = points[i - 1];
			p1 = points[i];
			p2 = points[i];
			denominator = dh;
		}
		// Рассчитываем производные:
		double dv_tau_dy = (p2.velocity_tau - p0.velocity_tau) / denominator;
		double dT_dy = (p2.temp - p0.temp) / denominator;
		vector<double> dy_dy(1);
		// Учёт граничных условий:
		if (i == 0 || i == numberOfCells - 1) {
			fill(dy_dy.begin(), dy_dy.end(), border->get_dyc_dy());
		}
		else {
			dy_dy[0] = (p2.fractionArray[0] - p0.fractionArray[0]) / denominator;
		}
		/*
		double inputs[1][2] = {
			(p1.pressure - P_MIN) / (P_MAX - P_MIN),
			(p1.temp - T_MIN) / (T_MAX - T_MIN)
		};

		double layersout[2][1] = { 0., 0. };
		for (int i = 0; i < 50; i++) {
			layersout[0][0] +=
				tanh(
					(inputs[0][0] * zllayers0weight[0][i])
					+ zllayers0bias[0][i]
				) * zllayers2weight[0][i];
			layersout[1][0] =
				tanh(
					(inputs[0][1] * zllayers0weight[1][i])
					+ zllayers0bias[0][i]
				) * zllayers2weight[1][i];
		}
		// Расчет потоковых членов (NN):
		double lambda = pow(10, -(layersout[1][0] + zllayers2bias[1][0]));
		double etta = coeffSolver->shareViscositySimple(p1);
		double bulk = pow(10, -(layersout[0][0]) + zllayers2bias[0][0]);
		*/
		// Расчет потоковых членов (Teor):
		double bulk = coeffSolver->bulkViscosityMultiAtom(p1);
		double etta = coeffSolver->shareViscositySimple(p1);
		double lambda = coeffSolver->lambdaMultiAtom(p1);
		// 1-е уравнение (однокомпонентная постановка) в векторе F с вязкими составляющими:
		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			Fv[j][i] = 0;
		}
		// Последние 2 уравнения в векторе F с вязкими составляющими:
		Fv[v_tau][i] = -(bulk + 4. / 3. * etta) * dv_tau_dy;
		Fv[energy][i] = -lambda * dT_dy - (bulk + 4. / 3. * etta) * dv_tau_dy * p1.velocity_tau;
	}
	// clock_t end = clock();
	// double seconds = (double)(end - start) / CLOCKS_PER_SEC;
	// printf("Time of Fv calculation: %f seconds\n", seconds);
}

void Shockwave2::calcAndRememberTemp()
{
#pragma omp parallel for schedule(static)
	for (int i = 0; i < numberOfCells; i++)
	{
		macroParam p0(mixture);
		p0.temp = temperature[i];
		p0.temp_vibr = temperature[i]; // one-temp model
		p0.densityArray.resize(1);
		p0.fractionArray.resize(1);
		p0.densityArray[0] = getDensity(i);
		p0.pressure = getPressure(i);
		p0.fractionArray[0] = p0.densityArray[0] / getDensity(i);
		p0.velocity = getVelocity(i);
		p0.density = getDensity(i);
		p0.gamma = getGamma(i);
		temperature[i] = eqSolver->solveEq(energyCalculator, p0, U[energy][i]);
	}
	return;
}
