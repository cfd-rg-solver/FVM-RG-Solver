#pragma once
#include "global.h"
#include "macroparam.h"
#include "energycalc.h"

struct BorderCondition
{
	virtual void updatePoints(vector<macroParam>& points) = 0;
	virtual double get_dyc_dy() { return 0; }; //затычка для более серьёзных условий
	virtual void setEnergyCalculator(EnergyCalc* energyCalculator_) { energyCalculator = energyCalculator_; };

protected:
	EnergyCalc* energyCalculator;

};

struct BorderConditionCouette : public BorderCondition
{
	void updatePoints(vector<macroParam>& points);
	void setWallParameters(double up_velocity_, double down_velocity_, double up_temp_, double down_temp_)
	{
		up_velocity = up_velocity_; down_velocity = down_velocity_; up_temp = up_temp_; down_temp = down_temp_;
	};
	double get_dyc_dy() { return 0; };
protected:
	double up_velocity, down_velocity, up_temp, down_temp;
};


struct BorderConditionSoda : public BorderCondition
{
	void updatePoints(vector<macroParam>& points);
	void setBorderParameters(double left_velocity_, double left_density_, double left_pressure_,
		double right_velocity_, double right_density_, double right_pressure_)
	{
		left_velocity = left_velocity_; left_density = left_density_; left_pressure = left_pressure_;
		right_velocity = right_velocity_; right_density = right_density_; right_pressure = right_pressure_;
	};
	double get_dyc_dy() { return 0; };
protected:
	double left_velocity, left_density, left_pressure, right_velocity, right_density, right_pressure; // these values wasn't used anywhere
};

struct BorderConditionPersonal : public BorderConditionCouette
{
	void updatePoints(vector<macroParam>& points);
	double get_dyc_dy() { return 0; };
};


struct BorderConditionShockwave : public BorderCondition // don't see the use of this class, setBorderParameters should set initial gas distribution
{
	void updatePoints(vector<macroParam>& points);
	void setBorderParameters(double left_velocity_, double left_density_, double left_temp_,
		double right_velocity_, double right_density_, double right_temp_)
	{
		left_velocity = left_velocity_; left_density = left_density_; left_temp = left_temp_;
		right_velocity = right_velocity_; right_density = right_density_; right_temp = right_temp_;
	};
	double get_dyc_dy() { return 0; };
protected:
	double left_velocity, left_density, left_temp, right_velocity, right_density, right_temp; // these values wasn't used anywhere
};