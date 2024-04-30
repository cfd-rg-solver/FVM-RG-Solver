#pragma once

#include "macroparam.h"
#include "bordercondition.h"
#include "energycalc.h"

struct StartCondition
{
	virtual void setStartDistribution(vector<macroParam>& points) = 0;
	virtual void setEnergyCalculator(EnergyCalc* energyCalculator_) { energyCalculator = energyCalculator_; };

protected:
	EnergyCalc* energyCalculator;
};

struct UniformDistribution : public StartCondition
{
	void setStartDistribution(vector<macroParam>& points);
	void setDistributionParameter(vector<macroParam>& example_)
	{
		exampleVec = example_; newSolving = 0;
	};
	void setDistributionParameter(macroParam example_)
	{
		example = example_; newSolving = 1;
	};
protected:
	bool newSolving;
	macroParam example;
	vector<macroParam> exampleVec;
};

struct UniformDistributionBorder : public UniformDistribution
{
	void setStartDistribution(vector<macroParam>& points);
	void setBorderCondition(BorderCondition* borderCondition_)
	{
		borderCondition = borderCondition_;
	};
protected:
	BorderCondition* borderCondition;
};

struct UniformDistributionBorderPersonal : public UniformDistributionBorder
{
	void setStartDistribution(vector<macroParam>& points)
	{
		UniformDistributionBorder::setStartDistribution(points); points[0].velocity_normal = start_velocity_normal;
	};
	void setNormalVelocity(double start_velocity_normal_)
	{
		start_velocity_normal = start_velocity_normal_;
	};
protected:
	BorderCondition* borderCondition;
	double start_velocity_normal;
};

struct GapDistribution : public StartCondition
{
	void setStartDistribution(vector<macroParam>& points);
	void setDistributionParameter(macroParam left_, macroParam right_) { left = left_; right = right_; };
protected:
	macroParam left, right;
};

struct ShockwaveGapDistribution : public GapDistribution
{
	void setStartDistribution(vector<macroParam>& points);
	void setDistributionParameter(macroParam left_, macroParam right_) { left = left_; right = right_; };
protected:
	macroParam left, right;
};