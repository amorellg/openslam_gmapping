#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>

namespace GMapping{

class RangeSensor: public Sensor{
	friend class Configuration;
	friend class CarmenConfiguration;
	friend class CarmenWrapper;
	public:
		struct Beam{
			OrientedPoint pose;	//pose relative to the center of the sensor
			double span;	//spam=0 indicates a line-like beam
			double maxRange;	//maximum range of the sensor
			double s,c;		//sinus and cosinus of the beam (optimization);
		};
		RangeSensor(std::string name);
		RangeSensor(std::string name, unsigned int beams, double res, int orientation=1, const OrientedPoint& position=OrientedPoint(0,0,0), double span=0, double maxrange=89.0);
    ~RangeSensor();
		inline const std::vector<Beam>& beams() const {return m_beams;}
		inline const int getOrientation() const {return m_orientation;}
		inline const double getAngleInc() const {return m_anginc;}
		inline const double getMaxRange() const {return m_maxRange;}
		inline double * getAngles() const {return m_angles;}
		inline std::vector<Beam>& beams() {return m_beams;}
		inline OrientedPoint getPose() const {return m_pose;}
		void updateBeamsLookup();
		bool newFormat;
	protected:
    double m_maxRange;
		OrientedPoint m_pose;
		std::vector<Beam> m_beams;
    int m_orientation;
    double m_anginc;
    double* m_angles;
};

};

#endif
