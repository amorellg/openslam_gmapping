#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#undef PLAIN_READINGS
#define PLAIN_READINGS

#include "icp.h"
#include "smmap.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>
#include <gmapping/sensor/sensor_range/rangereading.h>
#define LASER_MAXBEAMS 2048

namespace GMapping {

class ScanMatcher{
	public:
		typedef Covariance3 CovarianceMatrix;

		ScanMatcher();
		~ScanMatcher();
		double icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
#ifdef PLAIN_READINGS
		double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		double optimize(OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		double registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
#endif
    //new version which uses RangeReading directly
    double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const RangeReading& reading) const;
    //new version which uses RangeReading directly
    double registerScan(ScanMatcherMap& map, const OrientedPoint& p, const RangeReading& reading);

		void setLaserParameters
			(unsigned int beams, double* angles, const OrientedPoint& lpose);
		void setMatchingParameters
			(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations, double likelihoodSigma=1, unsigned int likelihoodSkip=0 );
		void invalidateActiveArea();
#ifdef PLAIN_READINGS
		void computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
#endif
    //new version which uses RangeReading directly
    void computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const RangeReading& reading);
		inline double icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
#ifdef PLAIN_READINGS
		inline double score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		inline unsigned int likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
#endif
    //new version which uses RangeReading directly
    inline unsigned int likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p,  const RangeReading& reading) const;
    //new version which uses RangeReading directly
    inline double score(const ScanMatcherMap& map, const OrientedPoint& p, const RangeReading& reading) const;

		double likelihood(double& lmax, OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		double likelihood(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, Gaussian3& odometry, const double* readings, double gain=180.);
		inline const double* laserAngles() const { return m_laserAngles; }
		inline unsigned int laserBeams() const { return m_laserBeams; }

		static const double nullLikelihood;
	protected:
		//state of the matcher
		bool m_activeAreaComputed;

		/**laser parameters*/
		unsigned int m_laserBeams;
		double       m_laserAngles[LASER_MAXBEAMS];
		//OrientedPoint m_laserPose;
		PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)
		PARAM_SET_GET(double, laserMaxRange, protected, public, public)
		/**scan_matcher parameters*/
		PARAM_SET_GET(double, usableRange, protected, public, public)
		PARAM_SET_GET(double, gaussianSigma, protected, public, public)
		PARAM_SET_GET(double, likelihoodSigma, protected, public, public)
		PARAM_SET_GET(int,    kernelSize, protected, public, public)
		PARAM_SET_GET(double, optAngularDelta, protected, public, public)
		PARAM_SET_GET(double, optLinearDelta, protected, public, public)
		PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)
		PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)
		PARAM_SET_GET(double, llsamplerange, protected, public, public)
		PARAM_SET_GET(double, llsamplestep, protected, public, public)
		PARAM_SET_GET(double, lasamplerange, protected, public, public)
		PARAM_SET_GET(double, lasamplestep, protected, public, public)
		PARAM_SET_GET(bool, generateMap, protected, public, public)
		PARAM_SET_GET(double, enlargeStep, protected, public, public)
		PARAM_SET_GET(double, fullnessThreshold, protected, public, public)
		PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)
		PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)
		PARAM_SET_GET(double, freeCellRatio, protected, public, public)
		PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)

		// allocate this large array only once
		IntPoint* m_linePoints;
};

inline double ScanMatcher::icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	std::list<PointPair> pairs;

	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange||*r==0.0) continue;
		if (skip) continue;
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		bool found=false;
		Point bestMu(0.,0.);
		Point bestCell(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
					Point mu=phit-cell.mean();
					if (!found){
						bestMu=mu;
						bestCell=cell.mean();
						found=true;
					}else
						if((mu*mu)<(bestMu*bestMu)){
							bestMu=mu;
							bestCell=cell.mean();
						}

				}
			//}
		}
		if (found){
			pairs.push_back(std::make_pair(phit, bestCell));
			//std::cerr << "(" << phit.x-bestCell.x << "," << phit.y-bestCell.y << ") ";
		}
		//std::cerr << std::endl;
	}

	OrientedPoint result(0,0,0);
	//double icpError=icpNonlinearStep(result,pairs);
	std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta << std::endl;
	pret.x=p.x+result.x;
	pret.y=p.y+result.y;
	pret.theta=p.theta+result.theta;
	pret.theta=atan2(sin(pret.theta), cos(pret.theta));
	return score(map, p, readings);
}

#ifdef PLAIN_READINGS
inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	double s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (skip||*r>m_usableRange||*r==0.0) continue;
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
					Point mu=phit-cell.mean();
					if (!found){
						bestMu=mu;
						found=true;
					}else
						bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
				}
			//}
		}
		if (found)
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
	}
	return s;
}

inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	using namespace std;
	l=0;
	s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	double noHit=nullLikelihood/(m_likelihoodSigma);
	unsigned int skip=0;
	unsigned int c=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange) continue;
		if (skip) continue;
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		Point pfree=lp;
		pfree.x+=(*r-freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-freeDelta)*sin(lp.theta+*angle);
		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )>m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
					Point mu=phit-cell.mean();
					if (!found){
						bestMu=mu;
						found=true;
					}else
						bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
				}
			//}
		}
		if (found){
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
			c++;
		}
		if (!skip){
			double f=(-1./m_likelihoodSigma)*(bestMu*bestMu);
			l+=(found)?f:noHit;
		}
	}
	return c;
}
#endif

//new version which uses RangeReading directly
inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p,   const RangeReading& reading) const{
  double s=0;
  OrientedPoint lp=p;

    //extracting data from de rangeSensor pointed by the RangeReading
    const RangeSensor *rangeSensor = dynamic_cast<const RangeSensor*>(reading.getSensor());
    double laserPoseX = rangeSensor->getPose().x;
    double laserPoseY = rangeSensor->getPose().y;
    double laserPoseTheta = rangeSensor->getPose().theta;

  lp.x+=cos(p.theta)*laserPoseX-sin(p.theta)*laserPoseY;
  lp.y+=sin(p.theta)*laserPoseX+cos(p.theta)*laserPoseY;
  lp.theta+=laserPoseTheta;
  unsigned int skip=0;
  double freeDelta=map.getDelta()*m_freeCellRatio;
    //for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
    std::vector<RangeSensor::Beam> beams = rangeSensor->beams();
    for (unsigned int i = m_initialBeamsSkip; i < beams.size(); i++) {
        //TODO implementar skip con incremento i += (m_likelihoodSkip+1)
    skip++;
    skip=skip>m_likelihoodSkip?0:skip; //to use only one out of m_likelihoodSkip beam
    if (skip) continue;

        double currAngle = beams[i].pose.theta;
        double currRay = reading[i];
        //const double * angle=m_laserAngles+m_initialBeamsSkip;
    if (currRay >m_usableRange|| currRay ==0.0) continue;

    Point phit=lp;
    phit.x+= currRay * cos(lp.theta+ currAngle );
    phit.y+= currRay * sin(lp.theta+ currAngle );
    IntPoint iphit=map.world2map(phit);
    Point pfree=lp;
    pfree.x+=( currRay -map.getDelta()*freeDelta)*cos(lp.theta+ currAngle );
    pfree.y+=( currRay -map.getDelta()*freeDelta)*sin(lp.theta+ currAngle );
     pfree=pfree-phit;
    IntPoint ipfree=map.world2map(pfree);
    bool found=false;
    Point bestMu(0.,0.);
    for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
    for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
      IntPoint pr=iphit+IntPoint(xx,yy);
      IntPoint pf=pr+ipfree;
            const PointAccumulator& cell=map.cell(pr);
            const PointAccumulator& fcell=map.cell(pf);
            if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
                Point mu=phit-cell.mean();
                if (!found){
                    bestMu=mu;
                    found=true;
                }else
                    bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
            }
    }
    if (found)
      s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
  }
  return s;
}

//new version which uses RangeReading directly
inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p,  const RangeReading& reading) const{
  using namespace std;
  l=0;
  s=0;
  //const double * angle=m_laserAngles+m_initialBeamsSkip;
  OrientedPoint lp=p;
    //extracting data from de rangeSensor pointed by the RangeReading
    const RangeSensor *rangeSensor = dynamic_cast<const RangeSensor*>(reading.getSensor());
    double laserPoseX = rangeSensor->getPose().x;
    double laserPoseY = rangeSensor->getPose().y;
    double laserPoseTheta = rangeSensor->getPose().theta;

  lp.x+=cos(p.theta)*laserPoseX-sin(p.theta)*laserPoseY;
  lp.y+=sin(p.theta)*laserPoseX+cos(p.theta)*laserPoseY;
  lp.theta+=laserPoseTheta;

  double noHit=nullLikelihood/(m_likelihoodSigma);
  unsigned int skip=0;
  unsigned int c=0;
  double freeDelta=map.getDelta()*m_freeCellRatio;
  //for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
    std::vector<RangeSensor::Beam> beams = rangeSensor->beams();
    for (unsigned int i = m_initialBeamsSkip; i < beams.size(); i++) {
        //TODO implementar skip con incremento i += (m_likelihoodSkip+1)
    skip++;
    skip=skip>m_likelihoodSkip?0:skip;
    if (skip) continue;

        double currAngle = beams[i].pose.theta;
        double currRay = reading[i];
        //const double * angle=m_laserAngles+m_initialBeamsSkip;

    //if (*r>m_usableRange) continue;
    //if (currRay >m_usableRange|| currRay ==0.0) continue;
    if (currRay >m_usableRange) continue;


    Point phit=lp;
    phit.x+= currRay * cos(lp.theta+ currAngle );
    phit.y+= currRay * sin(lp.theta+ currAngle );
    IntPoint iphit=map.world2map(phit);
    Point pfree=lp;
    pfree.x+=( currRay -freeDelta)*cos(lp.theta+ currAngle );
    pfree.y+=( currRay -freeDelta)*sin(lp.theta+ currAngle );
    pfree=pfree-phit;
    IntPoint ipfree=map.world2map(pfree);
    bool found=false;
    Point bestMu(0.,0.);
    for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
    for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
      IntPoint pr=iphit+IntPoint(xx,yy);
      IntPoint pf=pr+ipfree;
            const PointAccumulator& cell=map.cell(pr);
            const PointAccumulator& fcell=map.cell(pf);
            if (((double)cell )>m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
                Point mu=phit-cell.mean();
                if (!found){
                    bestMu=mu;
                    found=true;
                }else
                    bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
            }
    }
    if (found){
      s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
      c++;
    }
    if (!skip){
      double f=(-1./m_likelihoodSigma)*(bestMu*bestMu);
      l+=(found)?f:noHit;
    }
  }
  return c;
}

};

#endif
