/* This is the class for discretizing an integer value of a RV */

#ifndef CPMF_H
#define CPMF_H

//#define _USE_OLD_GAUSSIAN_METHOD_

#include <QVector>
#include <stdlib.h>
#include "fgcommon.h"
#include "utils.h"

class CPMF
{
public:
    CPMF();
    CPMF(FG::pdfType type, int nStates, double valRange, double extRatio = DEFAULT_EXT_RATIO, double gVar = DEFAULT_VAR_FACTOR, int convexness = 0);
    //convexness=0 --> non, convexness=1 --> concave, convexness=2 --> convex
    void setVarPdfType(FG::pdfType newType, int newNStates, double newValRange, double newExtRatio, double newGVar, int convexness = 0);

    bool setVarValue(int value);
    bool setVarValue(double value);
    bool setVarStates(QVector<double> newValue, bool normalize = false);
    int getiVarValue(bool recompute = true);
    double getfVarValue(bool recompute = true);
    //double getVarValue(bool recompute = true);
    int getVarCard() { return _cardinal; }
    QVector<double> getVarStates() {return _pmf;}
    /* for debugging */
    QVector<double> getVarMidVals() {return _midVals;}
    QVector<double> getVarOrgStates() {return _orgPmf;}
    QVector<double> getVarBinWidth() {return _binWidth;}

private:
    FG::pdfType _type;
    int _cardinal;  //cardinality of a variable, i.e. the number of states used for discrete representation
    //int _intValue;  //integer value of the variable
    double _value;  //double value of the variable
    double _valRange;
    double _extValRange;
    QVector<double> _pmf;   //the real probability mass function of each state
    /* for Single pdftype */
    double _m, _b;  //linearity parameters
    /* for Gaussian pdftype */
    QVector<double> _binWidth, _midVals;
    int _OrgStatesLength; //this represents the number of original bins for each state
    int _convexness; //0=non, 1=concave, 2=convex
    double _extRatio;   //to help to cope with the end effect of the slope
    double _gVar, _TwoVariance;
    /* for debugging, the following _orgPmf should be removed later */
    QVector<double> _orgPmf;
public:
    /*______________________ Utility Functions ________________________*/
    static void linspace(double x[], double a, double b, int n);
    double betapdf(double x, double alpa, double beta);
    static double errf(double x);  //compute error function
    double gar(double a, double b, double i);
    double tmf(int x, int nStates); //for computing trapezoidal membership function
};

#endif // CPMF_H
