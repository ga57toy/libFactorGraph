#include "cpmf.h"
#include <cmath>
#include <QtGlobal>
#include <stdio.h>

CPMF::CPMF()
{

}

CPMF::CPMF(FG::pdfType type, int nStates, double valRange, double extRatio, double gVar, int convexness)
{
    setVarPdfType(type, nStates, valRange, extRatio, gVar, convexness);
}

void CPMF::setVarPdfType(FG::pdfType newType, int newNStates, double newValRange, double newExtRatio, double newGVar, int convexness)
{
    _type = newType;
    _cardinal = newNStates;
    _valRange = newValRange;
    _extRatio = newExtRatio;
    _gVar = newGVar;
    _convexness = convexness; //0=non, 1=concave, 2=convex
    _pmf.resize(_cardinal);
    _midVals.resize(_cardinal);
    _binWidth.resize(_cardinal);

    switch(newType)
    {
    case FG::Single:
        _m = (_cardinal-1)/(2*_valRange);
        _b = (_cardinal-1)-(_m*_valRange);
        break;
    case FG::Gaussian:
    case FG::Beta:
        _extValRange = _extRatio*_valRange;
        double Variance;
        Variance = _gVar*_valRange;                      //standard variance of the gaussian distribution
        //_OrgStatesLength = 2*_extValRange;                //number of states in the discretized original value
        int Panjang;
        Panjang = int(round(2*_extValRange));
        _OrgStatesLength = Panjang > _cardinal ? Panjang:_cardinal;
        _TwoVariance = 2*Variance;                              //just to simplify the formula

        //now, the property of convexness determine the binWidth and midVals
        if(_convexness==0) //non-convex
        {
            /* beware with the following syntax:
            _bindWidth = round(_OrgStatesLength / _cardinal);
            * since it will produce wrong result because (_OrgStatesLength / _cardinal) will be treated as integer first, which is 0, then converted into double by the round()
            */
            for(int i=0; i<_cardinal; i++)
                _binWidth[i] = double(_OrgStatesLength) / double(_cardinal);
        } else
        {
            QVector<double> tBW; //temporary binWidth
            for(int i=0; i<_cardinal+2; i++)
                tBW.append(tmf(i,_cardinal+2));
            tBW.remove(_cardinal+1,1); //kemudian hilangkan bin di belakang
            tBW.remove(0,1); //dan di depan

            double tBWsum;
            tBWsum = sum(tBW,0,_cardinal-1);
            for(int i=0; i<_cardinal; i++) _binWidth[i] = tBW.at(i)*double(_OrgStatesLength)/tBWsum;
            tBWsum = sum(_binWidth, 0, _cardinal-1);
            QVector<double> newBinWidth;
            newBinWidth.resize(_binWidth.count());

            if(_convexness==1) //concave
            {
                //newBinWidth = _binWidth;
            } else //convex
            {
                int shiftPoint,n;
                shiftPoint = ceil(_cardinal/2);

                n = _cardinal%2;
                switch(n)
                {
                case 0: //jumlah state = genap
                    for(int i=0; i<shiftPoint; i++)
                    {
                        newBinWidth[i] = _binWidth.at(i+shiftPoint);
                        newBinWidth[i+shiftPoint] = _binWidth.at(i);
                    }
                    break;
                default: //jumlah state = ganjil
                    for(int i=0; i<shiftPoint; i++)
                    {
                        newBinWidth[i] = _binWidth.at(i+shiftPoint-1);
                        newBinWidth[i+shiftPoint-1] = _binWidth.at(i);
                    }
                }
                _binWidth = newBinWidth;
            }
        }
        double startPoint, lowerRange, upperRange;
        startPoint = -_extValRange;
        for(int i=0; i<_cardinal; i++)
        {
            lowerRange=startPoint;
            upperRange=startPoint+_binWidth.at(i);
            _midVals[i] = (lowerRange+upperRange)/2.0;
            startPoint = startPoint + _binWidth.at(i);
        }

        break;
    case FG::Binary:
        break;
    }
}

bool CPMF::setVarValue(int value)
{
    return setVarValue(double(value));
}

bool CPMF::setVarValue(double value)
{
    if(abs(value)<=_valRange)
    {
        _value=value;
        int ActiveState; //hanya untuk Single
        int i;
        double border1, border2; //untuk Gaussian
        double lowerPart, upperPart;
        /* then distribute the states */
        switch(_type)
        {
        case FG::Single:
        {
            ActiveState = int(round(_m*(double)value+_b));
            _pmf.fill(0.0); _pmf[ActiveState]=1.0;
        }
            break;
        case FG::Gaussian:
        {
#ifdef _USE_OLD_GAUSSIAN_METHOD_
            /* OLD METHOD: */
            int arg;
            QVector<double> OriginalBins(_OrgStatesLength);
            //Lesson to learn: hati-hati dengan pow(), jadi sebaiknya bikin variable untuk menampung sementara operand sebelum pow
            //for(i=0; i<_OrgStatesLength; i++) OriginalBins[i] = exp(-pow(i-_extValRange-value,2)/_TwoVariance)/sqrt(_PI*_TwoVariance);
            //int arg, pangkat2, pangkat2pow, pangkat2powarg;
            //double expo;
            //printf("_extValRange = %d\n", _extValRange);
            //for(i=0; i<_OrgStatesLength; i++)
            //{
            //    arg = i-_extValRange-value;
            //    pangkat2 = -1.0*(arg)*(arg);
            //    pangkat2pow = -1.0*pow(double(i-_extValRange-value),2.0);
            //    pangkat2powarg = -1.0*pow(arg,2);
            //    expo = exp(double(pangkat2)/_TwoVariance);
            //    printf("i = %d, arg = %d, Pangkat2 = %d, Pangkat2pow = %d, Pangkat2powarg = %d, expo = %lf\n", i, arg, pangkat2, pangkat2pow, pangkat2powarg, expo);
            //    OriginalBins[i] = expo/sqrt(_PI*_TwoVariance);
            //}
            for(i=0; i<_OrgStatesLength; i++)
            {
                //berikut cara matlab
                //arg = i-_extValRange-value; //lihat penjelasan di atas, akan bermasalah jika tidak dimasukkan dalam temporary variable
                //maka supaya hasilnya benar:
                arg = i-_extValRange-value+1; //+1 supaya indexingnya sama seperti yang di matlab
                OriginalBins[i] = exp(-1.0*pow(arg,2)/_TwoVariance)/sqrt(_PI*_TwoVariance);
            }

            /* for debugging, the following _orgPmf should be removed later *
            //_orgPmf = OriginalBins;
            /*______________________________________________________________*/

            for(i=0; i<_cardinal; i++)
            {
                //ada masalah dengan perhitungan border!!!
                //border1 = (i-1)*binWidth+1;                      %here plays the role of offset
                //border2 = min(i*binWidth,length(OriginalBins)); %TryToAvoidIndexExceedMatrixDimension_offset

                border1 = i*_binWidth;
                //berikut ini salah:
                //border2 = (i+1)*_binWidth > _OrgStatesLength-1 ? _OrgStatesLength-1:(i+1)*_binWidth;
                //yang benar:
                border2 = ((i+1)*_binWidth-1) > (_OrgStatesLength-1) ? _OrgStatesLength-1:(i+1)*_binWidth-1;
                _pmf[i] = sum(OriginalBins,border1,border2);
            }
#else
            //Lesson to learn: hati-hati dengan pow(), jadi sebaiknya bikin variable untuk menampung sementara operand sebelum pow...
            //lihat OLD METHOD untuk kelanjutannya
            if(_convexness!=0)
            {
                //ini adopsi dari jpdGetStates3.m
                double startPoint, border1, border2;
                startPoint = -_extValRange;
                for(int i=0; i<_cardinal; i++)
                {
                    border1 = startPoint;
                    border2 = border1+_binWidth.at(i);

                    _pmf[i] = gar(border1,border2,value);
                    startPoint = startPoint+_binWidth.at(i);
                }
            } else
            {
                //cara lama, tidak aku ubah: males ah!
                for(i=0; i<_cardinal; i++)
                {

                    border1 = i*_binWidth.at(i)-_extValRange;
                    border2 = border1+_binWidth.at(i)-1;

                    _pmf[i] = gar(border1,border2,value);
                }
            }
#endif
        }
            _pmf = normalizeStates(_pmf);
            break;
        case FG::Binary:
        {
            lowerPart = (value + _valRange) / (2.0*_valRange);
            upperPart = 1.0-lowerPart;
            _pmf[0] = lowerPart;
            _pmf[1] = upperPart;
        }
            break;
        case FG::Beta:
        {
            double a, b;
            if(value==0)
            {
                a = _extValRange;
                b = _extValRange;
            } else if(value > 0)
            {
                a = _extValRange;
                b = _extValRange - value;
            } else
            {
                a = _extValRange + value;
                b = _extValRange;
            }
            double x[_OrgStatesLength];
            linspace(x, 0, 1, _OrgStatesLength);
            _orgPmf.clear();
            for(int i=0; i<_OrgStatesLength; i++)
            {
                if(i==_OrgStatesLength-1)
                    int coba = i;
                _orgPmf.append(betapdf(x[i], a, b));
            }
            for(int i=0; i<_cardinal; i++)
            {
                border1 = i*_binWidth.at(i);
                border2 = ((i+1)*_binWidth.at(i)-1) > (_OrgStatesLength-1) ? _OrgStatesLength-1:(i+1)*_binWidth.at(i)-1;
                _pmf[i] = sum(_orgPmf,border1,border2);
            }
            _pmf = normalizeStates(_pmf);
        }
            break;
        }
        return true;
    } else return false;
}

/*
 * setVarStates() will set the states directly
 */
bool CPMF::setVarStates(QVector<double> newValue, bool normalize)
{
    if(newValue.count()!=_pmf.count())
        return false;
    else
    {
        if(normalize)
            _pmf = normalizeStates(newValue);
        else
            _pmf = newValue;
        return true;
    }
}

int CPMF::getiVarValue(bool recompute)
{
    double hasil = getfVarValue(recompute);
    int result = int(hasil);
    return result;
}

double CPMF::getfVarValue(bool recompute)
{
    double result = 0.0;
    /* if needs to recompute from the _pmf */
    if(recompute)
    {
        switch(_type)
        {
        case FG::Single:
            int idx; double mx;
            idx = 0; mx = _pmf.at(0);
            for(int i=1; i<_cardinal; i++) if(_pmf.at(i)>mx) { mx=_pmf.at(i); idx=i; }
            //Q_ASSERT_X(idx!=-1,"Invalid PMF","Cannot find active state!");
            result = (double(idx)-_b)/_m;
            break;
        case FG::Gaussian:
            //        for i=1:Nstates
            //            %use EXT_VALUE_RANGE for computing the middle point
            //            Rnum = Rnum + States(i) * midValues(i);
            //        end
            result = DEFAULT_ADJUSTMENT;
            //hati-hati, berikut ini menghasilkan nilai yang salah:
            //for(int i=0; i<_cardinal; i++) result += round(_pmf.at(i) * _midVals.at(i));
            for(int i=0; i<_cardinal; i++) result += _pmf.at(i) * _midVals.at(i);
            //        double angka1, angka2, kali;
            //        for(int i=0; i<_cardinal; i++)
            //        {
            //            angka1 = _pmf.at(i);
            //            angka2 = _midVals.at(i);
            //            kali = angka1*angka2;
            //            result = result + kali;
            //        }
            break;
        case FG::Binary:
            result = 2.0*double(_valRange)*_pmf[0] - double(_valRange);
            break;
        case FG::Beta:
        {
            double x[_cardinal], ev, vr, a, b, max_v, m, c;
            linspace(x,0,1,_cardinal);
            ev = 0;
            for(int i=0; i<_cardinal; i++)
                ev+=_pmf.at(i)*x[i];
            vr = 0;
            for(int i=0; i<_cardinal; i++)
                vr+=_pmf.at(i)*pow((x[i]-ev),2);
            a = ev*(ev*(1-ev)/vr - 1);
            b = (1-ev)*(ev*(1-ev)/vr - 1);
            max_v = a>b?a:b;
//            m = double(_extValRange)/max_v;
//            c = double(_extValRange) - m*max_v;
//            result = a>b?double(_extValRange) - (m*b+c):(m*a+c) - double(_extValRange);
            m = _valRange/max_v;
            c = _valRange - m*max_v;
            result = a>b?_valRange - (m*b+c):(m*a+c) - _valRange;
        }
            break;
        }
        //kalo yang ini salah:
        //return int(result);
    }
    return result;
}

/*______________________ Utility Functions ________________________*/

double CPMF::betapdf(double x, double alpha, double beta)
{
    double loga,logb,betaFunc,result;
    betaFunc = lgamma(alpha)+lgamma(beta)-lgamma(alpha+beta);
    loga = (alpha-1)*log(x);
    logb = (beta-1)*log(1-x);
    result = exp(loga+logb - betaFunc);
    return result;
}

void CPMF::linspace(double x[], double a, double b, int n)
{
    double interval = (b-a)/double(n-1);
    x[0] = a;
    for(int i=1; i<n; i++)
        x[i] = x[i-1]+interval;
    //we found a problem with betapdf when the last item greater than 1.0,
    //it happens when computing logb, so we have to make sure that the final item is equal the limit
    x[n-1] = b;
}

double CPMF::errf(double x)
{
    double result;
    double t = 1/(1+0.5*fabs(x));
    double x2 = pow(x,2);
    double t2 = t*t;
    double t3 = t*t2;
    double t4 = t*t3;
    double t5 = t*t4;
    double t6 = t*t5;
    double t7 = t*t6;
    double t8 = t*t7;
    double t9 = t*t8;
    double tau = t*exp(-x2 - 1.26551223 + 1.00002368*t + 0.37409196*t2 + 0.09678418*t3 - 0.18628806*t4 + 0.27886807*t5 - 1.13520398*t6 + 1.48851587*t7 - 0.82215223*t8 + 0.17087277*t9);
    if(x>=0) result = 1.0 - tau; else result = tau -1.0;
    return result;
}

double CPMF::gar(double a, double b, double i)
{
    double tv = _TwoVariance;
    double y = 0.5*(errf(3.0*(b-i)/sqrt(2.0*tv)) - errf(3.0*(a-i)/sqrt(2.0*tv)));
    return y;
}

double CPMF::tmf(int x, int nStates)
{
    //INGAT: indexing dimulai dari 0 hingga nStates-1
    int a,b,c,d,n;
    double y;

    a = 0;
    d = nStates-1;
    n = nStates%2;

    switch(n)
    {
    case 0: //jumlah state = genap
        b = nStates/2 - 1; c = b+1;
        break;
    default: //jumlah state = ganjil
        b = floor(nStates/2); c = b;
    }
    if(x<a || x>d) y = 0.0;
    else if(x<b) y = double(x-a)/double(b-a);
    else if(x<=c) y = 1.0;
    else y = double(d-x)/double(d-c);
    return y;
}

