#ifndef UTILS_H
#define UTILS_H

#include <QVector>
#include <QString>

int getJPDidx(QVector<int> vect, QVector<int> A, bool strideIncluded = false);      //A=assignment
QVector<int> getAssignment(QVector<int> card, int idx);
QVector<int> getAssignment(QVector<int> stride, QVector<int> card, int idx);
QVector<int> getUnion(QVector<int> S1, QVector<int> S2, bool sorted = true);
QVector<int> getIntersect(QVector<int> S1, QVector<int> S2, bool sorted = true);
QVector<int> getDiff(QVector<int> V1, QVector<int> V2);
QVector<double> normalizeStates(QVector<double> pmf);
QVector< QVector<int> > getAssignmentDenum(const QVector< QVector<int> > &A, int excludedVarID, bool validate = false, bool *ok = 0, QString denumFile = QString(), bool showDebugInfo = false);

template <typename T> T sum(QVector<T> values, int sp, int ep)
{
    T res=(T)0;
    for(int i=sp; i<ep+1; i++) res+=values.at(i);
    return res;
}

/* beware, avg on integer will user floor() */
template <typename T> T avg(QVector<T> values)
{
    T s = sum(values, 0, values.count()-1);
    T res = s / values.count();
    return res;
}

template <typename T> T getElementProd(QVector<T> vekt)
{
    T res = 1;
    for(int i=0; i<vekt.count(); i++) res*=vekt.at(i);
    return res;
}

template <typename T> T getElementSum(QVector<T> vekt)
{
    T res = 0;
    for(int i=0; i<vekt.count(); i++) res+=vekt.at(i);
    return res;
}

template <typename T> QVector<T> substVector(QVector<T> vekt, T element)
{
    QVector<T> res = vekt;
    int pos = res.indexOf(element);
    if(pos!=-1) res.remove(pos);
    return res;
}

template <typename T> QVector<T> addVectorElements(QVector<T> v1, QVector<T> v2)
{
    Q_ASSERT_X(v1.count()==v2.count(), "Size mismatch", "Both vector have different size!");
    QVector<T> res(v1.count());
    for(int i=0; i<v1.count(); i++)
        res[i] = v1.at(i)+v2.at(i);
    return res;
}

template <typename T> void swapVectorElement(QVector<T> &v, int idx1, int idx2)
{
    T t = v.at(idx1);
    v[idx1] = v.at(idx2);
    v[idx2] = t;
}

#endif // UTILS_H
