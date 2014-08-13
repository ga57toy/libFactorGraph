
#include "utils.h"
#include <QtAlgorithms>
#include <QString>
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QDebug>
#include <cmath>


/* This function will compute the index of a jpd vector, given it's assignment
   with cardinal is 1xn vector and assigment is also 1xn vector
   note: the left most variable is the most repetitive one */

/* Note: if strideIncluded is true, then vect is the stride, otherwise it is the card */
int getJPDidx(QVector<int> vect, QVector<int> A, bool strideIncluded) //A=assignment
{
    int Nvar = vect.count();
    int Sum_over_i = 0;
    if(!strideIncluded)
    {
        QVector<int> stride(Nvar, 1);
        /* First, build the stride vector */
        for(int i=1; i<Nvar; i++)
        {
            stride[i] = vect.at(i-1) * stride.at(i-1);
            Sum_over_i += A.at(i)*stride.at(i);
        }
        Sum_over_i+=A.at(0)*stride.at(0);
//        for(int i=0; i<Nvar; i++)
//            Sum_over_i += A.at(i)*stride.at(i);
    }
    else for(int i=0; i<Nvar; i++) Sum_over_i += A.at(i)*vect.at(i);
    return Sum_over_i;
}

/* This function will compute the assignment, given the index and the cardinal
   Rumus aslinya: assignment[i] = (index/stride[i]) mod cardinal[i] */
QVector<int> getAssignment(QVector<int> card, int idx)
{
    int Nvar = card.count();
    QVector<int> stride(Nvar, 1);

    /* First, build the stride vector */
    for(int i=1; i<Nvar; i++)
        stride[i] = card.at(i-1) * stride.at(i-1);

    QVector<int> A(Nvar);
    for(int i=0; i<Nvar; i++)
        A[i] = idx/stride.at(i) % card.at(i);
    return A;
}

QVector<int> getAssignment(QVector<int> stride, QVector<int> card, int idx)
{
    int Nvar = stride.count();
    QVector<int> A(Nvar);
    for(int i=0; i<Nvar; i++)
        A[i] = idx/stride.at(i) % card.at(i);
    return A;
}

QVector<int> getUnion(QVector<int> S1, QVector<int> S2, bool sorted)
{
    QVector<int> result = S1;
    for(int i=0; i<S2.count(); i++)
        if(!S1.contains(S2.at(i))) result.append(S2.at(i));

    if(sorted) qStableSort(result);
    return result;
}

QVector<int> getIntersect(QVector<int> S1, QVector<int> S2, bool sorted)
{
    QVector<int> result;
    for(int i=0; i<S1.count(); i++)
        if(S2.contains(S1.at(i))) result.append(S1.at(i));
    if(sorted) qStableSort(result);
    return result;
}

QVector<int> getDiff(QVector<int> V1, QVector<int> V2)
{
    QVector<int> result = V1;
    int pos;
    for(int i=0; i<V2.count(); i++)
    {
        pos = result.indexOf(V2.at(i));
        if(pos != -1) result.remove(pos, 1);
    }
    return result;
}

QVector<double> normalizeStates(QVector<double> pmf)
{
    QVector<double> result(pmf.count());
    double normFactor = sum(pmf,0, pmf.count()-1);
    for(int i=0; i<pmf.count(); i++) result[i]=pmf.at(i)/normFactor;
    return result;
}

/*
 * Create a denumerator table used for computing conditional probability, such as P(B|AC) = P(ABC)/[sum_all_B{P(ABC)}]
 **/
QVector< QVector<int> > getAssignmentDenum(const QVector< QVector<int> > &A, int excludedVarPos, bool validate, bool *ok, QString denumFile, bool showDebugInfo)
{
    bool OK = true;
    QVector< QVector<int> > denumTable;
    bool createNew = false;
    bool saveToFile = false;
    QString line;
    int cnt;

    if(!denumFile.isEmpty() || !denumFile.isNull())
    {
        QStringList items;
        QFile fid(denumFile);
        int cntr = 0;
        if(fid.exists())
        {
            fid.open(QIODevice::ReadOnly | QIODevice::Text);
            while(!fid.atEnd())
            {
                line = fid.readLine();
                items = line.split(',');
                if(items.count()>0)
                {
                    cntr++;
                    denumTable.resize(cntr);
                    for(int i=0; i<items.count(); i++)
                        denumTable[cntr-1].append(items.at(i).toInt());
                }
            }
            fid.close();
            /* if denumTable from file doesn't match the length of assignment table A, then create new */
            if(A.count()!=denumTable.count())
            {
                qDebug() << "DenumTable from file doesn't match assignment table A. Hence a new denumTable will be created!";
                createNew = true;
            }
        }
        else
        {
            createNew = true; //if denumFile is specified but doesn't exist
            saveToFile = true; //then we will create it and save it to a file
        }
    }
    else createNew = true;

    if(createNew)
    {
        cnt = A.count();            //how many rows of A?
        denumTable.resize(cnt);

        QVector<int> checkA, epochA;
        for(int i=0; i<cnt; i++)
        {
            if(showDebugInfo)
            {
                qDebug() << QString("Processing index %1 out of %2 (%3%%)").arg(i).arg(cnt).arg(i*100/cnt);
            }
            if(denumTable.at(i).isEmpty())
            {
                checkA = A.at(i);
                checkA.remove(excludedVarPos,1);
                for(int j=0; j<cnt; j++)
                {
                    epochA = A.at(j);
                    epochA.remove(excludedVarPos,1);
                    if(checkA == epochA)
                    {
                        denumTable[i].append(j);
                    }
                }
                /* then copy the resulted denumTable[i] to all equivalent position to save time */
                for(int k=0; k<cnt; k++)
                {
                    if(denumTable.at(i).indexOf(k)!=-1 && denumTable.at(k).isEmpty())
                        for(int l=0; l<denumTable[i].count(); l++)
                            denumTable[k].append(denumTable[i].at(l));
                }
            }
        }
        if(saveToFile)
        {
            QFile fid(denumFile);
            QTextStream fOut(&fid);
            fid.open(QIODevice::WriteOnly | QIODevice::Text);
            for(int i=0; i<denumTable.count(); i++)
            {
                line.clear();
                for(int j=0; j<denumTable[i].count(); j++)
                    line.append(QString("%1,").arg(denumTable.at(i).at(j)));
                line.remove(line.length()-1,1); //remove the last comma
                fOut << QString("%1\n").arg(line);
            }
            fid.close();
        }
    }

    if(validate)
        for(int i=1; i<cnt; i++)
            if(denumTable.at(i).count()!=denumTable.at(i-1).count())
            {
                OK = false;
                qFatal("Found mismatch!");
            }

    if(ok!=0) *ok = OK;
    return denumTable;
}
