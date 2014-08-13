#include "cfactor.h"
#include <QtGlobal>
#include <QDebug>

CFactor::CFactor()
{
    reset();
}

void CFactor::reset()
{
    _scopes.clear(); _card.clear(); _jpd.clear(); _jpd.append(1.0); _name.clear(); _valid = false;
}

void CFactor::normalize()
{
    QVector<double> tmp;
    tmp = _jpd;
    tmp = normalizeStates(tmp);
    _jpd = tmp;
}

bool CFactor::setJPD(QVector<double> newJPD)
{
    /* check consistency with cardinality */
    int jpdSize = getElementProd(_card);
    if(jpdSize == newJPD.count())
    {
        _jpd = newJPD; _valid = true; return true;
    }
    else
    {
        qFatal("Invalid values: Make sure the number of elements reflects the correct cardinality!");
        return false;
    }
}

/*
 * To see if the JPD is a valid probabilistic entry, which sums to 1
 */
bool CFactor::isJPDvalid()
{
    double TOLERANCE = 0.001;
    double sum = 0;
    for(int i=0; i<_jpd.count(); i++) sum+=_jpd.at(i);
    if(abs(sum-1.0)<=TOLERANCE) return true; else return false;
}

CFactor CFactor::prod(const CFactor &F1, const CFactor &F2, bool similarScope, int newID , QString newName)
{
    CFactor result;
    if(similarScope)
    {
        if(F1.getScope()!=F2.getScope())
            qFatal("Mismatch scope");
        result.setScope(F1.getScope());
        result.setCard(F1.getCard());
        QVector<double> elProd(result.getJPD().count());
        QVector<double> jpd1 = F1.getJPD(), jpd2 = F2.getJPD();
        for(int i=0; i<result.getJPD().count(); i++)
            elProd[i] = jpd1.at(i)*jpd2.at(i);
        result.setJPD(elProd);
        result.setID(newID);
        result.setName(newName);
    }
    else
    {
        result.setScope(F1.getScope()+F2.getScope());
        result.setCard(F1.getCard()+F2.getCard());
        QVector<double> elProd(result.getJPD().count());
        QVector<double> jpd1 = F1.getJPD(), jpd2 = F2.getJPD();
        int jpd1Length = jpd1.count();
        for(int i=0; i<jpd2.count(); i++)
            for(int j=0; j<jpd1.count(); j++)
                elProd[j+(i*jpd1Length)] = jpd1.at(j)*jpd2.at(i);
        result.setJPD(elProd);
        result.setID(newID);
        result.setName(newName);
    }
    return result;
}

/*
 *Multiply *this* factor with another factor F
 */
void CFactor::prod(const CFactor &F, bool similarScope)
{
    int myID = this->getID();
    QString myName = this->getName();
    CFactor me = prod(*this, F, similarScope, myID, myName);
    *this = me;
}

CFactor CFactor::cond(int varID, int state)
{
    CFactor result;
    int idx_varID = this->_scopes.indexOf(varID);                           //get the position of varID in the _VarNodesID
    if(idx_varID!=-1)
    {
        QVector<int> resVars = this->_scopes; resVars.remove(idx_varID);    //remove the varID in the result
        QVector<int> resCard = this->_card; resCard.remove(idx_varID);          //remove its cardinal as well
        result.setScope(resVars);
        result.setCard(resCard);
        QVector<int> A;
        QVector<double> resVal;
        for(int i=0; i<this->_jpd.count(); i++)                                 //scan through the entire table
        {
            A = getAssignment(this->_card, i);                                  //get the assignment formation for the current index
            if(A.at(idx_varID)==state)                                          //if the assignment contains the current state of varID
                resVal.append(this->_jpd.at(i));                                //then copy the content of the table to the result
        }
        result.setJPD(resVal);
    }
    return result;
}

CFactor CFactor::summ(CFactor F, int varID)
{
    CFactor result;
    int idx_varID = F.getScope().indexOf(varID);                                 //get the position of varID in the _VarNodesID
    if(idx_varID!=-1)
    {
        QVector<int> resVars = F.getScope(); resVars.remove(idx_varID);    //remove the varID in the result
        QVector<int> resCard = F.getCard(); resCard.remove(idx_varID);          //remove its cardinal as well
        int valSz = getElementProd(resCard);
        QVector<double> resVal(valSz);                                          //prepare the container of cpt's value
        result.setScope(resVars);
        result.setCard(resCard);
        int i, j, numOfCPT = F.getCard().at(idx_varID);
        QVector<CFactor> cpt(numOfCPT);
        for(i = 0; i<numOfCPT; i++)
            cpt[i] = F.cond(varID, i);                                         //create cpts
        for(j = 0; j<valSz; j++)
            for(i = 0; i<numOfCPT; i++)
                resVal[j]+=cpt.at(i).getJPD().at(j);                            //combine values from all cpts
        result.setJPD(resVal);
    }
    return result;
}

void CFactor::summ(int varID)
{
    CFactor me = summ(*this, varID);
    *this = me;
}

/*
 * marg() marginalize factor F over a variable with varID
 */
CFactor CFactor::marg(const CFactor &F, int varID, int method)
{
    CFactor result = F;                                                     //copy the original factor
    int idx_varID = _scopes.indexOf(varID);                                 //get the position of varID in the _VarNodesID
    if(idx_varID!=-1)
    {
        /* method 0 using iterative standard summ() */
        if(method==0)
        {
            QVector<int> resVars = this->_scopes; resVars.remove(idx_varID);    //remove the varID in the result
            for(int i=0; i<resVars.count(); i++) result = summ(result, resVars.at(i));
            //for(int i=resVars.count()-1; i>=0; i--) result = result.summ(resVars.at(i)); -> just the same!
        }
        /* method 1 using sorting mechanism */
        else if(method==1)
        {
            QVector<int> newScope;
            QVector<double> newJPD(result.getJPD().count()), resultJPD(result.getCard().at(idx_varID));
            QVector<int> card = result.getCard(), newCard = card;

            int lastPos = _scopes.count()-1, lastStride;
            int Nvar = lastPos+1;
            double total;
            QVector<int> stride(Nvar, 1), newStride(Nvar,1);

            // first, swap the card and re-compute the stride
            swapVectorElement(newCard, idx_varID, lastPos);
            /* then build the stride vector */
            for(int i=1; i<Nvar; i++)
            {
                stride[i] = card.at(i-1) * stride.at(i-1);
                newStride[i] = newCard.at(i-1) * newStride.at(i-1);
            }

            if(idx_varID!=lastPos) //if the input varID is not the right-most variable, then sort it first
            {
                QVector<int> A;
                for(int i=0; i<result.getJPD().count(); i++)
                {
                    A = getAssignment(stride, card, i);
                    swapVectorElement(A, idx_varID, lastPos);
                    newJPD[getJPDidx(newStride, A, true)] = result.getJPD().at(i);
                }
            }
            else
                newJPD = result.getJPD();
            //then do collecting based on the last column...
            lastStride = newStride.at(lastPos);
            for(int i=0; i<resultJPD.count(); i++)
            {
                total = 0.0;
                for(int j=0; j<lastStride; j++)
                    total+=newJPD.at(i*lastStride+j);
                resultJPD[i] = total;
            }
            newScope.append(varID);
            //result.setCard(card);
            //the following will not work: result.setCard(newCard.at(newCard.count()-1)); because the "append" feature of this function
            QVector<int> c; c << newCard.at(newCard.count()-1);
            result.setCard(c);
            result.setScope(newScope);
            result.setJPD(resultJPD);
        }
    }
    return result;
}

void CFactor::marg(int varID, int method)
{
    int orgID = this->getID();
    QString orgName = this->getName();
    CFactor me = marg(*this, varID, method);
    *this = me;
    this->setID(orgID);
    this->setName(orgName);
}

void CFactor::expand(const CFactor &parent, int missingVar)
{
    QVector<int> stride(parent.getScope().count());
    QVector<int> card = parent.getCard();
    stride[0] = 1;
    for(int i=1; i<parent.getCard().count(); i++)
        stride[i] = card.at(i-1) * stride.at(i-1);
    QVector<double> emptyJPD(getElementProd(card));
    CFactor result;
    result.setScope(parent.getScope());
    result.setCard(parent.getCard());
    QVector<double> myJPD = this->getJPD();
    int myID = this->getID();
    QString myName = this->getName();
    int myIdx = parent.getScope().indexOf(missingVar);
    int myStride = stride.at(myIdx);
    int myCard = card.at(myIdx);
    int resultIdx = 0;
    for(int strideLoop=0; strideLoop < this->getJPD().count(); strideLoop+=myStride)
        for(int cardLoop=0; cardLoop < myCard; cardLoop++)
        {
            for(int i=0; i<myStride; i++)
            {
                emptyJPD[resultIdx] = myJPD.at(strideLoop+i);
                resultIdx++;
            }
        }
    result.setJPD(emptyJPD);
    *this = result;
    this->setID(myID);
    this->setName(myName);
}

bool CFactor::createDenumFile(QVector<int> vars, QVector<int> cards, QVector<int> givenVarIDs, QString filename, bool showDebugInfo)
{
    /*
     * IMPORTANT: At this moment, we make conditioning on several givenVarIDs
     * BUT only for 1 remaining variable
     * Jadi, diff-nya harus cuman 1
     * Juga pastikan bahwa filename tidak exist, karena kalo exist maka getAssignmentDenum()
     * tidak akan menyimpannya
     * TODO: cari dokumentasi kertasku!
     */
    bool bOK = true;
    int rows = getElementProd(cards);
    QVector< QVector<int> > A(rows); //2D array: #rows = #vals, #cols = #vars
    int excludedVar;
    int excludedVarPos;
    excludedVar = getDiff(vars, givenVarIDs).at(0);   //just take 1 variable, see note above!
    excludedVarPos = vars.indexOf(excludedVar);

    /* 1st: build complete assignment table */
    QVector<int> ass;
    for(int i=0;i<rows; i++)
    {
        qDebug() << QString("Building complete assignment table: %1%%").arg(i*100/rows);
        A[i].resize(vars.count());
        ass = getAssignment(cards, i);
        for(int j=0; j<vars.count(); j++)
            A[i][j]=ass.at(j);
    }
    getAssignmentDenum(A, excludedVarPos, true, &bOK, filename, showDebugInfo);
    return bOK;
}

/*
 * At this moment, we make getCPT just for conditioning on several givenVarIDs but only for 1 remaining variable
 */
CFactor CFactor::getCPT(const CFactor &JPT, QVector<int> givenVarIDs, bool *ok, QString denumFile)
{
    bool bOK = true;
    QVector<int> cards = JPT.getCard();
    QVector<double> currentJPD = JPT.getJPD();
    int rows = currentJPD.count();
    //int rows = getElementProd(cards);
    QVector< QVector<int> > A(rows); //2D array: #rows = #vals, #cols = #vars
    QVector<int> vars = JPT.getScope();
    int excludedVar;
    int excludedVarPos;
    excludedVar = getDiff(JPT.getScope(), givenVarIDs).at(0);   //just take 1 variable, see note above!
    excludedVarPos = vars.indexOf(excludedVar);

    /* 1st: build complete assignment table */
    QVector<int> ass;
    for(int i=0;i<rows; i++)
    {
        A[i].resize(vars.count());
        ass = getAssignment(cards, i);
        for(int j=0; j<vars.count(); j++)
            A[i][j]=ass.at(j);
    }

    /* 2nd: build the second list of denumerator table (reduced assignment table) */
    QVector< QVector<int> > denumTable = getAssignmentDenum(A, excludedVarPos, true, &bOK, denumFile);

    /* perform vector division */
    double denum;
    QVector<double> resultedJPD(rows);
    QVector<int> denumIdx;
    for(int i=0; i<rows; i++)
    {
        denumIdx = denumTable.at(i);
        denum = 0.0;
        for(int j=0; j<denumIdx.count(); j++)
            denum += currentJPD.at(denumIdx.at(j));
        if(denum!=0)
            resultedJPD[i] = currentJPD.at(i)/denum;
        else
            resultedJPD[i] = 0;
    }
    CFactor result;
    result.setCard(cards);
    result.setScope(vars);
    result.setJPD(resultedJPD);
    if(ok!=0) *ok = bOK; //always check if it is not a NULL pointer!
    return result;
}

void CFactor::getCPT(QVector<int> givenVarIDs, QString denumFile)
{
    bool ok;
    int myID = getID();
    QString myName = getName();
    CFactor me = getCPT(*this, givenVarIDs, &ok, denumFile);
    if(!ok)
        qFatal("Wrong Conditional Probability Table computation!");
    else
    {
        *this = me;
        setID(myID);
        setName(myName);
    }
}

void CFactor::getCPT(int givenVarID, QString denumFile)
{
    QVector<int> givenVarIDs;
    givenVarIDs << givenVarID;
    getCPT(givenVarIDs, denumFile);
}
