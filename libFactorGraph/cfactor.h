/*
 *Aku menemukan satu problematik: CFactor membutuhkan variable dengan indexing mulai dari 0, sedangkan
 *aku ingin FG nanti untuk variable mulai dari 1
 */


#ifndef CFACTOR_H
#define CFACTOR_H

#include <QVector>
#include <QString>
#include "utils.h"

/*
 *stuct chn will be used to simplify multiplication of two factors, so that it will not try to read
 *the cardinalities of both factors, but read from the precomputed channel (chn)
 */

class CFactor
{

public:
    //for empty CFactor
    CFactor();                          //this is for creating a dummy factor which contains only value '1', for example for product/multiplication function
    void reset();
    void setID(int newID) {_ID = newID;}
    void setName(QString newName) {_name = newName;}
    void setScope(QVector<int> newScopes) {_scopes = newScopes;}
    void setScope(int newScope) {_scopes.append(newScope);}
    void setCard(QVector<int> newCard) {_card = newCard; _jpd.resize(getElementProd(_card));}
    void setCard(int newCard) {_card.append(newCard);}
    bool setJPD(QVector<double> newJPD);
    int getID() const {return _ID;}
    QString getName() const {return _name;}
    QVector<int> getCard() const {return _card;}
    QVector<double> getJPD() const {return _jpd;}
    QVector<int> getScope() const {return _scopes;}

    void normalize();
    CFactor cond(int varID, int state); //get the conditioned factor on specific state of variable varID
    CFactor getCPT(const CFactor &JPT, QVector<int> givenVarIDs, bool *ok = 0, QString denumFile = QString()); //get conditional probability table on varIDs, ex. Pr(C|AB)
    void getCPT(QVector<int> givenVarIDs, QString denumFile = QString());
    void getCPT(int givenVarID, QString denumFile = QString());
    CFactor summ(CFactor F, int varID);
    void summ(int varID);
    CFactor marg(const CFactor &F, int varID, int method = 1);
    void marg(int varID, int method = 1);
    bool isValid() const { return _valid;} //send the status of JPD, it is exist or uninitialized
    bool  isJPDvalid();
    CFactor prod(const CFactor &F1, const CFactor &F2, bool similarScope = true, int newID = 0, QString newName = QString());
    void prod(const CFactor &F, bool similarScope = true);
    void expand(const CFactor &parent, int missingVar);

private:
    int _ID;
    QString _name;
    QVector<int> _scopes;            //which nodes connected to this factor
    QVector<int> _card;              //when first created, it might be empty. it should be used to accelerate the computation, instead of reading the value from _VarNodes
    QVector<double> _jpd;               //holds the joint probability values for each joint states
    bool _valid;                        //valid means the _jpd is initialized or has values

public:
    /*___________________________________ Public Utilities __________________________________*/
    static bool createDenumFile(QVector<int> vars, QVector<int> cards, QVector<int> givenVarIDs, QString filename, bool showDebugInfo = true);
};

#endif // CFACTOR_H

/* maybe later?
 *
    void addVars(QVector<CVar> varNodes);
    void addVar(CVar varNode);
 *
 **/
