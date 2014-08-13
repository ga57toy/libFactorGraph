#ifndef CNODE_H
#define CNODE_H

/*
 * Convention: Positive-ID for Variable, Negative-ID for Factor
 * Note: in the beginning, every node has _localF with value 1
 * the _localF has also the same ID as this node
 *
 */

#include <QObject>
#include <QVector>
#include "fgcommon.h"
#include "cmessage.h"

class CNode : public QObject
{
    Q_OBJECT
public:
    CNode();
    CNode(FG::NodeType type, int ID, int card = DEFAULT_CARD, QObject *parent = 0);
    ~CNode();
    void init(FG::NodeType type, int ID, int card = 0);
    void addNeighbor(int nID);
    FG::NodeType getType() const {return _type;}
    int getID() const {return _ID;}
    /* for factor type node */
    void initFactor(QVector<int> vars, QVector<int> cards);
    bool setFactor(QVector<double> newJPD) {return _localF.setJPD(newJPD);}
    void setFactor(QVector<int> vars, QVector<int> cards, QVector<double> vals);
    int getCard() const {return _card;}
    CFactor getFactor() const {return _localF;}
    bool checkScopeOrder();
signals:

public slots:
    bool putMessage(CMessage msg);
    CMessage computeMessage(bool *success = 0);             //for a node as leaf
    CMessage computeMessage(int destID, bool *success = 0); //for a node which is not as leaf

private:
    FG::NodeType _type;
    int _ID;
    QVector<int> _neighbor;
    //the idea: for every neighbor, there exist n-vectors containing input messages
    QVector< QVector<CMessage> > msgTable;  //this is 2D array: row reflect the destination, while cols are the source nodes
    //the first msgTable element belongs to the first element in _neighbor, etc
    CFactor _localF;    //for node as a factor node
    int _card; //only for variable node
    QVector<CMessage> sortMsg(const QVector<CMessage> &msg, const QVector<int> &nb);
};

#endif // CNODE_H
