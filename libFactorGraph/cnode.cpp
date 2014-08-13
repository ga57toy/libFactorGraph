#include "cnode.h"

/* constructor for factor type node */

CNode::CNode()
{

}

CNode::CNode(FG::NodeType type, int ID, int card, QObject *parent) :
    QObject(parent)
{
    if((type==FG::Variable && ID>0) || (type==FG::Factor && ID<0))
    {
        _ID = ID;
        _type = type;
        _localF.setID(ID);
        _card = card; //only for variable node
    } else qFatal("Invalid Node: Check node type and its ID!");
}

CNode::~CNode()
{
}

void CNode::init(FG::NodeType type, int ID, int card)
{
    //for variable node: its ID must positive and card must be positive nonzero
    if(((type==FG::Variable && ID>0) && (type==FG::Variable && card>0)) || (type==FG::Factor && ID<0))
    {
        _ID = ID;
        _type = type;
        _localF.setID(ID);
        _card = card; //only for variable node
    } else qFatal("Invalid Node: Check node type and its ID!");
}

void CNode::addNeighbor(int nID)
{
    /* check: Factor node must have only variables neighbor and vice versa */
    if((_type==FG::Variable && nID < 0) || (_type==FG::Factor && nID > 0))
        if(!_neighbor.contains(nID))
        {
            _neighbor.append(nID);  //add to the list if it is not yet there
            //QVector<CMessage> emptyMsg;
            //then resize the msgTable
            msgTable.resize(_neighbor.count());
            for(int i=0; i<_neighbor.count(); i++)
            {
                msgTable[i].resize(_neighbor.count());
            }
        }
}

/* _localF only used by Factor-node
   _localF has ID the same as the node itself! */

void CNode::initFactor(QVector<int> vars, QVector<int> cards)
{
    _localF.setScope(vars);
    _localF.setCard(cards);
    int length = getElementProd(cards);
    double initValue = 1.0/double(length);      //uniform distribution
    QVector<double> initJPD(length, initValue);
    _localF.setJPD(initJPD);
}

void CNode::setFactor(QVector<int> vars, QVector<int> cards, QVector<double> vals)
{
    Q_ASSERT_X(vars==_neighbor, "Invalid scope", "Check the order of scope in the neighborhood");
    if(cards.count()!=_neighbor.count())
        qFatal("Invalid scope: Make sure the scope is the same as the neighbor nodes!");
    if(getElementProd(cards)!=vals.count())
        qFatal("Invalid values: Make sure the number of elements reflects the correct cardinality!");
    _localF.setScope(_neighbor);
    _localF.setCard(cards);
    _localF.setJPD(vals);
    if(!checkScopeOrder())
        qFatal("The scope and the neighbor are mismatch!");
}

/*
 *Put a message into the msgTable
 *msgTable is a 2D array, where row = destination, col = sources
 */
bool CNode::putMessage(CMessage msg)
{
    //important check: msg._ID must be the same with the sourceID because the factor product relies on this ID

    //we will put the msg in all rows, means it can be delivered later to all neighbors
    //first get the source index in the list of neighborhood
    int idx = _neighbor.indexOf(msg.getSourceID()); //msg.getID will be the same as msg._ID which is also, apparently, msg._sourceID
    if(idx!=-1) //if it is a valid neighbor, then add the msg to the table
    {
        for(int i=0; i<msgTable.count(); i++)   //scan for all neighbors
        {
            //but skip the source neighbor
            if(i!=idx)  //except from the message source itself :)
            {
                msgTable[i][idx] = msg;
            }
        }
        return true;
    } else return false;
}

/*
 *compute a message as a leaf node: only valid for a factor node
 */
CMessage CNode::computeMessage(bool *success)
{
    CMessage result;
    bool bSuccess;
    if(_localF.isValid())   //isValid means if the _localF has already contained JPD
    {
        result.setCard(_localF.getCard());
        result.setJPD(_localF.getJPD());
        result.setScope(_localF.getScope());
        result.setSourceID(this->getID());    //don't you think this is redundant
        //here we don't specify the varID and it MUST BE assigned outside!
        bSuccess = true;
    } else bSuccess = false;
    if(success!=0) *success = bSuccess; //always check if it is not a NULL pointer!
    return result;
}

CMessage CNode::computeMessage(int destID, bool *success)
{
    CMessage result;    //NOTE: result won't have any ID -> important for product operation
    bool msgValid;
    QVector<int> includedNeighbor = substVector(_neighbor, destID);
    QVector<int> checkNeighbor;
    if(includedNeighbor.isEmpty()) //for example, the case where the node is a leaf node
    {
        qFatal("Invalid neighborhood detected!");
    }
    else
    {
        //for both variable and factor nodes, first we have to do factor product of all neighboring messages
        //first, check the availability of all required input messages
        bool inputMsgsExist = true;
        int col, row;   //of msg_table
        row = _neighbor.indexOf(destID); //get the index of the destID in the neighborhood matrix
        //QVector<CMessage> msgIn(includedNeighbor.count());
        QVector<CMessage> msgIn;
        //for(col=0; col<includedNeighbor.count(); col++)
        for(col=0; col<_neighbor.count(); col++)
        {
            if(col!=row)
            {
                inputMsgsExist = inputMsgsExist & msgTable[row][col].isValid();
                //msgIn[col] = msgTable.at(row).at(col);
                msgIn.append(msgTable.at(row).at(col));
                checkNeighbor.append(msgTable[row][col].getSourceID()); //here, CMessage::_sourceID is important
            }
        }
        if(inputMsgsExist) //if all input messages are present
        {
            if(this->_type==FG::Variable)
            {
                if(checkNeighbor!=includedNeighbor)
                    sortMsg(msgIn, includedNeighbor);
                result = msgIn.at(0);
                for(int i=1; i<msgIn.count(); i++)
                    result.prod(msgIn.at(i));
            }
            else
            {
                //TODO: short the msgIn according to the _localF.getScope()
                result = msgIn.at(0);
                for(int i=1; i<msgIn.count(); i++)
                    result.prod(msgIn.at(i), false);
                result.expand(_localF, destID);
                result.prod(_localF, true);
                //for(int i=0; i<includedNeighbor.count(); i++)
                //    result.summ(includedNeighbor.at(i));
                result.marg(destID);
            }
            msgValid = true;
        }
    }
    if(success!=0) *success = msgValid;
    return result;
}

/*
 *checkScopeOrder will return false, if:
 *- it is a variable node
 *- the number of neighboor and _localF._scopes is mismatch
 */
bool CNode::checkScopeOrder()
{
    if(_type==FG::Variable || _localF.getScope().count()!=_neighbor.count()) return false;
    if(_localF.getScope()==_neighbor) return true;
    else
    {
        _neighbor = _localF.getScope();
        return true;
    }
}

QVector<CMessage> CNode::sortMsg(const QVector<CMessage> &msg, const QVector<int> &nb)
{
    QVector<CMessage> result;
    for(int i=0; i<nb.count(); i++)                 //scan through all neighborhod nb
        for(int j=0; j<msg.count(); j++)            //scan through all message list
            if(msg.at(j).getSourceID()==nb.at(i))
                result.append(msg.at(j));
}
