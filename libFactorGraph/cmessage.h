#ifndef CMESSAGE_H
#define CMESSAGE_H
#include <QString>
#include "cfactor.h"

class CMessage : public CFactor
{
public:
    CMessage();                                 //the most basic form of a message is without ID
    void setVarID(int varID);
    void setName(QString strname) {_name = strname;}  //Optional, might be useful in addition to msgID
    void setSourceID(int src) {_sourceID=src;}
    QString getName() const {return _name;}
    int getSourceID() const {return _sourceID;}
    void updateFactor(CMessage newMessage, int msgID);
    void updateFactor(CMessage newMessage);
private:
    int _sourceID;
    QString _name;
};

#endif // CMESSAGE_H
