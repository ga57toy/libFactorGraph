#include "cmessage.h"

CMessage::CMessage()
{
}

void CMessage::updateFactor(CMessage newMessage, int msgID)
{
    //setScope(newMessage.getScope());
    //setCard(newMessage.getCard());
    setJPD(newMessage.getJPD());
    setID(msgID);
}

void CMessage::updateFactor(CMessage newMessage)
{
    setJPD(newMessage.getJPD());
}
