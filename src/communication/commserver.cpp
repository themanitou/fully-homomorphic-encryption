/******************************************************************
 *
 *   Fully-Homomorphic Cryptography library,
 *   based on Gentry-Halevi ideal lattice scheme.
 *
 *   Author: Quan Nguyen (https://github.com/themanitou)
 *
 *   This library is open-source software distributed under the
 *   terms of the GNU Lesser General Public License (LGPL) version
 *   2.1 or later.  See the file doc/copying.txt for complete
 *   details on the licensing of this library.
 *
 *******************************************************************/

#include "commserver.h"

#include <QDataStream>
#include <QLocalServer>
#include <QtNetwork>

namespace Fhe
{

CommServerSocket::CommServerSocket(CommServer& server, QLocalSocket* socket) :
    QObject(&server),
    server_(server),
    socket_(socket),
    blockSize_(0)
{
    connect(socket_, &QLocalSocket::readyRead, this, &CommServerSocket::onRead);
    connect(socket_, static_cast<void (QLocalSocket::*)(QLocalSocket::LocalSocketError)>(&QLocalSocket::error), this, &CommServerSocket::onError);
}


///
///
///
void CommServerSocket::onRead()
{
    QDataStream in(socket_);
    in.setVersion(QDataStream::Qt_4_0);

    while (socket_->bytesAvailable() && blockSize_ <= socket_->bytesAvailable())
    {
        if (blockSize_ == 0)
        {
            // this is the first read of this socket, check how many bytes the server should expect
            if (socket_->bytesAvailable() < static_cast<int>(sizeof(quint32)))
            {
                return;
            }

            in >> blockSize_;
            if (blockSize_ == 0)
            {
                qDebug() << QString("[CommServerSocket::onRead] received a possibly corrupted message, blockSize_ = 0");
                return;
            }
        }

        // just wait (return) in case not all data has arrived yet
        if (socket_->bytesAvailable() < blockSize_ || in.atEnd())
        {
            return;
        }

        QByteArray message;
        in >> message;

        if (message.size())
        {
            emit server_.onMessageReceived(message);
        }

        // reset the blockSize_ to get ready for the next message block
        blockSize_ = 0;
    }
}


///
///
///
void CommServerSocket::onError(QLocalSocket::LocalSocketError socketError)
{
    qDebug() << QString("[CommServerSocket] LocalSocketError: %1 (errorCode=%2)").arg(socket_->errorString()).arg(socketError);
}


///
///
///
bool CommServerSocket::sendMessage(const QByteArray& message)
{
    QByteArray block;
    QDataStream out(&block, QIODevice::WriteOnly);
    quint32 blockSize = message.size();

    out.setVersion(QDataStream::Qt_4_0);
    out << blockSize;
    out << message;

    return (socket_->write(block) != -1);
}

//***********************************


///
///
///
CommServer::CommServer(const QString& commName)
{
    server_ = std::make_unique<QLocalServer>(this);
    connect(server_.get(), &QLocalServer::newConnection, this, &CommServer::onConnect);

    if (!server_->listen(commName))
    {
        qDebug() << QString("[CommServer::CommServer] listen failed with name %1: %2").arg(commName).arg(server_->errorString());
        return;
    }
    qDebug() << QString("[CommServer::CommServer] listening with name %1").arg(commName);
}


///
///
///
CommServer::~CommServer() = default;


///
///
///
void CommServer::onConnect()
{
    QLocalSocket* rawSocket_ = server_->nextPendingConnection();
    auto* socket = new CommServerSocket(*this, rawSocket_);
    connect(rawSocket_, SIGNAL(disconnected()), socket, SLOT(deleteLater()));
}


///
/// \brief send message to all connected clients
///
bool CommServer::broadcast(const QByteArray& message)
{
    bool allSuccessful = true;
    for (auto socket : findChildren<CommServerSocket*>())
    {
        if (!socket->sendMessage(message))
        {
            allSuccessful = false;
        }
    }
    return allSuccessful;
}


} // namespace Fhe
