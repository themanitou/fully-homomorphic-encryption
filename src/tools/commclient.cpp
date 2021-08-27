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

#include "commclient.h"

#include <QDataStream>
#include <QtNetwork>

namespace Fhe
{

CommClient::CommClient(const QString& wsLabel) :
    blockSize_(0)
{
    socket_ = std::make_unique<QLocalSocket>(this);
    connect(socket_.get(), SIGNAL(error(QLocalSocket::LocalSocketError)), this, SLOT(onError(QLocalSocket::LocalSocketError)));
    connect(socket_.get(), &QLocalSocket::readyRead, this, &CommClient::onRead);
    socket_->connectToServer(wsLabel);
}


///
///
///
bool CommClient::sendMessage(const QByteArray& message)
{
    QByteArray block;
    QDataStream out(&block, QIODevice::WriteOnly);
    quint32 blockSize = message.size();

    out.setVersion(QDataStream::Qt_4_0);
    out << blockSize;
    out << message;

    return (socket_->write(block) != -1);
}

///
///
///
bool CommClient::flush()
{
    return socket_->flush();
}

///
///
///
void CommClient::onError(QLocalSocket::LocalSocketError socketError)
{
    qDebug() << QString("[CommClient] LocalSocketError: %1 (errorCode=%2)").arg(socket_->errorString()).arg(socketError);
}


///
///
///
void CommClient::onRead()
{
    QDataStream in(socket_.get());
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
                qDebug() << QString("[CommClient::onRead] received a possibly corrupted message, blockSize_ = 0");
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
            emit onMessageReceived(message);
        }

        blockSize_ = 0;
    }
}

} // namespace Fhe
