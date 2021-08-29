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

#ifndef COMMSERVER_H_
#define COMMSERVER_H_

#include <memory>

#include <QByteArray>
#include <QLocalSocket>
#include <QObject>

class QLocalServer;

namespace Fhe
{

class CommServer : public QObject
{
    Q_OBJECT

    friend class CommServerSocket;

    CommServer() = delete;
    CommServer(const CommServer&) = delete;
    CommServer& operator=(const CommServer&) = delete;

public:
    CommServer(const QString& commName);
    ~CommServer() override;

    bool broadcast(const QByteArray& message);

signals:
    void onMessageReceived(QByteArray message);

private slots:
    void onConnect();

private:
    std::unique_ptr<QLocalServer> server_;
};


class CommServerSocket : public QObject
{
    Q_OBJECT

    CommServer& server_;
    QLocalSocket* socket_;
    quint32 blockSize_;

    CommServerSocket() = delete;
    CommServerSocket(const CommServerSocket&) = delete;
    CommServerSocket& operator=(const CommServerSocket&) = delete;

public:
    CommServerSocket(CommServer& server, QLocalSocket* socket);

    bool sendMessage(const QByteArray& message);

private slots:
    void onRead();
    void onError(QLocalSocket::LocalSocketError socketError);
};

} // namespace Fhe
#endif
