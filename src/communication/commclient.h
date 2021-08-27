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

#ifndef COMMCLIENT_H_
#define COMMCLIENT_H_

#include <memory>

#include <QLocalSocket>
#include <QObject>

namespace Fhe
{

class CommClient : public QObject
{
    Q_OBJECT

    CommClient() = delete;
    CommClient(const CommClient&) = delete;
    CommClient& operator=(const CommClient&) = delete;

public:
    CommClient(const QString& wsLabel);
    bool sendMessage(const QByteArray& message);
    bool flush();

signals:
    void onMessageReceived(QByteArray message);

private slots:
    void onError(QLocalSocket::LocalSocketError socketError);
    void onRead();

private:
    std::unique_ptr<QLocalSocket> socket_;
    quint32 blockSize_;
};

} // namespace Fhe
#endif
