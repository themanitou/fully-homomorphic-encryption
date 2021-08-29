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

#include <QDataStream>
#include <QString>

#include "calculation/fullhom.h"
#include "communication/commserver.h"
#include "jobmanager.h"

namespace Fhe
{

    JobManager::JobManager() :
        server_(std::make_unique<CommServer>("FullyHomomorphicEncryption"))
    {
        connect(server_.get(), &CommServer::onMessageReceived, this, &JobManager::jobReceived);
        qDebug() << "[JobManager] started and waiting for jobs...";
    }


    JobManager::~JobManager()
    {
        qDebug() << "[JobManager] end.";
    }


    void JobManager::sendResult(QUuid uuid, bool ok, int data)
    {
        QByteArray result;
        QDataStream out(&result, QIODevice::WriteOnly);

        out << uuid << ok << data;
        server_->broadcast(result);
    }


    void JobManager::jobReceived(QByteArray message)
    {
        QDataStream in(message);
        QUuid opId;
        QString jobId;
        in >> opId >> jobId;
        qDebug() << "[jobReceived] opId=" << opId.toString() << ", jobId=" << jobId;

        if (jobId == "Encrypt Bit")
        {
            QUuid bitId;
            int bit;
            in >> bitId >> bit;

            ZZ ciphertext;
            EncryptBit(ciphertext, bit);
            zzMap_[bitId] = ciphertext;

            sendResult(opId, true);
        }
        else if (jobId == "Decrypt Bit")
        {
            QUuid bitId;
            long bit;
            in >> bitId;

            ZZ ciphertext = zzMap_[bitId];
            DecryptBit(bit, ciphertext);

            sendResult(opId, true, bit);
        }
    }


} // namespace Fhe
