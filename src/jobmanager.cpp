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


    void JobManager::sendResult(QUuid uuid, bool ok, QString data)
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
            long ret = EncryptBit(ciphertext, bit);
            zzMap_[bitId] = ciphertext;

            sendResult(opId, ret == 0);
        }
        else if (jobId == "Decrypt Bit")
        {
            QUuid bitId;
            long bit;
            in >> bitId;

            ZZ ciphertext = zzMap_[bitId];
            long ret = DecryptBit(bit, ciphertext);

            sendResult(opId, ret == 0, QString::number(bit));
        }
        else if (jobId == "Encrypt Byte")
        {
            QUuid byteId;
            int byte;
            in >> byteId >> byte;

            vec_ZZ ciphertext;
            long ret = EncryptByte(ciphertext, byte);
            vecZZMap_[byteId] = ciphertext;

            sendResult(opId, ret == 0);
        }
        else if (jobId == "Decrypt Byte")
        {
            QUuid byteId;
            long byte;
            in >> byteId;

            vec_ZZ ciphertext = vecZZMap_[byteId];
            long ret = DecryptByte(byte, ciphertext);

            sendResult(opId, ret == 0, QString::number(byte));
        }
        else if (jobId == "And Bit")
        {
            QUuid inBitId1, inBitId2;
            QUuid outBitId;
            in >> inBitId1 >> inBitId2 >> outBitId;

            ZZ andBit;
            long ret = AndBit(andBit, zzMap_[inBitId1], zzMap_[inBitId2]);
            zzMap_[outBitId] = andBit;

            sendResult(opId, ret == 0);
        }
        else if (jobId == "Xor Bit")
        {
            QUuid inBitId1, inBitId2;
            QUuid outBitId;
            in >> inBitId1 >> inBitId2 >> outBitId;

            ZZ xorBit;
            long ret = XorBit(xorBit, zzMap_[inBitId1], zzMap_[inBitId2]);
            zzMap_[outBitId] = xorBit;

            sendResult(opId, ret == 0);
        }
        else if (jobId == "Flip Bit")
        {
            QUuid inBitId;
            QUuid outBitId;
            in >> inBitId >> outBitId;

            ZZ flipBit;
            long ret = FlipBit(flipBit, zzMap_[inBitId]);
            zzMap_[outBitId] = flipBit;

            sendResult(opId, ret == 0);
        }
        else if (jobId == "Choose Bit")
        {
            QUuid inBitIdChoice, inBitId1, inBitId2;
            QUuid outBitId;
            in >> inBitIdChoice >> inBitId1 >> inBitId2 >> outBitId;

            ZZ choosenBit;
            long ret = ChooseBit(choosenBit, zzMap_[inBitIdChoice], zzMap_[inBitId1], zzMap_[inBitId2]);
            zzMap_[outBitId] = choosenBit;

            sendResult(opId, ret == 0);
        }
    }


} // namespace Fhe
