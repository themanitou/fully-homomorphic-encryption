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

#ifndef JOBMANAGER_H_
#define JOBMANAGER_H_

#include <memory>

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <QByteArray>
#include <QObject>
#include <QMap>
#include <QUuid>

NTL_CLIENT

namespace Fhe
{

    class CommServer;

    class JobManager : public QObject
    {
        Q_OBJECT

        std::unique_ptr<CommServer> server_;
        QMap<QUuid,ZZ> zzMap_;
        QMap<QUuid,vec_ZZ> vecZZMap_;

    public:
        JobManager();
        ~JobManager() override;

    private:
        void sendResult(QUuid uuid, bool ok, int data = 0);

    private slots:
        void jobReceived(QByteArray jobMessage);
    };

} // namespace Fhe

#endif
