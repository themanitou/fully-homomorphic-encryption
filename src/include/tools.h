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


#ifndef TOOLS_H_
#define TOOLS_H_


#include <string>


namespace Fhe
{
    using namespace std;

    class Tools
    {
    public:
        Tools(const Tools&) = delete;
        Tools& operator=(const Tools&) = delete;

        ~Tools() = default;

        static Tools& getInstance();
        string convBase (const unsigned long&, const long&);

    private:
        Tools() = default;
    };
}

#endif /* TOOLS_H_ */
