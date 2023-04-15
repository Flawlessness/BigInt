#include "BigInt.h"

int32_t main()
{
    BigInt a("-623492");
    BigInt b("-859038539875934578347");
    BigInt c = 99;
    a+=b;
    std :: cout << a;
    return 0;
}
