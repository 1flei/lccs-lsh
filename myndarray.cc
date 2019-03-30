#include "myndarray.h"


NDArray<1, double> make_darray(size_t sz0)
{
    return NDArray<1, double>({sz0});
}
NDArray<2, double> make_darray(size_t sz0, size_t sz1)
{
    return NDArray<2, double>({sz0, sz1});
}
NDArray<3, double> make_darray(size_t sz0, size_t sz1, size_t sz2)
{
    return NDArray<3, double>({sz0, sz1, sz2});
}
NDArray<4, double> make_darray(size_t sz0, size_t sz1, size_t sz2, size_t sz3)
{
    return NDArray<4, double>({sz0, sz1, sz2, sz3});
}

NDArray<1, float> make_farray(size_t sz0)
{
    return NDArray<1, float>({sz0});
}
NDArray<2, float> make_farray(size_t sz0, size_t sz1)
{
    return NDArray<2, float>({sz0, sz1});
}
NDArray<3, float> make_farray(size_t sz0, size_t sz1, size_t sz2)
{
    return NDArray<3, float>({sz0, sz1, sz2});
}
NDArray<4, float> make_farray(size_t sz0, size_t sz1, size_t sz2, size_t sz3)
{
    return NDArray<4, float>({sz0, sz1, sz2, sz3});
}

NDArray<1, int> make_iarray(size_t sz0)
{
    return NDArray<1, int>({sz0});
}
NDArray<2, int> make_iarray(size_t sz0, size_t sz1)
{
    return NDArray<2, int>({sz0, sz1});
}
NDArray<3, int> make_iarray(size_t sz0, size_t sz1, size_t sz2)
{
    return NDArray<3, int>({sz0, sz1, sz2});
}
NDArray<4, int> make_iarray(size_t sz0, size_t sz1, size_t sz2, size_t sz3)
{
    return NDArray<4, int>({sz0, sz1, sz2, sz3});
}