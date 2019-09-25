#include "srp.h"

using namespace std;

SRP::SRP(int d, int K, int M)
    : dim(d)
    , K(K)
    , M(M)
    , sigdim(K)
    , p(K * M * d)
{
    assert(d > 0 && K > 0 && M > 0);

    std::normal_distribution<double> normalDistribution(0.);
    std::random_device rd;
    std::default_random_engine rng(rd());
    for (int i = 0; i < K * M * d; i++) {
        p[i] = normalDistribution(rng);
    }
}

SRP::~SRP()
{
}

std::vector<SigType> SRP::getSig(const Scalar* data)
{
    std::vector<SigType> ret(K);
    for (int i = 0; i < K; i++) {
        int sigi = 0;
        for (int j = 0; j < M; j++) {
            int idx = i * M * dim + j * dim;
            float sum = calc_inner_product(dim, &p[idx], data);
            if (sum >= 0) {
                sigi = (sigi << 1) + 1;
            } else {
                sigi = sigi << 1;
            }
        }
        ret[i] = sigi;
    }
    return ret;
}

void SRP::getSig(const Scalar* data, SigType* ret)
{
    for (int i = 0; i < K; i++) {
        int sigi = 0;
        for (int j = 0; j < M; j++) {
            int idx = i * M * dim + j * dim;
            float sum = calc_inner_product(dim, &p[idx], data);
            if (sum >= 0) {
                sigi = (sigi << 1) + 1;
            } else {
                sigi = sigi << 1;
            }
        }
        ret[i] = sigi;
    }
}

SRPCompact::SRPCompact(int d, int K)
    : dim(d)
    , K(K)
    , sigdim((K+63)/64)
    , p(K * d)
{
    assert(d > 0 && K > 0);

    std::normal_distribution<double> normalDistribution(0.);
    std::random_device rd;
    std::default_random_engine rng(rd());
    for (int i = 0; i < K * d; i++) {
        p[i] = normalDistribution(rng);
    }
}

SRPCompact::~SRPCompact()
{
}

std::vector<uint64_t> SRPCompact::getSig(const Scalar* data)
{
    std::vector<uint64_t> ret(sigdim);
    for (int i = 0; i < sigdim; i++) {
        uint64_t sigi = 0;
        for (int j = 0; j < 64; j++) {
            int idx = (i * 64 + j) * dim;
            float sum = calc_inner_product(dim, &p[idx], data);
            if (sum >= 0) {
                sigi = (sigi << 1) + 1;
            } else {
                sigi = sigi << 1;
            }
        }
        ret[i] = sigi;
    }
    return ret;
}

void SRPCompact::getSig(const Scalar* data, uint64_t* ret)
{
    for (int i = 0; i < sigdim; i++) {
        uint64_t sigi = 0;
        for (int j = 0; j < 64; j++) {
            int idx = (i * 64 + j) * dim;
            float sum = calc_inner_product(dim, &p[idx], data);
            if (sum >= 0) {
                sigi = (sigi << 1) + 1;
            } else {
                sigi = sigi << 1;
            }
        }
        ret[i] = sigi;
    }
}