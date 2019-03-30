#pragma once

#include <vector>
#include <array>

template<int D, class Scalar> 
class NDArrayView;

// template<int D, class Scalar>
// struct NDArrayViewGivenD
// {
//     typedef NDArrayView<D, Scalar> ValType;
// };

// template<class Scalar>
// struct NDArrayViewGivenD<0, Scalar>
// {
//     typedef float& ValType;
// };

//a light-weighted compact n-dimensional array
template<int D, class Scalar> 
class NDArray
{
public:
    NDArray(const std::array<size_t, D> &dims):lens(dims), buffer_raw(nullptr){
        elemLens[D] = 1;
        for(int i=D-1;i>=0;--i){
            elemLens[i] = elemLens[i+1]*lens[i];
        }
        buffer = new Scalar[elemLens[0]];
    }

    ~NDArray(){
        delete[] buffer;
        if(buffer_raw!=nullptr){
            delete[] buffer_raw;
        }
    }

    void resize(const std::array<size_t, D> &dims){
        for(int i=0;i<D;i++){
            lens[i] = dims[i];
        }
        elemLens[D] = 1;
        for(int i=D-1;i>=0;--i){
            elemLens[i] = elemLens[i+1]*lens[i];
        }

        delete[] buffer;
        buffer = new Scalar[elemLens[0]];
        if(buffer_raw!=nullptr){
            delete[] buffer_raw;
            buffer_raw = nullptr;
        }
    }

    NDArrayView<D-1, Scalar> operator[](int i){
        return NDArrayView<D-1, Scalar>(buffer+i*elemLens[1], &elemLens[2]);
    } 

    std::array<size_t, D> lens;
    std::array<size_t, D+1> elemLens;

    template<class F, class ...Args>
    void for_enumerate2(const F& f);

    Scalar* begin(){
        return buffer;
    }
    Scalar* end(){
        return buffer+elemLens[0];
    }

    typedef typename NDArrayView<D-1, Scalar>::PType* PType;
    typedef typename NDArrayView<D-1, Scalar>::CPType* CPType;

    PType to_ptr(){
        if(buffer_raw==nullptr){
            int size = 0;
            int curp = 1;
            for(int i=0;i<D-1;i++){
                size += curp*lens[i];
                curp *= lens[i];
            }
            buffer_raw = new void*[size];

            size = 0;
            curp = 1;
            for(int i=0;i<D-2;i++){
                void** basei = buffer_raw+size;
                size += curp*lens[i];
                curp *= lens[i];
                void** baseip1 = buffer_raw+size;

                for(int j=0;basei+j<baseip1;j++){
                    basei[j] = (void*)(baseip1+lens[i+1]*j);
                }
            }

            //for second last dimension
            void** basei = buffer_raw+size;
            void** baseip1 = buffer_raw+size+curp*lens[D-2];
            for(int j=0;basei+j<baseip1;j++){
                basei[j] = (void*)(buffer+lens[D-1]*j);
            }
        }
        return (PType)buffer_raw;
    }

    inline CPType to_cptr(){
        return (CPType) to_ptr();
    }

protected:
    Scalar* buffer;
    void** buffer_raw;       //will be used only when convert to raw_ptr is called
};

template<int D, class Scalar> 
class NDArrayView
{
public:
    NDArrayView<D-1, Scalar> operator[](int i){
        return NDArrayView<D-1, Scalar>(p+(*elemLensp)*i, elemLensp+1);
    }
    NDArrayView(Scalar *p, size_t *elemLensp):p(p), elemLensp(elemLensp){}

    typedef typename NDArrayView<D-1, Scalar>::PType* PType;
    typedef typename NDArrayView<D-1, Scalar>::CPType* CPType;

    Scalar* begin(){
        return p;
    }
    Scalar* end(){
        return p+*(elemLensp-1);
    }
protected:
    Scalar* p;
    size_t* elemLensp;
};

template<class Scalar> 
class NDArrayView<1, Scalar>
{
public:
    Scalar& operator[](int i){
        return p[i];
    }
    NDArrayView(Scalar *p, size_t *elemLensp):p(p), elemLensp(elemLensp){}

    typedef Scalar* PType;
    typedef const Scalar* CPType;
    operator PType() const{
        return p;
    }

    Scalar* begin(){
        return p;
    }
    Scalar* end(){
        return p+*(elemLensp-1);
    }
protected:
    Scalar* p;
    size_t* elemLensp;
};

template<class Scalar> 
class NDArrayView<0, Scalar>
{
public:
    typedef Scalar PType;
    typedef const Scalar CPType;
};

NDArray<1, double> make_darray(size_t sz0);
NDArray<2, double> make_darray(size_t sz0, size_t sz1);
NDArray<3, double> make_darray(size_t sz0, size_t sz1, size_t sz2);
NDArray<4, double> make_darray(size_t sz0, size_t sz1, size_t sz2, size_t sz3);

NDArray<1, float> make_farray(size_t sz0);
NDArray<2, float> make_farray(size_t sz0, size_t sz1);
NDArray<3, float> make_farray(size_t sz0, size_t sz1, size_t sz2);
NDArray<4, float> make_farray(size_t sz0, size_t sz1, size_t sz2, size_t sz3);

NDArray<1, int> make_iarray(size_t sz0);
NDArray<2, int> make_iarray(size_t sz0, size_t sz1);
NDArray<3, int> make_iarray(size_t sz0, size_t sz1, size_t sz2);
NDArray<4, int> make_iarray(size_t sz0, size_t sz1, size_t sz2, size_t sz3);

