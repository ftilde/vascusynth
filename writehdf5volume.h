#include "H5Cpp.h"
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

template<class T>
struct TemplateVec3 {
    T x, y, z;

    TemplateVec3<T>(T x, T y, T z) : x(x), y(y), z(z) {}
    bool operator==(TemplateVec3<T> other) {
        return x == other.x && y == other.y && z == other.z;
    }

    TemplateVec3<T> operator*(TemplateVec3<T> other) {
        return { x * other.x, y * other.y, z * other.z};
    }
    TemplateVec3<T> operator*(T other) {
        return { x * other, y * other, z * other};
    }
    TemplateVec3<T> operator-(TemplateVec3<T> other) {
        return { x - other.x, y - other.y, z - other.z};
    }
    TemplateVec3<T> operator-(T other) {
        return { x - other, y - other, z - other};
    }
    TemplateVec3<T> operator+(TemplateVec3<T> other) {
        return { x + other.x, y + other.y, z + other.z};
    }
    TemplateVec3<T> operator+(T other) {
        return { x + other, y + other, z + other};
    }
    TemplateVec3<T> operator/(TemplateVec3<T> other) {
        return { x / other.x, y / other.y, z / other.z};
    }
    TemplateVec3<T> operator/(T other) {
        return { x / other, y / other, z / other};
    }

    TemplateVec3<T> max(TemplateVec3<T> other) {
        return { std::max(x, other.x), std::max(y, other.y), std::max(z, other.z), };
    }

    TemplateVec3<T> min(TemplateVec3<T> other) {
        return { std::min(x, other.x), std::min(y, other.y), std::min(z, other.z), };
    }

    TemplateVec3<T> ceil() {
        return { std::ceil(x), std::ceil(y), std::ceil(z), };
    }

    TemplateVec3<T> floor() {
        return { std::floor(x), std::floor(y), std::floor(z), };
    }
    T hsum() {
        return { x + y + z };
    }
    T hmul() {
        return { x * y * z };
    }

    template<typename O>
    TemplateVec3<T>(TemplateVec3<O> other)
        : x(other.x)
        , y(other.y)
        , z(other.z)
    {
    }

    TemplateVec3<T>(T* vals)
        : x(vals[0])
        , y(vals[1])
        , z(vals[2])
    {
    }
    TemplateVec3<T>(T v)
        : x(v)
        , y(v)
        , z(v)
    {
    }
};

typedef TemplateVec3<size_t> svec3;
typedef TemplateVec3<int> ivec3;
typedef TemplateVec3<float> vec3;
typedef TemplateVec3<double> dvec3;

class HDF5WriteException : public std::exception {
public:
    HDF5WriteException(const std::string& what = "") : what_(what) {}
    virtual ~HDF5WriteException() throw() {}

    virtual const char* what() const throw() {
        return what_.c_str();
    }
protected:
    std::string what_;
};

#if H5_VERSION_GE(1, 10, 1) || H5_VERSION_GE(1, 8, 21) && !H5_VERSION_GE(1, 10, 0)
#define H5_STUPID_LOCATION_API_CHANGES
#endif

static void writeVec3Attribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const vec3& vec) {
    // Convert tgt vector to hdf style
    float h5vec[3];
    h5vec[2] = vec.x;
    h5vec[1] = vec.y;
    h5vec[0] = vec.z;

    H5::DataType h5type = H5::PredType::IEEE_F32LE;

    // Create Attribute
    hsize_t numberOfValues = 3;
    const H5::DataSpace dataspace(1, &numberOfValues);
    H5::Attribute attribute = loc.createAttribute(name, h5type, dataspace);

    attribute.write(h5type, h5vec);
}


struct ByteVolume {
    ByteVolume(svec3 dim)
        : dim(dim)
        , data(dim.x*dim.y*dim.z)
    { }
    svec3 dim;
    std::vector<unsigned char> data;

    size_t data_index(svec3 pos) {
        return pos.x + dim.x * (pos.y + pos.z*dim.y);
    }

    void set_voxel(svec3 pos, unsigned char val) {
        data.at(data_index(pos)) = val;
    }
    char get_voxel(svec3 pos) {
        return data.at(data_index(pos));
    }
};

struct HDF5Volume {
    static HDF5Volume create(const std::string& fileName, svec3 size, float voxelsize, int deflateLevel=1);
    void writeSlice(const ByteVolume& vol, size_t z);

    std::unique_ptr<H5::DataSet> dataSet;
    svec3 size;
};
H5::DataType HDF5_VOL_TYPE = H5::PredType::STD_U8LE;

HDF5Volume HDF5Volume::create(const std::string& fileName, svec3 size, float voxelsize, int deflateLevel) {
    bool shuffle = false;
    svec3 chunkSize = svec3(0,0,0);
    bool truncateFile = true;

    hsize_t dimensions_hdf5[3];
    dimensions_hdf5[2] = size.x;
    dimensions_hdf5[1] = size.y;
    dimensions_hdf5[0] = size.z;
    // Open the file ======================================================================================
    std::unique_ptr<H5::H5File> file = nullptr;

    try {
        // Flags for opening/creating the file
        // Although the documentation states that the constructor of H5File "opens or creates" the file
        // it only decides what to to based on the flags given as second argument...
        // Because of that if we explicitly want to truncate the file, the file does not exist yet or is
        // not a hdf5 file, we set the truncate flag.
        unsigned int flags = 0;
        if(truncateFile) {
            flags |= H5F_ACC_TRUNC;
        } else {
            // If we do not truncate the file, we have to specify that we open it in Read/Write mode.
            flags |= H5F_ACC_RDWR;
        }

        // Setup file access property list for file creation.
        H5::FileAccPropList accList(H5::FileAccPropList::DEFAULT);

        // We require at least format version 1.8 to be able to write arbitrarily large attribute data.
        // Wrapper not available in 1.8.13, so we use the c version:
        //accList.setLibverBounds(H5F_LIBVER_18, H5F_LIBVER_LATEST);
        // ...
        // Unfortunately libhdf5 v.1.10 removed H5F_LIBVER_18 (which for 1.8 was an alias for H5F_LIBVER_LATEST)
        // so for now we set the minimum version to H5F_LIBVER_LATEST and hope that there are not incompatiblitities
        // => TODO: Set the min version accordingly if there is (backwarts compatible) support in a future version
        H5Pset_libver_bounds(accList.getId(), H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

        // Create the hdf5 file using the propertylists created earlier.
        file = std::unique_ptr<H5::H5File>(new H5::H5File(fileName, flags, H5::FileCreatPropList::DEFAULT, accList));
    } catch(H5::Exception error) {
        throw HDF5WriteException("Error opening/creating hdf5 file:" + error.getFuncName() + ": " + error.getDetailMsg());
    }


    // Create the dataSet =================================================================================
    std::unique_ptr<H5::DataSet> dataSet = nullptr;

    try {
        // Prepare to create the volume inside the file
        H5::DataSpace dataSpace(3, dimensions_hdf5);


        // Set up properties for the volume data set to be created.
        H5::DSetCreatPropList propList;

        // Set Shuffle if wanted.
        if(shuffle) {
            propList.setShuffle();
        }

        // Set chunking
        // We need to set chunking even if we want to use the default chunk size
        // because otherwise the chunk size of propList will be uninitiated
        // apparently.
        // We also cannot create a copy of the default property list because the
        // HDF5 lib does not really create a copy (although it states that it does),
        // but still references the original interal c-style property list...
        // Because of that we:

        // ... Set the chunk size to be an xy-slice of the volume manually, if
        // none was specified.
        if(chunkSize == svec3(0,0,0)) {
            chunkSize.x = size.x;
            chunkSize.y = size.y;
            chunkSize.z = 1;
        }
        // ... And set the chunk size of the previously created property list.
        hsize_t chunkdimHDF5[3];
        // If we have a multidimensionale Volume, we definitely want the channels separated
        chunkdimHDF5[2] = chunkSize.x;
        chunkdimHDF5[1] = chunkSize.y;
        chunkdimHDF5[0] = chunkSize.z;
        propList.setChunk(3, chunkdimHDF5);

        // Set compression level
        // setDeflate(0) still adds a filter altough it does not compress.
        if(deflateLevel != 0) {
            propList.setDeflate(deflateLevel);
        }

        // Finally try to create the data set inside the file
        dataSet = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet("/vol", HDF5_VOL_TYPE, dataSpace, propList)));
    } catch(H5::Exception error) {
        throw HDF5WriteException("Error constructing dataset:" + error.getFuncName() + ": " + error.getDetailMsg());
    }
    //Write metadata to dataSet
    float voxelSizeUm = voxelsize*1000;
    writeVec3Attribute(*dataSet, "element_size_um", vec3(voxelSizeUm, voxelSizeUm, voxelSizeUm));

    return HDF5Volume {
        std::move(dataSet),
        size,
    };
}

void HDF5Volume::writeSlice(const ByteVolume& vol, size_t z) {

    //Write data to dataset
    try {
        hsize_t start[3];
        start[2] = 0;
        start[1] = 0;
        start[0] = z;

        hsize_t dimensions_hdf5[3];
        dimensions_hdf5[2] = vol.dim.x;
        dimensions_hdf5[1] = vol.dim.y;
        dimensions_hdf5[0] = vol.dim.z;

        H5::DataSpace fileSpace(dataSet->getSpace());
        fileSpace.selectHyperslab(H5S_SELECT_SET, dimensions_hdf5, start);

        //Memory space does not have an offset
        H5::DataSpace memSpace(3, dimensions_hdf5);

        //Write the volume to disk.
        dataSet->write(vol.data.data(), HDF5_VOL_TYPE, memSpace, fileSpace);

    } catch(H5::Exception error) { // catch HDF5 exceptions
        std::cout << error.getFuncName() << ": " << error.getDetailMsg() << std::endl;
        throw HDF5WriteException("An Error occured while writing volume to file");
    }
}

