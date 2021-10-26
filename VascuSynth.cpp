/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: VascuSynth.cpp,v $
    Language: C++
    Date: $Date: 2011/02/08 10:43:00 $
    Version: $Revision: 1.0 $

Copyright (c) 2011 Medical Imaging Analysis Lab, Simon Fraser University,
British Columbia, Canada.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

 * The name of the Insight Consortium, nor the names of any consortium members,
nor of any contributors, may be used to endorse or promote products derived
from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
misrepresented as being the original software.

 * Free for non-commercial use only.  For commercial use, explicit approval
must be requested by contacting the Authors.

 * If you use the code in your work, you must acknowledge it

 * Modifications of the source code must also be released as open source

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

//commands for mkdir
#ifdef _WIN32
#include "direct.h"
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "writehdf5volume.h"
#include "intervalwalker.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <iterator>
#include <cmath>
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

/**
 * Utility function to read a text file and store the lines into a vector
 * that is returned.  This way we do not have to have a bunch of files open at
 * the same time.
 * @param const char* filename the filename to read
 * @return vector<string> a vector of file lines from the file
 * @throws string exception if the file cannot be read
 */
vector<string> * readFileLines(const char * filename){

    ifstream oFile;
    oFile.open(filename, ios::in);
    vector<string> * lines = new vector<string>;
    string line;

    if(oFile.is_open()){

        while(!oFile.eof()){
            getline(oFile, line);
            if(line.empty()) {
                break;
            }
            lines->push_back(line);
        }

        oFile.close();

    } else {
        throw "Could not open file: " + ( (string) filename);
    }

    return lines;

}


#include "OxygenationMap.h"
#include "SupplyMap.h"
#include "TreeDrawer.h"
#include "VascularTree.h"



/**
 * make dir that is cross platform
 */
int mmkdir(const char * dirname) {
#ifdef _WIN32
    return mkdir(dirname);
#else
    return mkdir(dirname, 0777);
#endif
}

/**
 * itoa is non standard so define it and use it,
 * converts an integer to a string
 */
string itoa(int value, int base) {

    string buf;

    // check that the base if valid
    if (base < 2 || base > 16) return buf;

    enum { kMaxDigits = 35 };
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.

    int quotient = value;

    // Translating number to string with base:
    do {
        buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
        quotient /= base;
    } while ( quotient );

    // Append the negative sign
    if ( value < 0) buf += '-';

    reverse( buf.begin(), buf.end() );
    return buf;
}

/**
 *  Reads the parameters from the parameter file and then builds
 *  the vascular structure in the form of a tree
 *
 *    Parameter File Entries:
 *
 *    SUPPLY_MAP: supply map file
 *    OXYGENATION_MAP: oxygenation map file
 *    PERF_POINT: perf_x perf_y perf_z
 *    PERF_PRESSURE: perfussion pressure
 *    TERM_PRESSURE: termination pressure
 *    PERF_FLOW:    perfusion flow
 *    RHO: rho
 *    GAMMA: gamma
 *    LAMBDA: lambda
 *    MU: mu
 *    MIN_DISTANCE: minDistance
 *    NUM_NODES: numNodes
 *    VOXEL_WIDTH: voxelWidth
 *    CLOSEST_NEIGHBOURS: closestNeighbours
 */
VascularTree * buildTree(const char * filename){

    SupplyMap * sm = NULL;
    OxygenationMap * om = NULL;

    double* perf = new double[3];
    double pperf;
    double pterm;
    double qperf;
    double rho;
    double gamma;
    double lambda;
    double mu;
    double minDistance;
    int numNodes;
    double voxelWidth;
    int closestNeighbours;
    int randomSeed = -1;
    string line;
    string supplyMapFileName;
    string oxygenMapFileName;

    //c++ doesn't allow us to check for undefined variables
    //so we need a boolean so that we know that the variables
    //are defined.
    bool perfSet = false;
    bool pperfSet = false;
    bool ptermSet = false;
    bool qperfSet = false;
    bool rhoSet = false;
    bool gammaSet = false;
    bool lambdaSet = false;
    bool muSet = false;
    bool minDistanceSet = false;
    bool numNodesSet = false;
    bool voxelWidthSet = false;
    bool closestNeighboursSet = false;
    bool supplyMapFileNameSet = false;
    bool oxygenMapFileNameSet = false;

    vector<string> *mapFilesLines = readFileLines(filename);
    int size = (int) mapFilesLines->size();

    for (int i=0; i < size; i++) {

        line = mapFilesLines->at(i);

        if (line.compare("") == 0) {
            break;
        }

        int colonPosition = line.find(":");
        string name = line.substr(0, colonPosition);
        string value = line.substr(colonPosition+2);

        if (name.compare("SUPPLY_MAP") == 0) {

            //store the supplyMapFileName for later
            //the ordering of things matters so just save the file name
            //then initialize the map later so that the user can type the fields
            //in whatever order
            supplyMapFileName = value;
            supplyMapFileNameSet = true;

        } else if (name.compare("OXYGENATION_MAP") == 0) {

            //same as supply map since we need random seed before
            //we init the oxygenation map
            oxygenMapFileName = value;
            oxygenMapFileNameSet = true;

        } else if (name.compare("PERF_POINT") == 0) {

            int spacePosition = value.find(" ");
            string pointvalue = value.substr(0, spacePosition);
            perf[0] = atof(pointvalue.c_str());

            int spacePosition2 = value.find(" ", spacePosition+1);
            pointvalue = value.substr(spacePosition+1, spacePosition2);
            perf[1] = atof(pointvalue.c_str());

            pointvalue = value.substr(spacePosition2+1);
            perf[2] = atof(pointvalue.c_str());

            perfSet = true;

        } else if (name.compare("PERF_PRESSURE") == 0){

            pperf = atof(value.c_str());
            pperfSet = true;

        } else if (name.compare("TERM_PRESSURE") == 0){

            pterm = atof(value.c_str());
            ptermSet = true;

        } else if (name.compare("PERF_FLOW") == 0){

            qperf = atof(value.c_str());
            qperfSet = true;

        } else if (name.compare("RHO") == 0){

            rho = atof(value.c_str());
            rhoSet = true;

        } else if (name.compare("GAMMA") == 0){

            gamma = atof(value.c_str());
            gammaSet = true;

        } else if (name.compare("LAMBDA") == 0){

            lambda = atof(value.c_str());
            lambdaSet = true;

        } else if (name.compare( "MU") == 0){

            mu = atof(value.c_str());
            muSet = true;

        } else if (name.compare("MIN_DISTANCE") == 0){

            minDistance = atof(value.c_str());
            minDistanceSet = true;

        } else if (name.compare("NUM_NODES") == 0){

            numNodes = atoi(value.c_str());
            numNodesSet = true;

        } else if (name.compare("VOXEL_WIDTH") == 0){

            voxelWidth = atof(value.c_str());
            voxelWidthSet = true;

        } else if (name.compare("CLOSEST_NEIGHBOURS") == 0){

            closestNeighbours = atoi(value.c_str());
            closestNeighboursSet = true;

        } else if (name.compare("RANDOM_SEED") == 0){

            randomSeed = atoi(value.c_str());

        } else {

        }


    }

    //make sure that we have everything defined
    if (perfSet && pperfSet && ptermSet && qperfSet && rhoSet && gammaSet && lambdaSet && muSet && minDistanceSet && numNodesSet && voxelWidthSet && closestNeighboursSet && supplyMapFileNameSet && oxygenMapFileNameSet) {

        //load the supply map
        sm = new SupplyMap();

        try {
            sm->loadMap(supplyMapFileName);
        } catch (char * str) {
            throw (string) str;
        }

        //load the oxygenation map, rand seed will be -1 if there
        //is no randomSeed specified or it will be whatever the user specifies
        om = new OxygenationMap(sm, randomSeed);

        try {
            om->loadMap(oxygenMapFileName);
        } catch (char * str) {
            throw (string) str;
        }
        om->supply = sm;

        //TODO: should probably check and make sure everything is defined
        //and throw an error if it is not (and catch the error in the main function)
        VascularTree *vt = new VascularTree(om, perf, pperf, pterm, qperf, rho, gamma, lambda, mu, minDistance, numNodes, voxelWidth, closestNeighbours);
        vt->buildTree();

        return vt;

    } else {

        string errorStr = "Error while parsing the parameter file, not all parameters have been defined.";
        throw errorStr;

    }
}

/**
 * draws the tree to a matrix
 */
TreeDrawer *drawTree(VascularTree * vt, double * corner1, double * corner2, double imageVoxelWidth){
    TreeDrawer * td = new TreeDrawer(vt, imageVoxelWidth, corner1, corner2);
    td->drawImage();
    return td;
}

static bool inSegment(vec3 point, vec3 p1, vec3 p2, float radius){
    vec3 diff = p1-point;
    vec3 diff2 = diff*diff;
    float d2 = diff2.hsum();
    if(d2 < radius*radius) {
        return true;
    }

    vec3 pdiff = p2-p1;
    float t = -(pdiff*diff).hsum()/(pdiff*pdiff).hsum();

    if (t < 0 || t > 1) {
        return false;
    } else {
        vec3 p = (pdiff*t+p1)-point;
        return (p*p).hsum() < radius*radius;
    }
}

struct NoNoise {
    unsigned char apply(unsigned char val) {
        return val;
    }
    const static int deflateLevel = 1;
};

struct GaussianNoise {
    double median;
    double sigma;
    MTRand rand;
    unsigned char apply(unsigned char c) {
        int x = c + (int)(rand.randNorm(median, sigma));
        if(x > 255) {
            return (unsigned char)255;
        } else if(x < 0) {
            return (unsigned char)0;
        } else {
            return (unsigned char) x;
        }
    }
    const static int deflateLevel = 0;
};

/**
 * draws the tree from the TreeDrawer into a volumetric 3D
 * image as a series of 2D png slices
 */
template<typename Noise>
void drawImage(VascularTree& td, svec3 mapSize, svec3 size, float voxelWidth, const char* rootName, Noise noise){

    // VascularTree components are apparently in "physical" coordinates (with
    // applied spacing, but no offset). So in order to go from normalized
    // ([0,1]) to physical we go normalized->voxel->physical.
    vec3 normalizedToVT = vec3(mapSize);

    vec3 sampleToNormalized = vec3(1.0) / vec3(size);

    vec3 sampleToVT = normalizedToVT * sampleToNormalized;
    vec3 vtToSample = vec3(1.0) / sampleToVT;

    auto vol = HDF5Volume::create(std::string(rootName) + ".h5", size, voxelWidth, Noise::deflateLevel);

    NodeTable& nt = td.nt;
    size_t numNodes = nt.nodes.size();

    std::vector<Interval<int, size_t>> cylinderIntervals;

    struct Cylinder {
        vec3 p1;
        vec3 p2;
        float radius;
        svec3 llf;
        svec3 urb;
    };
    std::vector<Cylinder> cylinders;
    for(int i = 1; i < numNodes; i++){
        dvec3 p1(nt.getPos(i));
        int parent = nt.getParent(i);
        dvec3 p2(nt.getPos(parent));
        float radius = nt.getRadius(i)/td.mapVoxelWidth;

        svec3 llf = (vtToSample*(p1.min(p2) - vec3(radius))).floor().max(vec3(0.0));
        svec3 urb = (vtToSample*(p1.max(p2) + vec3(radius))).ceil();

        cylinders.push_back(Cylinder {
            p1,
            p2,
            radius,
            llf,
            urb,
        });
    }

    size_t numCylinders = cylinders.size();
    std::vector<Interval<int, size_t>> zIntervals;
    for(int i = 0; i < numCylinders; ++i) {
        auto& c = cylinders[i];
        zIntervals.emplace_back(
                    c.llf.z,
                    c.urb.z,
                    i);
    }
    IntervalWalker<int, size_t> zWalker(0, std::move(zIntervals));

    svec3 sliceDim { size.x, size.y, 1 };
    for(int z = 0; z < size.z; z++){
        ByteVolume slice(sliceDim);

        auto cylinderIt = zWalker.next();
        std::vector<Interval<int, size_t>> yIntervals;
        for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
            auto& i = it->value;
            auto& c = cylinders[i];
            yIntervals.emplace_back(c.llf.y, c.urb.y, i);
        }
        IntervalWalker<int, size_t> yWalker(0, std::move(yIntervals));

        for(int y = 0; y < size.y; y++){

            auto cylinderIt = yWalker.next();
            std::vector<Interval<int, size_t>> xIntervals;
            for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
                auto& i = it->value;
                auto& c = cylinders[i];
                xIntervals.emplace_back(c.llf.x, c.urb.x, i);
            }
            IntervalWalker<int, size_t> xWalker(0, std::move(xIntervals));

            for(int x = 0 ; x < size.x; x++){
                auto cylinderIt = xWalker.next();

                const vec3 offsets[8] {
                    vec3(-0.25, -0.25, -0.25),
                    vec3( 0.25, -0.25, -0.25),
                    vec3(-0.25,  0.25, -0.25),
                    vec3( 0.25,  0.25, -0.25),
                    vec3(-0.25, -0.25,  0.25),
                    vec3( 0.25, -0.25,  0.25),
                    vec3(-0.25,  0.25,  0.25),
                    vec3( 0.25,  0.25,  0.25),
                };
                const unsigned char outsideBase = 96;
                const unsigned char insideBase = 160;
                const unsigned char insideIncrement = (insideBase-outsideBase)/8;

                unsigned char val = outsideBase;
                for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
                    auto& i = it->value;
                    auto& c = cylinders[i];

                    vec3 imagePos(x, y, z);

                    unsigned char v = outsideBase;
                    for(auto& offset : offsets) {
                        vec3 subVoxelPos = imagePos + offset;
                        vec3 treePos = sampleToVT * subVoxelPos;
                        if(inSegment(treePos, c.p1, c.p2, c.radius)) {
                            v += insideIncrement;
                        }
                    }
                    val = std::max(v, val);
                }
                val = noise.apply(val);

                //char val = td->imageAt(i, j, k);
                svec3 slicePos(x,y,0);
                slice.set_voxel(slicePos, val);
            }
        }
        vol.writeSlice(slice, z);
    }

    return;
}

/**
 *  applies noise to the volumetric image
 *
 * noise.txt format:
 *
 * SHADOW: numShadows
 * GAUSSIAN: median sigma
 * UNIFORM: lb ub
 * SALTPEPER: valSalt probsalt valpepper probpepper
 *
 * noise will be added in the order specified by the file
 *
 * one image for each noise file will be generated
 */
void drawWithNoise(VascularTree& td, svec3 mapSize, svec3 size, float voxelWidth, const char* rootName, const char* noiseFile) {
    ifstream mapFile;
    mapFile.open(noiseFile);
    string line;

    double lb, ub, median, sigma, probSalt, probPepper;
    int numShadows;
    char valSalt, valPepper;
    bool first = true;

    if(mapFile.is_open()){
        while(!mapFile.eof()){
            getline(mapFile, line);
            char* tok = new char[line.size()];
            strcpy(tok, line.c_str());

            if(line.length() == 0)
                continue;

            char * field = strtok(tok, ":");
            if(!first) {
                throw "Sorry, combining multiple noise types supported in this fork";
            }

            if(strcmp(field, "SHADOW") == 0){

                //apply shadow noise
                numShadows = atoi(strtok(NULL, " "));

                throw "Sorry, shadow is not yet supported in this fork";

            } else if(strcmp(field, "GAUSSIAN") == 0){

                //apply gaussian noise
                median = atof(strtok(NULL, " "));
                sigma = atof(strtok(NULL, " "));

                drawImage(td, mapSize, size, voxelWidth, rootName, GaussianNoise { median, sigma });
            } else if(strcmp(field, "UNIFORM") == 0){

                //applying uniform noise
                lb = atof(strtok(NULL, " "));
                ub = atof(strtok(NULL, " "));

                throw "Sorry, uniform noise is not yet supported in this fork";

            } else {

                if (strcmp(field, "SALTPEPPER") == 0) {

                    //apply salt and pepper noise
                    valSalt = (char)atoi(strtok(NULL, " "));
                    probSalt = atof(strtok(NULL, " "));
                    valPepper = (char)atoi(strtok(NULL, " "));
                    probPepper = atof(strtok(NULL, " "));
                    valPepper = (char)valPepper;

                    throw "Sorry, salt and pepper noise is not yet supported in this fork";
                }
            }
            first = false;
        }

    } else {

        throw "Could not read the noise file";

    }

    mapFile.close();
}

/**
 * prints a node into XML/GXL format from the NodeTable
 */
void subPrint_node(NodeTable *nt, int segment, ofstream &os){
    os<<"  <node id=\"n"<<segment<<"\">"<<endl;
    os<<"    <attr name=\" nodeType\">"<<endl;
    if(nt->getType(segment) == NodeTable::ROOT){
        os<<"      <string> root node </string>"<<endl;
    } else if(nt->getType(segment) == NodeTable::TERM){
        os<<"      <string> terminal node </string>"<<endl;
    } else if(nt->getType(segment) == NodeTable::BIF){
        os<<"      <string> bifurication </string>"<<endl;
    } else {
        os<<"      <string> unknown type </string>"<<endl;
    }
    os<<"    </attr>"<<endl;

    os<<"    <attr name=\" position\">"<<endl;
    os<<"      <tup>"<<endl;
    double *pos = nt->getPos(segment);
    os<<"        <float>"<<pos[0]<<"</float>"<<endl;
    os<<"        <float>"<<pos[1]<<"</float>"<<endl;
    os<<"        <float>"<<pos[2]<<"</float>"<<endl;
    os<<"      </tup>"<<endl;
    os<<"    </attr>"<<endl;
    os<<"  </node>"<<endl;

    if(nt->getType(segment) != NodeTable::TERM){
        subPrint_node(nt, nt->getLeftChild(segment), os);
        if(nt->getType(segment) != NodeTable::ROOT)
            subPrint_node(nt, nt->getRightChild(segment), os);
    }
}

/**
 * prints an edge into XML/GXL format from a node table
 */
void subPrint_edge(NodeTable *nt, int segment, ofstream &os){

    if(nt->getType(segment) != NodeTable::ROOT){
        os<<"  <edge id=\"e"<<segment<<"\" to=\"n"<<segment<<"\" from=\"n"<<nt->getParent(segment)<<"\">"<<endl;
        os<<"    <attr name=\" flow\">"<<endl;
        os<<"      <float>"<<nt->getFlow(segment)<<"</float>"<<endl;
        os<<"    </attr>"<<endl;

        os<<"    <attr name=\" radius\">"<<endl;
        os<<"      <float>"<<nt->getRadius(segment)<<"</float>"<<endl;
        os<<"    </attr>"<<endl;

        os<<"  </edge>"<<endl;
    }

    if(nt->getType(segment) != NodeTable::TERM){
        subPrint_edge(nt, nt->getLeftChild(segment), os);
        if(nt->getType(segment) != NodeTable::ROOT)
            subPrint_edge(nt, nt->getRightChild(segment), os);
    }
}

/**
 * Writes the tree structure to a GXL file and stores information
 * about the nodes, edges, their heirarchy, radii and bifurcation locations
 */
void printTreeStructure(VascularTree * vt, const char * filePath){

    ofstream output;

    //writing the tree structure as GXL to the filePath specified
    output.open(filePath);
    output<<"<gxl><graph id=\""<<filePath<<"\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">"<<endl;

    //this seems really really stupid to do, why would we
    //loop through the entire structure to find the root, and then
    //recursively the nodes and edges?
    //TODO: fix this - should have the root always stored and easily
    //accessible and then recursively output the node/edges
    for(int i = 0; i < vt->nt.length; i++){
        if(vt->nt.getType(i) == NodeTable::ROOT){
            subPrint_node(&vt->nt, i, output);
            subPrint_edge(&vt->nt, i, output);
            output<<"</graph></gxl>"<<endl;
            output.close();
            return;
        }
    }

    output.close();

    throw "Unable to find root node.  The GXL file has not been generated.";

}


/**
 * VascuSynth: takes a series of parameter files, image names and (optionally) noise files
 * and generates a vascular structure based on the parameters.  The 3d volume is saved
 * as a series of 2D png slices.  Information about the vascular structure is saved as
 * a GXL file that can be visualized using software such as GraphViz.
 *
 * Arguments: paramFile.txt imageNames.txt voxelWidth noiseFiles.txt
 * noiseFiles.txt is an optional parameter.
 *
 * For each imageName/parameterFile/noiseFile, a folder is generated with the name specified
 * which will contain the 2D slices and the GXL file.
 */
int main(int argc, char** argv){

    int numFixedArgs = 5;
    if (argc < numFixedArgs) {
        //not enough parameters specified
        cout << "An error has occured: incorrect number of arguments" << endl;
        cout << "Usage: VascuSynth [paramFile] [imageNameFile] [voxelWidth] [outputSize]" << endl;
        return 0;
    }

    try {

        //read the param files and image name files
        vector<string> *paramFiles = readFileLines(argv[1]);
        vector<string> *imageNameFiles = readFileLines(argv[2]);

        //voxel widths
        string voxelWidth = argv[3];
        string outputSize = argv[4];
        string *noiseFiles = new string[argc-numFixedArgs];

        int paramFilesSize = (int) paramFiles->size();

        //go through each param file and build tree an spit it out
        for(int m = 0; m < paramFilesSize; m++){

            string paramFile = paramFiles->at(m);
            string rootDirectory = imageNameFiles->at(m);

            for(int i = numFixedArgs; i < argc; i++) {
                noiseFiles[i-numFixedArgs] = argv[i];
            }

            int numNoise = argc-numFixedArgs;

            cout << "Reading parameters and building the tree..." << endl;

            //build the tree
            VascularTree * vt = buildTree(paramFile.c_str());

            cout << "The vascular tree has been built sucessfully..." << endl;

            //filter out the /r that appears at times and messes up the directory name
            if (rootDirectory[ rootDirectory.length() - 1 ] == '\r') {
                rootDirectory = rootDirectory.substr(0, rootDirectory.length()-1);
            }

            //create the directory for the output
            rootDirectory = "./"+rootDirectory;
            mmkdir(rootDirectory.c_str());

            //create the GXL file
            string treeStructure = rootDirectory+"/tree_structure.xml";
            printTreeStructure(vt, treeStructure.c_str());

            cout << "The directory for the image has been created..." << endl;
            cout << "Information about the vascular structure has been saved in the gxl file tree_structure.xml..." << endl;

            double voxelWidthD = atof(voxelWidth.c_str());

            //create the subdirectory for the images
            string imageName = rootDirectory+"/original_image";
            mmkdir(imageName.c_str());

            //output the images
            imageName = imageName + "/image";
            svec3 outputSizeVec(atoi(outputSize.c_str()));

            drawImage(*vt, svec3(ivec3(vt->oxMap->dim)), outputSizeVec, voxelWidthD, imageName.c_str(), NoNoise{});

            cout << "The volumetric image has been saved..." << endl;

            if (numNoise > 0) {
                cout << "The images are being degraded by noise..." << endl;
            }

            //apply noise to the images - creating niose_images
            for(int i = 0; i < numNoise; i++){
                string noiseImage = rootDirectory+"/noise_image_"+itoa(i, 10);
                mmkdir(noiseImage.c_str());

                noiseImage = noiseImage+"/image";
                drawWithNoise(*vt, svec3(ivec3(vt->oxMap->dim)), outputSizeVec, voxelWidthD, noiseImage.c_str(), noiseFiles[i].c_str());
            }

            if (numNoise > 0) {
                cout << "Images have been succesfully degraded by noise and saved..." << endl;
            }

            //clean up
            delete vt;
        }

    } catch (string str) {
        cout << "ERROR: " << str << endl;
        cout << "Exiting VascuSynth" << endl;
    }

    return 0;
}

