#ifndef SURFACEMODEL_H
#define SURFACEMODEL_H

#include "global.h"
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
/*Assimp Open Asset Import Librairy*/
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla
/*local*/
#include "SurfaceMesh.h"



using namespace std;


GLint TextureFromFile(const char* path, string directory);

class Model {
    public:
        /*  Functions   */
        Model(GLchar* path);
        Model();
        ~Model(){}
        void loadModel(string path);
        void loadReceiverLightPoints(QString path);
//        void Draw(Shader shader);
        int findSurfaceMesh();
        void printAllVertices();
        void exportModel(std::string filename);
        vector<glm::vec3> getLightRayPositions();
        glm::vec3 getFocalPlaneScale() { return this->targetPlaneScale; }
        glm::vec3 getFocalPlanePos() { return this->targetPlanePosition; }
        glm::quat getFocalPlaneQuat() { return this->targetPlaneRotationQuaternion; }

        void setFocalPlanePosX(float x);
        void setFocalPlanePosY(float y);
        void setFocalPlanePosZ(float z);

        void setFocalPlaneRotY(float y);
        void setFocalPlaneRotZ(float z);

        void setFocalPlaneScaleY(float y);
        void setFocalPlaneScaleZ(float z);

        void updateTargetPlaneRotationMatrix();

        void rescaleMeshes(float newScale);
        void modifyMesh();
        void shootRay(vector<glm::highp_dvec3> & direction, vector<glm::highp_dvec3> & redirect, vector<glm::highp_dvec3> & endpoint);
        void setNormals(bool edge);
        void computeLightDirectionsScreenSurface();
        void fresnelMapping(); //compute the position on the surface to for the desired normals
        /*  Model Data  */
        vector<Mesh> meshes; //in our case just one Mesh -> this is for more complex models
        vector<glm::vec3> desiredNormals;
        vector<glm::vec3> screenDirections;
        vector<glm::vec3> currentNormals;
        vector<glm::vec3> incidentNormals;
        vector<glm::vec3> receiverLightPositions;

        float surfaceSize;

    protected:
        aiScene* scene;
        float focalLength;

        glm::vec3 targetPlanePosition;
        glm::vec3 targetPlaneRotation;
        glm::vec3 targetPlaneScale;

        glm::quat targetPlaneRotationQuaternion;

        Mesh SurfaceMesh;
        //Mesh mesh;
        string directory;
        /*  Functions   */
        virtual void loadToSurface(int index){}
        virtual void processNode(aiNode* node, const aiScene* scene);
        virtual Mesh processMesh(aiMesh* mesh, const aiScene* scene);
        virtual vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

        void populateMesh(){}
        void compute(){}
        void clean(){
        meshes.clear();
        }
};

#endif // SURFACEMODEL_H
