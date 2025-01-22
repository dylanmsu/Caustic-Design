#include "targetoptimization.h"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;



/********************************************************/
/************** TargetOptimization Class ****************/
/********************************************************/

TargetOptimization::TargetOptimization(){}

TargetOptimization::~TargetOptimization(){}

/**
 * @brief TargetOptimization::gatherVertexInformation Gathers information about neighbors of current vertex.
 * @param vertex The current vertex
 * @param vertexIndex The index of the current vertex
 * @param neighborList [out] A list of indices to the neighbors of the current vertex. Each neighbor is only included once
 * @param neighborMap [out] A list of indices to indices to neighbors. The value is an index to the neighborList. Is used to handle several references to one neighbor (e.g. for neighboring faces)
 * @param gridNeighbors [out] A list indices to of horizontal and vertical neighbors
 */
void TargetOptimization::gatherVertexInformation(Vertex *vertex, uint vertexIndex, vector<int> &neighborList, vector<int> &neighborMap, vector<int> & gridNeighbors)
{
    Mesh * mesh = &model->meshes[0];
    vector<uint> adjacentFaces = mesh->edgeAdjacentFaces[vertexIndex];

    for(uint faceIndex = 0; faceIndex < adjacentFaces.size(); faceIndex++)
    {
        glm::uvec3 face = model->meshes[0].edgeIndices[adjacentFaces[faceIndex]];
        int thisVertexIndex = 0;
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            Vertex * v = mesh->faceVerticesEdge[face[vIndex]];

            if(v == vertex){
                thisVertexIndex = vIndex;
                break;
            }
        }

        int other[] = {
            (thisVertexIndex+1)%3,
            (thisVertexIndex+2)%3
        };


        for(uint j=0; j<2; j++)
        {
            std::vector<int>::iterator it;
            it = std::find(neighborList.begin(), neighborList.end(), face[other[j]]);
            if(it != neighborList.end() ) {
                // push back index of the neighbor
                neighborMap.push_back(it - neighborList.begin());
            } else {
                // push back neighbor and the index of it
                neighborList.push_back(face[other[j]]);
                neighborMap.push_back(neighborList.size()-1);
            }
        }
    }


    // add grid-neighbors (4 neighbors max, less possible if vertex is at edge)
    int row = mesh->vertexRowMap[vertexIndex];
    int col = mesh->vertexColMap[vertexIndex];

    if(row > 0)
        gridNeighbors.push_back(mesh->frontFaceMatrix[row-1][col]);

    if(col > 0)
        gridNeighbors.push_back(mesh->frontFaceMatrix[row][col-1]);
    if((col+1) < mesh->frontFaceMatrix[row].size())
        gridNeighbors.push_back(mesh->frontFaceMatrix[row][col+1]);

    if((row+1) < mesh->frontFaceMatrix.size())
        gridNeighbors.push_back(mesh->frontFaceMatrix[row+1][col]);


}

/**
 * @brief TargetOptimization::addResidualBlocks Adds residuals blocks for the given vertex to the given problem.
 * @param problem The problem to be solved
 * @param vertexIndex The vertex index (in the edgeVertices-vector of the mesh)
 * @param neighbors The neighboring-vertices of the current vertex. These include only the ones that share a face with the current vertex
 * @param neighborMap The index of the neighbors. Two neighbors make one face (with the vertex itself). So vertex + neighbors[neighborMap[0]] + neighbors[neighborMap[1]] forms one face
 * @param gridNeighbors The horizontal and vertical neighbors of the current vertex
 * @param vertices The reference to the vertices. Each vertex uses three doubles (x,y,z).
 */
void TargetOptimization::addResidualBlocks(Problem *problem, uint vertexIndex, vector<int> &neighbors, vector<int> &neighborMap, vector<int> & gridNeighbors, double *vertices)
{
    //std::cout << "vertices[vertexIndex*3] = {" 
    //    << vertices[vertexIndex*3 + 0] << ", " 
    //    << vertices[vertexIndex*3 + 1] << ", " 
    //    << vertices[vertexIndex*3 + 2] << "}" << std::endl;

    // EBar only calculates the distance from the vertex to the receiver-plane. So passing x-coordiate of receiving plane is sufficient
    //glm::vec3 receiverPos = glm::vec3(model->meshes[0].getMaxX() + model->getFocalLength(), 0, 0);
    glm::vec3 receiverPos = model->getFocalPlanePos();
    CostFunction* cost_function_ebar =
            new AutoDiffCostFunction<CostFunctorEbar2, 1, 3>(new CostFunctorEbar2(&receiverPos));

    problem->AddResidualBlock(
                cost_function_ebar,
                NULL, // no loss function
                &vertices[vertexIndex*3]
                );

    float weightMult = 1.0;
    if(model->meshes[0].edgeToNoEdgeMapping[vertexIndex] == -1) // we have an edge. Set weight for edir extremely high
        weightMult = 10000;

    // EDir depends on the original position
    CostFunction* cost_function_edir =
            new AutoDiffCostFunction<CostFunctorEdir2, 3, 3>(new CostFunctorEdir2(&x_sources[vertexIndex], weightMult));

    problem->AddResidualBlock(
                cost_function_edir,
                NULL, // no loss function
                &vertices[vertexIndex*3]
                );

    if(model->meshes[0].edgeToNoEdgeMapping[vertexIndex] != -1){ //not an edge we optimize the normals

        // For eint we have several functors, each for a different amount of neighbors
        switch(neighbors.size()){
        case 2:
        {

            CostFunction* cost_function_eint2 =
                new AutoDiffCostFunction<CostFunctorEint2Neighbors, 3, 3, 3, 3>(new CostFunctorEint2Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock(cost_function_eint2, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3]);


            break;
            }
        case 3:
        {

            CostFunction* cost_function_eint3 =
                new AutoDiffCostFunction<CostFunctorEint3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEint3Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock(cost_function_eint3, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3]);


            break;
            }

        case 4:
        {

            CostFunction* cost_function_eint4 =
                new AutoDiffCostFunction<CostFunctorEint4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEint4Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock(cost_function_eint4, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3],
                                       &vertices[neighbors[3]*3]);


            break;
            }

        case 5:
        {
            CostFunction* cost_function_eint5 =
                new AutoDiffCostFunction<CostFunctorEint5Neighbors, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint5Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock(cost_function_eint5, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3]);



            break;
        }

        case 6:
        {
            CostFunction* cost_function_eint6 =
                new AutoDiffCostFunction<CostFunctorEint6Neighbors, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint6Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock( cost_function_eint6, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3]);



            break;
        }

        case 7:
        {
            CostFunction* cost_function_eint7 =
                new AutoDiffCostFunction<CostFunctorEint7Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint7Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock( cost_function_eint7, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3],
                       &vertices[neighbors[6]*3]);


            break;
        }

        case 8:
        {
            CostFunction* cost_function_eint8 =
                new AutoDiffCostFunction<CostFunctorEint8Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint8Neighbors(&model->desiredNormals[model->meshes[0].edgeToNoEdgeMapping[vertexIndex]], neighborMap));

            problem->AddResidualBlock( cost_function_eint8, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3],
                           &vertices[neighbors[5]*3],
                           &vertices[neighbors[6]*3],
                           &vertices[neighbors[7]*3]);

            break;

        }

        } // switch end
    }

    // regularization term with 4-neighborhood
    CostFunction* ereg;
    switch(gridNeighbors.size())
    {
    case 2:

        ereg = new AutoDiffCostFunction<CostFunctorEreg2Neighbors, 3, 3, 3, 3>(new CostFunctorEreg2Neighbors(model, gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3]
                    );
        break;
    case 3:

        ereg = new AutoDiffCostFunction<CostFunctorEreg3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEreg3Neighbors(model, gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3],
                    &vertices[gridNeighbors[2]*3]
                    );
        break;
    case 4:

        ereg = new AutoDiffCostFunction<CostFunctorEreg4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEreg4Neighbors(model, gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3],
                    &vertices[gridNeighbors[2]*3],
                    &vertices[gridNeighbors[3]*3]
                    );
        break;
    }

}

bool Contains(const std::vector<int> &list, int x)
{
	return std::find(list.begin(), list.end(), x) != list.end();
}

void TargetOptimization::runOptimization(Model* m, Renderer* renderer){

    model = m;
    Mesh * mesh = &model->meshes[0];

    // make a copy of the original positions of the vertices
    for (int i = 0; i < mesh->faceVerticesEdge.size(); i++) {
        glm::vec3 v = mesh->faceVerticesEdge[i]->Position;
        x_sources.push_back(v);
    }

    vector<vector<int> > neighborsPerVertex;
    neighborsPerVertex.resize(mesh->faceVerticesEdge.size());

    vector<vector<int> > neighborMapPerVertex;
    neighborMapPerVertex.resize(mesh->faceVerticesEdge.size());

    vector<vector<int> > eightNeighborsPerVertex;
    eightNeighborsPerVertex.resize(mesh->faceVerticesEdge.size());

    // gather information for each vertex to optimize
    for(uint i = 0; i < mesh->faceVerticesEdge.size(); i++)
    {
        Vertex * v = mesh->faceVerticesEdge[i];
        vector<int> neighbors;
        vector<int> neighborMap;
        vector<int> eightNeighbors;

        gatherVertexInformation(v, i, neighbors, neighborMap, eightNeighbors);

        neighborsPerVertex[i] = neighbors;
        neighborMapPerVertex[i] = neighborMap;
        eightNeighborsPerVertex[i] = eightNeighbors;
    }

    // prepare model and mesh
    mesh->calculateVertexNormals();
    model->computeLightDirectionsScreenSurface();
    model->fresnelMapping();

    // optimize until converged or maximum step amount reached
    for(uint loop = 0; loop < 10; loop++)
    {
        // put all positions in one big list that we access later
        double* vertices = new double[3*mesh->faceVerticesEdge.size()];
        for(uint i=0; i<mesh->faceVerticesEdge.size(); i++)
        {
            glm::vec3 * pos = &mesh->faceVerticesEdge[i]->Position;

            vertices[3*i + 0] = pos->x;
            vertices[3*i + 1] = pos->y;
            vertices[3*i + 2] = pos->z;
        }

        bool duplicate = false;

        // check for duplicate items 
        /*for (int i=0; i<mesh->faceVerticesEdge.size(); i++) {
            for (int j=0; j<mesh->faceVerticesEdge.size(); j++) {
                if (i != j) {
                    if (
                        floor(vertices[3*i + 0]*100000) == floor(vertices[3*j + 0]*100000) &&
                        floor(vertices[3*i + 1]*100000) == floor(vertices[3*j + 1]*100000) &&
                        floor(vertices[3*i + 2]*100000) == floor(vertices[3*j + 2]*100000)) {
                        std::cout << "has same" << std::endl;
                    }
                }
            }
        }*/

        Problem prob;

        // iterate over all vertices and add the corresponding residual blocks
        for(uint i=0; i<neighborsPerVertex.size(); i++)
        {
            //std::cout << "block " << i << std::endl;
            //std::cout << "neighborsPerVertex " << neighborsPerVertex[i].size() <<  std::endl;
            //std::cout << "neighborMapPerVertex " << neighborMapPerVertex[i].size() <<  std::endl;
            //std::cout << "eightNeighborsPerVertex " << eightNeighborsPerVertex[i].size() <<  std::endl;
            addResidualBlocks(&prob, i, neighborsPerVertex[i], neighborMapPerVertex[i], eightNeighborsPerVertex[i], vertices);
        }

        Solver::Options options;
        options.minimizer_progress_to_stdout = true;
        options.linear_solver_type = ceres::ITERATIVE_SCHUR;
        options.max_num_iterations = 200;
        options.dense_linear_algebra_library_type = ceres::LAPACK;
        options.num_threads = 16;

        string error;
        if(!options.IsValid(&error))
        {
            std::cout << "Options not valid: " << error << std::endl;
        }

        Solver::Summary summary;
        Solve(options, &prob, &summary);

        std::cout << summary.FullReport() << std::endl;

        glm::vec3 * pos;
        for(uint i=0; i<mesh->faceVerticesEdge.size(); i++)
        {
            pos = &mesh->faceVerticesEdge[i]->Position;
            pos->x = vertices[3*i + 0];
            pos->y = vertices[3*i + 1];
            pos->z = vertices[3*i + 2];
        }

        delete[] vertices;

        mesh->calculateVertexNormals();
        model->computeLightDirectionsScreenSurface();
        model->fresnelMapping();
        renderer->repaint();

        // check for convergence
        if((summary.num_successful_steps + summary.num_unsuccessful_steps) == 0){
            std::cout << std::endl << "Outer Loop Converged" << std::endl << std::endl;
            break;
        }

    }
    model->meshes[0].calculateVertexNormals();

    renderer->repaint();
}


