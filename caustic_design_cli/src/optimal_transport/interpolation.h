#ifndef INTERPOLATION_H
#define INTERPOLATION_H
// STL
#include <map>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>

// local
#include "types.h"
#include "scene.h"

class Scene;
class VoronoiCreator;
class MainWindow;

class Interpolation{

public:
    Scene** computeScenes;
    Scene* source_scene;
    Scene* m_scene;
    Scene* compute_scene;
    MainWindow* win;
    std::vector<Point> Xo;
    std::map<Point, Vertex_handle> map;

    Interpolation(Scene* sc, Scene* tsc, Scene* csc, int sitesAmount, MainWindow* win);
    ~Interpolation(){}


    void runInterpolation(std::string imageFile, std::string datFile);

};

void multiThreadInterpolation(Interpolation* inter, uint id, uint from, uint to);
void findNeighborr(Point oP, Interpolation* inter, uint id, int index);

#endif // INTERPOLATION_H

