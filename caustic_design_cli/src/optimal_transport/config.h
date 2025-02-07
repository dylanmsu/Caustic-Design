#ifndef CONFIG_H
#define CONFIG_H
    /*Voronoi Generation*/
    #define EPSILON 0.1
    #define LLOYD_STEPS 200

    /*Optimal Transport*/
    //#define LEVEL_MAX 6
    //int LEVEL_MAX = 6;

    /*Interpolation*/
    #define MESH_AMOUNT 50
    #define SCALING_X 0.5
    #define SCALING_Y 0.5

    /* debug and demo */
    #define LIVE_DEMO

    /* manual optimization implementation */
    //#define DESCENT_GRADIENT

    #define DOMAIN_WIDTH 0.5
#endif // CONFIG_H
