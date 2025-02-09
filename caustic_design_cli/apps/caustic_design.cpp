// This file is part of otmap, an optimal transport solver.
//
// Copyright (C) 2017-2018 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2017 Georges Nader
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#include "normal_integration/normal_integration.h"
#include "normal_integration/mesh.h"

#include "optimal_transport/scene.h"
//#include "voronoi_creation.h"
#include "optimal_transport/optimal_transport.h"
#include "optimal_transport/interpolation.h"
#include "optimal_transport/config.h"

template<typename T,typename S>
T lerp(S u, const T& a0, const T& a1)
{
  return (1.-u)*a0 + u*a1;
}

// Function to calculate the gradient of f(y, z)
void gradient(  std::vector<double> source,
                std::vector<double> interf,
                std::vector<double> target,
                double n1, double n2,
                double & grad_x, double & grad_y) {
    double d1 = std::sqrt((interf[0] - source[0]) * (interf[0] - source[0]) + (interf[1] - source[1]) * (interf[1] - source[1]) + (interf[2] - source[2]) * (interf[2] - source[2]));
    double d2 = std::sqrt((target[0] - interf[0]) * (target[0] - interf[0]) + (target[1] - interf[1]) * (target[1] - interf[1]) + (target[2] - interf[2]) * (target[2] - interf[2]));

    grad_x = n1 * (interf[0] - source[0]) / d1 - n2 * (target[0] - interf[0]) / d2;
    grad_y = n1 * (interf[1] - source[1]) / d1 - n2 * (target[1] - interf[1]) / d2;
}

void scaleAndTranslatePoints(std::vector<std::vector<double>>& points, double MAX_X, double MAX_Y, double margin) {
    double scaleFactorX = (MAX_X - 2 * margin) / MAX_X;
    double scaleFactorY = (MAX_Y - 2 * margin) / MAX_Y;

    for (auto& point : points) {
        double& x = point[0];
        double& y = point[1];

        x = x * scaleFactorX;
        y = y * scaleFactorY;

        x += margin;
        y += margin;
    }
}

void export_grid_to_svg(std::vector<std::vector<double>> &points, double width, double height, int res_x, int res_y, std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * (height / width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    for (int j = 0; j < res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < res_x; i++) {
            int idx = i + j * res_x;

            const auto& point = points[idx];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));
            if (i < res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    for (int j = 0; j < res_x; j++) {
        std::string path_str = "M";
        for (int i = 0; i < res_y; i++) {
            int idx = j + i * res_x;

            const auto& point = points[idx];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));

            if (i < res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

void export_triangles_to_svg(std::vector<std::vector<double>> &points, std::vector<std::vector<unsigned int>> &triangles, double width, double height, int res_x, int res_y, std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * (height / width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
    
    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Draw polygons
    for (const auto& polygon : triangles) {
        std::vector<std::vector<double>> poly_points;
        for (const auto& point_idx : polygon) {
            poly_points.push_back(points[point_idx]);
        }

        std::string path_str = "M";
        for (std::size_t j = 0; j < poly_points.size(); ++j) {
            const auto& point = poly_points[j];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));

            if (j < poly_points.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

void scalePoints(std::vector<std::vector<double>>& trg_pts, const std::vector<double>& scale, const std::vector<double>& origin) {
    for (std::vector<double>& point : trg_pts) {
        for (size_t j = 0; j < point.size(); ++j) {
            // Scale each dimension relative to the origin
            point[j] = origin[j] + (point[j] - origin[j]) * scale[j];
        }
    }
}

void translatePoints(std::vector<std::vector<double>>& trg_pts, std::vector<double> position_xyz) {
  for (int i = 0; i < trg_pts.size(); i++)
  {
    trg_pts[i][0] += position_xyz[0];
    trg_pts[i][1] += position_xyz[1];
    trg_pts[i][2] += position_xyz[2];
  }
}

// Define the rotation function
void rotatePoints(std::vector<std::vector<double>>& trg_pts, std::vector<double> angle_xyz) {
    double PI = 3.14159265358979323846;

    // Convert angles from degrees to radians
    angle_xyz[0] = angle_xyz[0] * PI / 180.0;
    angle_xyz[1] = angle_xyz[1] * PI / 180.0;
    angle_xyz[2] = angle_xyz[2] * PI / 180.0;

    // Precompute sine and cosine values for each rotation angle
    double cos_x = std::cos(angle_xyz[0]);
    double sin_x = std::sin(angle_xyz[0]);
    double cos_y = std::cos(angle_xyz[1]);
    double sin_y = std::sin(angle_xyz[1]);
    double cos_z = std::cos(angle_xyz[2]);
    double sin_z = std::sin(angle_xyz[2]);

    // Define the rotation matrices for each axis
    std::vector<std::vector<double>> rot_x = {
        {1, 0, 0},
        {0, cos_x, -sin_x},
        {0, sin_x, cos_x}
    };

    std::vector<std::vector<double>> rot_y = {
        {cos_y, 0, sin_y},
        {0, 1, 0},
        {-sin_y, 0, cos_y}
    };

    std::vector<std::vector<double>> rot_z = {
        {cos_z, -sin_z, 0},
        {sin_z, cos_z, 0},
        {0, 0, 1}
    };

    // Apply rotation to each point in the point cloud
    for (std::vector<double>& point : trg_pts) {
        // Extract x, y, z for clarity
        double x = point[0];
        double y = point[1];
        double z = point[2];

        // Rotate around x-axis
        double new_y = rot_x[1][1] * y + rot_x[1][2] * z;
        double new_z = rot_x[2][1] * y + rot_x[2][2] * z;
        y = new_y;
        z = new_z;

        // Rotate around y-axis
        double new_x = rot_y[0][0] * x + rot_y[0][2] * z;
        new_z = rot_y[2][0] * x + rot_y[2][2] * z;
        x = new_x;
        z = new_z;

        // Rotate around z-axis
        new_x = rot_z[0][0] * x + rot_z[0][1] * y;
        new_y = rot_z[1][0] * x + rot_z[1][1] * y;
        x = new_x;
        y = new_y;

        // Update the point with the rotated coordinates
        point[0] = x;
        point[1] = y;
        point[2] = z;
    }
}

std::vector<double> normalize_vec(std::vector<double> p1) {
    std::vector<double> vec(3);
    double squared_len = 0;
    for (int i=0; i<p1.size(); i++) {
        squared_len += p1[i] * p1[i];
    }

    double len = std::sqrt(squared_len);

    for (int i=0; i<p1.size(); i++) {
        vec[i] = p1[i] / len;
    }

    return vec;
}

//compute the desired normals
std::vector<std::vector<double>> fresnelMapping(std::vector<std::vector<double>> &vertices, std::vector<std::vector<double>> &target_pts, double refractive_index) {
    std::vector<std::vector<double>> desiredNormals;

    //double boundary_z = -0.1;

    //vector<std::vector<double>> boundary_points;

    bool use_point_src = false;
    bool use_reflective_caustics = false;

    //std::vector<double> pointLightPosition(3);
    //pointLightPosition[0] = 0.5;
    //pointLightPosition[1] = 0.5;
    //pointLightPosition[2] = 0.5;

    // place initial points on the refractive surface where the light rays enter the material
    /*if (use_point_src && !use_reflective_caustics) {
        for(int i = 0; i < vertices.size(); i++) {
            std::vector<double> boundary_point(3);

            // ray to plane intersection to get the initial points
            double t = ((boundary_z - pointLightPosition[2]) / (vertices[i][2] - pointLightPosition[2]));
            boundary_point[0] = pointLightPosition[0] + t*(vertices[i][0] - pointLightPosition[0]);
            boundary_point[1] = pointLightPosition[1] + t*(vertices[i][1] - pointLightPosition[1]);
            boundary_point[2] = boundary_z;
            boundary_points.push_back(boundary_point);
        }
    }*/

    // run gradient descent on the boundary points to find their optimal positions such that they satisfy Fermat's principle
    /*if (!use_reflective_caustics && use_point_src) {
        for (int i=0; i<boundary_points.size(); i++) {
            for (int iteration=0; iteration<100000; iteration++) {
                double grad_x;
                double grad_y;
                gradient(pointLightPosition, boundary_points[i], vertices[i], 1.0, refractive_index, grad_x, grad_y);

                boundary_points[i][0] -= 0.1 * grad_x;
                boundary_points[i][1] -= 0.1 * grad_y;

                // if magintude of both is low enough
                if (grad_x*grad_x + grad_y*grad_y < 0.000001) {
                    break;
                }
            }
        }
    }*/

    for(int i = 0; i < vertices.size(); i++) {
        std::vector<double> incidentLight(3);
        std::vector<double> transmitted = {
            target_pts[i][0] - vertices[i][0],
            target_pts[i][1] - vertices[i][1],
            target_pts[i][2] - vertices[i][2]
        };

        /*if (use_point_src) {
            incidentLight[0] = vertices[i][0] - pointLightPosition[0];
            incidentLight[1] = vertices[i][1] - pointLightPosition[1];
            incidentLight[2] = vertices[i][2] - pointLightPosition[2];
        } else {*/
            incidentLight[0] = 0;
            incidentLight[1] = 0;
            incidentLight[2] = -1;
        //}

        transmitted = normalize_vec(transmitted);
        incidentLight = normalize_vec(incidentLight);

        std::vector<double> normal(3);
        /*if (use_reflective_caustics) {
            normal[0] = ((transmitted[0]) - incidentLight[0]) * 1.0f;
            normal[1] = ((transmitted[1]) - incidentLight[1]) * 1.0f;
            normal[2] = ((transmitted[2]) - incidentLight[2]) * 1.0f;
        } else {*/
            normal[0] = ((transmitted[0]) - (incidentLight[0]) * refractive_index) * -1.0f;
            normal[1] = ((transmitted[1]) - (incidentLight[1]) * refractive_index) * -1.0f;
            normal[2] = ((transmitted[2]) - (incidentLight[2]) * refractive_index) * -1.0f;
        //}

        normal = normalize_vec(normal);

        desiredNormals.push_back(normal);
    }

    return desiredNormals;
}


std::unordered_map<std::string, std::string> parse_arguments(int argc, char const *argv[]) {
    // Define a map to store the parsed arguments
    std::unordered_map<std::string, std::string> args;

    // Iterate through command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string key, value;

        // Check if argument starts with '--'
        if (arg.substr(0, 2) == "--") {
            // Split argument by '=' to separate key and value
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                key = arg.substr(2, pos - 2);
                value = arg.substr(pos + 1);
            }
            else {
                key = arg.substr(2);
                value = ""; // No value provided
            }
        }
        // Check if argument starts with '-'
        else if (arg[0] == '-') {
            // The next argument is the value
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                key = arg.substr(1);
                value = argv[++i];
            }
            else {
                key = arg.substr(1);
                value = ""; // No value provided
            }
        }
        // Invalid argument format
        else {
            //std::cerr << "Invalid argument format: " << arg << std::endl;
            //return 1;
        }

        // Store key-value pair in the map
        args[key] = value;
    }

    return args;
}

// interpolate target mesh into a rectangular grid
std::vector<double> interpolate_point(Mesh &mesh, std::vector<std::vector<double>>& positions, std::vector<double> &point, bool &triangle_miss) {
    Hit hit;
    bool intersection = false;
    mesh.source_bvh->query(point, hit, intersection);
    if (intersection) {
        double interpolation_x = positions[mesh.triangles[hit.face_id][0]][0]*hit.barycentric_coords[0] + positions[mesh.triangles[hit.face_id][1]][0]*hit.barycentric_coords[1] + positions[mesh.triangles[hit.face_id][2]][0]*hit.barycentric_coords[2];
        double interpolation_y = positions[mesh.triangles[hit.face_id][0]][1]*hit.barycentric_coords[0] + positions[mesh.triangles[hit.face_id][1]][1]*hit.barycentric_coords[1] + positions[mesh.triangles[hit.face_id][2]][1]*hit.barycentric_coords[2];
        triangle_miss = false;
        return {interpolation_x, interpolation_y, 0.0f};
    } else {
        printf("interpolation miss!\r\n");
        printf("x: %f, y: %f\r\n", point[0], point[1]);
        //exit(0);
        triangle_miss = true;
        return point;
    }

    //printf("interpolation miss!\r\n");

    return {NAN, NAN};
}

std::vector<double> calculate_centroid_vector(std::vector<Point> vertices) {
    std::vector<double> centroid;
    centroid.push_back(0.0);
    centroid.push_back(0.0);

    double signed_area = 0;

    for (int i = 0; i < vertices.size(); i++) {
        double x0 = vertices[i].x();
        double y0 = vertices[i].y();
        double x1 = vertices[(i + 1) % vertices.size()].x();
        double y1 = vertices[(i + 1) % vertices.size()].y();

        // Shoelace formula
        double area = (x0 * y1) - (x1 * y0);
        signed_area += area;
        centroid[0] += (x0 + x1) * area;
        centroid[1] += (y0 + y1) * area;
    }

    signed_area *= 0.5;
    centroid[0] /= 6 * signed_area;
    centroid[1] /= 6 * signed_area;

    return centroid;
}

int main(int argc, char const *argv[])
{
    setlocale(LC_ALL,"C");

    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    std::vector<Eigen::Vector2d> vertex_positions;
    normal_integration normal_int;

    int mesh_res = 100;
    int n_sites = 10000;
    int n_levels = 4;
    float focal_l = 1.5 / 5;
    std::string source_image_filename = "";
    std::string target_image_filename = "";

    n_sites = std::stoi(args["n_sites"]);
    n_levels = std::stoi(args["n_levels"]);
    mesh_res = std::stoi(args["mesh_res"]);
    focal_l = std::stod(args["focal_l"]) / 5;
    source_image_filename = args["source_png"].c_str();
    target_image_filename = args["target_png"].c_str();

    Scene* source_scene;
    source_scene = new Scene;
    source_scene->load_image(source_image_filename);

    Scene* m_scene;
    m_scene = new Scene;
    m_scene->load_image(target_image_filename);

    //Mesh mesh(1.0, 1.0/2, opts.resolution, (int)(opts.resolution/2));
    Mesh mesh(1.0, 1.0, mesh_res, mesh_res);
    
    mesh.build_vertex_to_triangles();

    mesh.calculate_vertex_laplacians();

    normal_int.initialize_data(mesh);
    
    OptimalTransport ot = OptimalTransport(m_scene, source_scene, n_levels, n_sites);
    ot.runOptimalTransport(false);

    std::vector<std::vector<double>> points;
    std::vector<std::vector<double>> pd_centroids;
    std::vector<std::vector<unsigned int>> triangles;
    std::vector<std::vector<double>> target_points;

    for (auto fit = source_scene->m_rt.finite_faces_begin(); fit != source_scene->m_rt.finite_faces_end(); ++fit) {
        std::vector<unsigned int> triangle = {
            static_cast<unsigned int>(fit->vertex(0)->get_index()),
            static_cast<unsigned int>(fit->vertex(1)->get_index()),
            static_cast<unsigned int>(fit->vertex(2)->get_index())
        };
        triangles.push_back(triangle);
    }

    for (int i = 0; i < source_scene->m_vertices.size(); i++)
    {
        points.push_back({
            source_scene->m_vertices[i]->get_position().x() + 0.5,
            source_scene->m_vertices[i]->get_position().y() + 0.5
        });
    }

    for (unsigned i = 0; i < source_scene->m_vertices.size(); ++i) 
    {
        Vertex_handle vi = source_scene->m_vertices[i];
        std::vector<Point> polygon;

        if (vi->is_hidden()) continue;
        source_scene->m_rt.build_polygon(vi, polygon);

        pd_centroids.push_back(calculate_centroid_vector(polygon));
        pd_centroids[i][0] += 0.5;
        pd_centroids[i][1] += 0.5;
    }

    Mesh interpolation_mesh(points, triangles);

    interpolation_mesh.build_source_bvh(5, 30);
    for (int i = 0; i < mesh.source_points.size(); i++)
    {
        bool triangle_miss = false;
        target_points.push_back(interpolate_point(interpolation_mesh, pd_centroids, mesh.source_points[i], triangle_miss));
        target_points[i].push_back(0.0f);
    }
  
    //export_grid_to_svg(mesh.source_points, 1, 1, opts.resolution, opts.resolution, "../grid.svg", 0.5);

    //scaleAndTranslatePoints(mesh.source_points, 1.0, 1.0, 1.0 / 100);
    
    for (int i=0; i<mesh.source_points.size(); i++)
    {
        Eigen::Vector2d point = {mesh.source_points[i][0], mesh.source_points[i][1]};
        vertex_positions.push_back(point);
    }

    //export_grid_to_svg(trg_pts, 1, 0.5, opts.resolution, opts.resolution, "../grid.svg", 0.5);

    std::vector<std::vector<double>> desired_normals;

    //scalePoints(trg_pts, {8, 8, 0}, {0.5, 0.5, 0});
    rotatePoints(target_points, {0, 0, 0});
    translatePoints(target_points, {0, 0, focal_l});

    for (int i=0; i<10; i++)
    {
        std::vector<std::vector<double>> normals = fresnelMapping(mesh.source_points, target_points, 1.49);

        normal_int.perform_normal_integration(mesh, normals);
    }

    mesh.save_solid_obj_source(0.2, "output.obj");
}
