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

std::vector<double> cross(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result(3);
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return result;
}

double dot(std::vector<double> a, std::vector<double> b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::vector<double> mult(double a, std::vector<double> b) {
  return {a*b[0], a*b[1], a*b[2]};
}

std::vector<double> add(std::vector<double> a, std::vector<double> b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double magnitude(std::vector<double> a) {
  return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

double cot(const std::vector<double>& a, const std::vector<double>& b) {
    auto cross_product = cross(a, b);
    double cross_magnitude = magnitude(cross_product);

    if (cross_magnitude < 1e-12) {
        //throw std::invalid_argument("Vectors are parallel or one is a zero vector, cotangent undefined.");
        cross_magnitude = 1e-12;
    }

    return dot(a, b) / cross_magnitude;
}

bool is_boundary_vertex(Mesh &mesh, std::vector<std::pair<int, int>> &adjacent_edges, std::vector<int> &adjacent_triangles, int vertex_index, std::vector<std::pair<int, int>>& boundary_edges) {
    std::unordered_map<std::pair<int, int>, int, HashPair> edge_triangle_count;
    for (int triangle_index : adjacent_triangles) {
        const std::vector<unsigned int>& triangle = mesh.triangles[triangle_index];
        for (int j = 0; j < 3; ++j) {
            int v1 = triangle[j];
            int v2 = triangle[(j + 1) % 3];
            std::pair<int, int> edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            edge_triangle_count[edge]++;
        }
    }

    bool is_boundary = false;
    for (const auto& edge : adjacent_edges) {
        if (edge_triangle_count[edge] == 1) { // Boundary edge
            boundary_edges.push_back(edge);
            is_boundary = true;
        }
    }

    return is_boundary;
}

void project_onto_boundary(std::vector<double> &point) {
  point[0] -= 0.5;
  point[1] -= 0.5;

  double dist = sqrt(pow(point[0], 2) + pow(point[1], 2))*2;

  point[0] /= dist;
  point[1] /= dist;

  point[0] += 0.5;
  point[1] += 0.5;
}

std::vector<double> compute_laplacian(Mesh &mesh, std::vector<int> adjacent_triangles, std::vector<int> neighboring_vertices, int i) {
    std::vector<double> laplacian(neighboring_vertices.size(), 0.0f);

    for (int j_index = 0; j_index < neighboring_vertices.size(); ++j_index) {
        int j = neighboring_vertices[j_index];

        // Find triangles shared between `i` and `j`
        std::vector<int> shared_triangles;
        for (int triangle : adjacent_triangles) {
          // Check if `j` is one of the vertices in this triangle
          const auto& vertices = mesh.triangles[triangle];
          if (std::find(vertices.begin(), vertices.end(), j) != vertices.end()) {
              shared_triangles.push_back(triangle);
          }
        }

        // Handle cases based on the number of shared triangles
        if (shared_triangles.size() == 2) {
            // Interior case: Two triangles are connected
            std::vector<int> k_vertices;
            for (int triangle : shared_triangles) {
                for (int vertex : mesh.triangles[triangle]) {
                    if (vertex != i && vertex != j) {
                        k_vertices.push_back(vertex);
                        break; // Only one `k` per triangle
                    }
                }
            }

            // Ensure we found two `k` vertices
            if (k_vertices.size() != 2) {
                throw std::runtime_error("Error identifying k vertices in triangles.");
            }

            int k1 = k_vertices[0];
            int k2 = k_vertices[1];

            std::vector<double> edge1;
            std::vector<double> edge2;

            edge1 = sub(mesh.source_points[k1], mesh.source_points[j]);
            edge2 = sub(mesh.source_points[k1], mesh.source_points[i]);
            double cot_k1 = cot(edge1, edge2);

            edge1 = sub(mesh.source_points[k2], mesh.source_points[j]);
            edge2 = sub(mesh.source_points[k2], mesh.source_points[i]);
            double cot_k2 = cot(edge1, edge2);

            laplacian[j_index] += cot_k1;
            laplacian[j_index] += cot_k2;

            //std::cout << "k1=" << k1 << ", k2=" << k2 << std::endl;

        } else if (shared_triangles.size() == 1) {
            // Boundary case: Only one triangle is connected
            int triangle = shared_triangles[0];
            int k = -1;

            // Find the single `k` vertex
            for (int vertex : mesh.triangles[triangle]) {
                if (vertex != i && vertex != j) {
                    k = vertex;
                    break;
                }
            }

            if (k == -1) {
                throw std::runtime_error("Error identifying k vertex in boundary triangle.");
            }

            std::vector<double> edge1;
            std::vector<double> edge2;

            edge1 = sub(mesh.source_points[k], mesh.source_points[j]);
            edge2 = sub(mesh.source_points[k], mesh.source_points[i]);
            double cot_k = cot(edge1, edge2);

            laplacian[j_index] += cot_k;

            //std::cout << "k=" << k << std::endl;

        } else {
            throw std::runtime_error("No shared triangles between i and j; invalid mesh or disconnected vertex.");
        }
    }

    return laplacian;
}

//compute the desired normals
std::vector<std::vector<double>> fresnelMapping(
  std::vector<std::vector<double>> &vertices, 
  std::vector<std::vector<double>> &target_pts, 
  double refractive_index
) {
    std::vector<std::vector<double>> desiredNormals;

    //double boundary_z = -0.1;

    //vector<std::vector<double>> boundary_points;

    bool use_point_src = false;
    bool use_reflective_caustics = false;

    std::vector<double> pointLightPosition(3);
    pointLightPosition[0] = 0.5;
    pointLightPosition[1] = 0.5;
    pointLightPosition[2] = 0.5;

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

        if (use_point_src) {
            incidentLight[0] = vertices[i][0] - pointLightPosition[0];
            incidentLight[1] = vertices[i][1] - pointLightPosition[1];
            incidentLight[2] = vertices[i][2] - pointLightPosition[2];
        } else {
            incidentLight[0] = 0;
            incidentLight[1] = 0;
            incidentLight[2] = -1;
        }

        //transmitted = normalize_vec(transmitted);
        //incidentLight = normalize_vec(incidentLight);

        std::vector<double> normal(3);
        if (use_reflective_caustics) {
            normal[0] = ((transmitted[0]) - incidentLight[0]) * 1.0f;
            normal[1] = ((transmitted[1]) - incidentLight[1]) * 1.0f;
            normal[2] = ((transmitted[2]) - incidentLight[2]) * 1.0f;
        } else {
            normal[0] = ((transmitted[0]) - (incidentLight[0]) * refractive_index) * -1.0f;
            normal[1] = ((transmitted[1]) - (incidentLight[1]) * refractive_index) * -1.0f;
            normal[2] = ((transmitted[2]) - (incidentLight[2]) * refractive_index) * -1.0f;
        }

        //normal = normalize_vec(normal);

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
std::vector<double> interpolate_raster_source(Mesh *mesh, std::vector<std::vector<double>>& positions, std::vector<double> &point, bool &triangle_miss) {
    mesh->build_source_bvh(5, 30);
    
    Hit hit;
    bool intersection = false;
    mesh->source_bvh->query(point, hit, intersection);
    if (intersection) {
        double interpolation_x = positions[mesh->triangles[hit.face_id][0]][0]*hit.barycentric_coords[0] + positions[mesh->triangles[hit.face_id][1]][0]*hit.barycentric_coords[1] + positions[mesh->triangles[hit.face_id][2]][0]*hit.barycentric_coords[2];
        double interpolation_y = positions[mesh->triangles[hit.face_id][0]][1]*hit.barycentric_coords[0] + positions[mesh->triangles[hit.face_id][1]][1]*hit.barycentric_coords[1] + positions[mesh->triangles[hit.face_id][2]][1]*hit.barycentric_coords[2];
        triangle_miss = false;
        return {interpolation_x, interpolation_y};
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

  Scene* source_scene;
  source_scene = new Scene;
  source_scene->load_image(args["source_png"].c_str());

  Scene* m_scene;
  m_scene = new Scene;
  m_scene->load_image(args["target_png"].c_str());

  //Mesh mesh(1.0, 1.0/2, opts.resolution, (int)(opts.resolution/2));
  Mesh mesh(1.0, 1.0, 100, 100);
  
  mesh.build_vertex_to_triangles();

  normal_int.initialize_data(mesh);
  
  OptimalTransport ot = OptimalTransport(m_scene, source_scene, 4, 5000);
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

    printf("test\r\n");

    Mesh interpolation_mesh(points, triangles);

    printf("test2\r\n");

    for (int i = 0; i < mesh.source_points.size(); i++)
    {
        bool triangle_miss = false;
        target_points.push_back(interpolate_raster_source(&interpolation_mesh, pd_centroids, mesh.source_points[i], triangle_miss));
    }

    printf("testt\r\n");
  
  //export_grid_to_svg(mesh.source_points, 1, 1, opts.resolution, opts.resolution, "../grid.svg", 0.5);

  scaleAndTranslatePoints(mesh.source_points, 1.0, 1.0, 1.0 / 100);
  
  for (int i=0; i<mesh.source_points.size(); i++)
  {
    Eigen::Vector2d point = {mesh.source_points[i][0], mesh.source_points[i][1]};
    vertex_positions.push_back(point);
  }

  printf("testtt\r\n");

  //export_grid_to_svg(trg_pts, 1, 0.5, opts.resolution, opts.resolution, "../grid.svg", 0.5);

  std::vector<std::vector<double>> desired_normals;

  //scalePoints(trg_pts, {8, 8, 0}, {0.5, 0.5, 0});
  rotatePoints(target_points, {0, 0, 0});
  translatePoints(target_points, {0, 0, 1.5});

  double r = 1.55;

  printf("testttt\r\n");

  for (int i=0; i<10; i++)
  {
      double max_z = -10000;

      for (int j = 0; j < mesh.source_points.size(); j++) {
        if (max_z < mesh.source_points[j][2]) {
          max_z = mesh.source_points[j][2];
        }
      }

      for (int j = 0; j < mesh.source_points.size(); j++) {
          mesh.source_points[j][2] -= max_z;
      }
      

      std::vector<std::vector<double>> normals = fresnelMapping(mesh.source_points, target_points, r);

      normal_int.perform_normal_integration(mesh, normals);

      printf("testtttt\r\n");

      //std::vector<std::vector<double>> vertex_normals = normal_int.calculate_vertex_normals(mesh);

      //std::vector<double> incidentLight(3);
      //incidentLight[0] = 0;
      //incidentLight[1] = 0;
      //incidentLight[2] = -1;

      //std::vector<int32_t> plane_triangle = mesh.triangles[0];

      //std::vector<double> plane_normal = calc_plane_normal(trg_pts[plane_triangle[0]], trg_pts[plane_triangle[0]], trg_pts[plane_triangle[0]]);

      //std::vector<std::vector<double>> intersections(vertex_normals.size());

      /*std::vector<double> pointLightPosition(3);
      pointLightPosition[0] = 0.5;
      pointLightPosition[1] = 0.5;
      pointLightPosition[2] = 0.5;

      for (int i = 0; i < vertex_normals.size(); i++)
      {
        vertex_normals[i][0] *= -1.0;
        vertex_normals[i][1] *= -1.0;
        vertex_normals[i][2] *= -1.0;

        std::vector<double> incidentLight(3);
        incidentLight[0] = mesh.source_points[i][0] - pointLightPosition[0];
        incidentLight[1] = mesh.source_points[i][1] - pointLightPosition[1];
        incidentLight[2] = mesh.source_points[i][2] - pointLightPosition[2];

        std::vector<double> intersectionPoint(3);
        std::vector<double> refracted = refract(vertex_normals[i], incidentLight, r, 1.0);

        refracted[0] /= -(refracted[2] + mesh.source_points[i][2]);
        refracted[1] /= -(refracted[2] + mesh.source_points[i][2]);
        refracted[2] /= -(refracted[2] + mesh.source_points[i][2]);

        refracted[0] += mesh.source_points[i][0];
        refracted[1] += mesh.source_points[i][1];
        refracted[2] += mesh.source_points[i][2];

        intersections[i] = refracted;
      }*/

      //export_triangles_to_svg(intersections, mesh.triangles, 1, 1, opts.resolution, opts.resolution, "../triangles.svg", 0.5);
      //export_grid_to_svg(intersections, 1, 1, opts.resolution, opts.resolution, "../intersections.svg", 0.5);

      //mesh.save_solid_obj_source(0.4, "../output.obj");
  }

  mesh.save_solid_obj_source(0.2, "../output.obj");
}
