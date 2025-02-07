#include "scene.h"
#include "util.h"
#include "random.h"

bool Scene::is_valid() const
{
    if (!m_domain.is_valid()) return false;
    if (m_vertices.empty()) return false;
    return true;
}

void Scene::load_image(const std::string& filename, const int width)
{
    bool ok = m_domain.load(filename, width);
    if (!ok) {
        std::cerr << "failed to load image: " << filename << std::endl;
    }
    
    m_rt.set_boundary(m_domain.get_dx(), m_domain.get_dy());
    std::cout << "Dx vs Dy: " << m_domain.get_dx() << " ; " << m_domain.get_dy() << std::endl;
}

Scene::Scene(const Scene& sc){
    srand(0);
    m_tau = 1.0;
    m_timer_on = false;
    //m_fixed_connectivity = false;
    m_rt = sc.m_rt;
    m_domain = sc.m_domain;
    m_capacities = sc.m_capacities;
    m_tau = sc.m_tau;
    m_ratio = sc.m_ratio;
    m_r = sc.m_r; m_g=sc.m_g; m_b=sc.m_b;
    m_vertices= sc.m_vertices;
}

Scene* Scene::operator=(const Scene& sc){
    srand(0);
    m_tau = 1.0;
    m_timer_on = false;
    //m_fixed_connectivity = false;
    m_rt = sc.m_rt;
    m_domain = sc.m_domain;
    m_capacities = sc.m_capacities;
    m_tau = sc.m_tau;
    m_ratio = sc.m_ratio;
    m_r = sc.m_r; m_g=sc.m_g; m_b=sc.m_b;
    m_vertices= sc.m_vertices;
}

void Scene::generate_random_sites_based_on_image(const unsigned nb)
{
    if (!m_domain.is_valid()) return;
    std::vector<Point> points;
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();
    while (points.size() != nb)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        
        Point point(x, y);
        double prob  = random_double(0.0, 1.0);
        double value = m_domain.get_value(point, true) - PIXEL_EPS;
        if (prob < value) points.push_back(point);
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
    init_colors(points.size());
}

void Scene::compute_position_gradient(std::vector<Vector>& gradient, FT coef)
{
    gradient.clear();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];        
        if (vi->is_hidden()) continue;
 
        FT Vi = vi->compute_area();
        Point xi = vi->get_position();
        Point ci = vi->compute_centroid();
        Vector gi = -2.0*Vi*(xi - ci);
        gradient.push_back(coef * gi);
    }
}

FT Scene::compute_position_threshold(FT epsilon) const
{
    // reference: 1e-4 for 1000 sites
    FT A = compute_value_integral();
    unsigned N = count_visible_sites();
    return (0.1*epsilon) * (std::sqrt(A*A*A)) / FT(N);
}

void Scene::init_colors(const unsigned nb)
{
    m_r.clear();
    m_g.clear();
    m_b.clear();
    for (unsigned i = 0; i < nb; ++i)
    {
        m_r.push_back(random_double(0.0, 1.0));
        m_g.push_back(random_double(0.0, 1.0));
        m_b.push_back(random_double(0.0, 1.0));
    }
}

unsigned Scene::count_visible_sites() const
{
    unsigned nb = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        nb++;
    }
    return nb;   
}

void Scene::collect_visible_points(std::vector<Point>& points) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        points.push_back(vi->get_position());
    }
}

void Scene::collect_visible_weights(std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        weights.push_back(vi->get_weight());
    }
}

void Scene::collect_sites(std::vector<Point>& points,
                          std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        Point pi = vi->get_position();
        points.push_back(pi);
        
        FT wi = 0.0;
        wi = vi->get_weight();
        weights.push_back(wi);
    }
}

void Scene::collect_singularities(std::vector<PointSingularity> &pointSingularities,
                                  std::vector<CurveSingularity> &curveSingularities) const
{
    for (unsigned i = 0; i < m_point_singularities.size(); i++)
    {
        pointSingularities.push_back(m_point_singularities[i]);
    }

    for (unsigned i = 0; i < m_curve_singularities.size(); i++)
    {
        curveSingularities.push_back(m_curve_singularities[i]);
    }
}

void Scene::clear_triangulation()
{
    m_ratio.clear();
    m_vertices.clear();
    m_capacities.clear();
    m_rt.clear();
}

bool Scene::update_triangulation(bool skip)
{
    std::vector<FT> weights;
    std::vector<Point> points;
    collect_sites(points, weights);
    return construct_triangulation(points, weights, skip);
}

bool Scene::construct_triangulation(const std::vector<Point>& points,
                                    const std::vector<FT>& weights,
                                    bool skip)
{
    clear_triangulation();
    bool ok = populate_vertices(points, weights);

    if(!ok){
        std::cerr << "Warning, some vertices are hidden" << std::endl;
    }

    if (ok || !skip)
    {
        pre_build_dual_cells();
        assign_pixels();
        assign_singularites();
        pre_compute_area();
        //compute_capacities(m_capacities);
    }
    
    return (ok || !skip);
}

/*void Scene::update_triangulation_values()
{
    pre_build_dual_cells();
    assign_pixels();
    assign_singularites();
    pre_compute_area();
}*/

bool Scene::populate_vertices(const std::vector<Point>& points,
                              const std::vector<FT>& weights)
{    
    unsigned nb = 0;
    unsigned nsites = points.size();
    FT weight_sum = 0.0;
    for (unsigned i = 0; i < nsites; ++i)
    {
        weight_sum += weights[i];
        Vertex_handle vertex = insert_vertex(points[i], weights[i], nb);
        if (vertex == Vertex_handle())
        {
            std::cerr << "did not insert vertex with index " << nb << std::endl;
            continue;
        }
        m_vertices.push_back(vertex);
        nb++;
    }

    std::cout << "weight sum is " << weight_sum << std::endl;

    bool none_hidden = true;
    if (count_visible_sites() != m_vertices.size())
        none_hidden = false;
    
    return none_hidden;
}

Vertex_handle Scene::insert_vertex(const Point& point,
                                   const FT weight,
                                   const unsigned index)
{
    Weighted_point wp(point, weight);
    Vertex_handle vertex = m_rt.insert(wp);

    if (vertex->get_index() != -1)
    {
        std::cerr << "Did not insert point @ (" << point.x() << ", " << point.y() << ") with weight = " << weight << std::endl;
        return Vertex_handle();
    }else{
        //std::cout << "inserted point @ (" << point.x() << ", " << point.y() << ") with weight = " << weight << std::endl;
    }

    vertex->set_index(index);
    return vertex;
}

/*void Scene::delete_vertex(Vertex_handle vd)
{
    int i= -1;
    i=findIndexVerticeBySite(vd);

    if (i == -1) std::cout << "Not inside the m_vertices, can't be deleted" << std::endl;
    return;
    //m_rt.delete_vertex(vd);
    //int i=findIndexVertice(vd);
    //m_vertices.erase(m_vertices.begin()+i);

}*/

FT Scene::compute_mean_capacity() const
{
    FT domain_area = compute_value_integral();
    unsigned nb = count_visible_sites();
    return (domain_area / nb);
}

void Scene::compute_capacities(std::vector<FT>& capacities) const
{
    FT C = compute_mean_capacity();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        FT Ci = 0.0;
        Vertex_handle vi = m_vertices[i];
        if (!vi->is_hidden()) Ci = C;
        capacities.push_back(Ci);
    }
}

void Scene::update_positions(const std::vector<Point>& points, bool clamp, bool hidden)
{
    unsigned j = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;
        
        Point pi = points[j++];
        if (clamp) pi = m_domain.clamp(pi);
        vi->set_position(pi);
    }
}

void Scene::update_weights(const std::vector<FT>& weights, bool hidden)
{    
    unsigned j = 0;
    FT mean = compute_mean(weights);
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;
        //vi->set_weight(weights[j++]);
        vi->set_weight(weights[j++] - mean);
    }
}

void Scene::update_singularities(const std::vector<PointSingularity> &pointSingularities, const std::vector<CurveSingularity> &curveSingularities)
{
    m_point_singularities = pointSingularities;
    m_curve_singularities = curveSingularities;
}

void Scene::reset_weights()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->set_weight(0.0);
    }
    update_triangulation();
}

FT Scene::compute_value_integral() const
{
	return m_domain.integrate_intensity();
}

void Scene::pre_build_dual_cells()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        if (vertex->is_hidden()) continue;
        
        bool ok = m_rt.pre_build_polygon(vertex, vertex->dual().points());
        /*
        if (!ok) 
            std::cout << "Vertex " << vertex->get_index() 
            << ": pre_build_dual_cell failed" << std::endl;
        */
    }
}

void Scene::pre_compute_area()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->pre_compute_area();
    }

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->pre_compute_centroid();
    }
}

FT Scene::integrate_singularities()
{
    FT val = 0.0;

    std::vector<PointSingularity>::iterator it;

    for (it = m_point_singularities.begin(); it != m_point_singularities.end(); it++)
    {
        PointSingularity ps = *it;
        val += ps.get_value();
    }

    return val;
}

void export_cells_as_svg(std::vector<std::vector<Point>> &cells, std::vector<std::vector<double>> &colors, std::string filename) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)1 / (double)1) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    for (int i=0; i<cells.size(); i++) {
        std::vector<Point> cell = cells[i];
        std::string path_str = "M";
        for (std::size_t j = 0; j < cell.size(); ++j) {
            const auto& point = cell[j];
            path_str += std::to_string((point[0] / (double)1) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)1) * 1000.0f * ((double)1 / (double)1));

            if (j < cell.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"" << "rgb(" << colors[i][0] * 255 << ", " << colors[i][1] * 255 << ", " << colors[i][2] * 255 << ")\" stroke=\"none\" stroke-width=\"" << 5.0 << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

void Scene::draw_bounded_dual() const
{
    std::vector<std::vector<Point>> polygons;
    std::vector<std::vector<double>> colors;

    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        std::vector<Point> polygon;

        if (vi->is_hidden()) continue;
        m_rt.build_polygon(vi, polygon);

        for (int j = 0; j < polygon.size(); j++)
        {
            double new_x = (polygon[j].x() + 1.0) * 0.5;
            double new_y = (polygon[j].y() + 1.0) * 0.5;

            Point new_point(new_x, new_y);

            polygon[j] = new_point;
        }
        

        polygons.push_back(polygon);
        colors.push_back({1, 1, 1});
    }

    export_cells_as_svg(polygons, colors, "bounded_dual.svg");
}

void Scene::draw_gradient()  const
{
    std::vector<std::vector<Point>> polygons;
    std::vector<std::vector<double>> colors;

    float g_max = 0.00001;
    float g_min = -0.000001;

    for (uint i=0; i<gradient.size(); i++)
    {
        if(gradient[i] > g_max)
            g_max = gradient[i];
        else if (gradient[i] < g_min)
            g_min = gradient[i];
    }

    for (uint i=0; i<gradient.size(); i++)
    {
        if(m_vertices[i]->is_hidden()) continue;

        Vertex_handle vertex = m_vertices[i];

        double r,g,b=0.8;

        r = gradient[i] / g_min;
        b = gradient[i] / g_max;
        g = (1.0 - r - b);

        colors.push_back({r, g, b});

        for (unsigned i = 0; i < vertex->nb_pixels(); ++i)
        {
            const Pixel& pixel = vertex->get_pixel(i);
            const ConvexPolygon& shape = pixel.get_shape();
            polygons.push_back(shape.get_points());

            for (int j = 0; j < polygons[i].size(); j++)
            {
                double new_x = (polygons[i][j].x() + 1.0) * 0.5;
                double new_y = (polygons[i][j].y() + 1.0) * 0.5;

                Point new_point(new_x, new_y);

                polygons[i][j] = new_point;
            }

            colors.push_back({r, g, b});
        }
    }

    export_cells_as_svg(polygons, colors, "../gradient.svg");
}