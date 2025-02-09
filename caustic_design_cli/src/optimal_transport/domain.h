#ifndef _DOMAIN_H_
#define _DOMAIN_H_


// local
#include "util.h"
#include "pgm.h"
#include "config.h"

#include <CImg/CImg.h>

#define PIXEL_EPS 0.0 //1.0e-6

template <class Kernel>
class CDomain
{
public:
    typedef typename Kernel::FT        FT;
    typedef typename Kernel::Point_2   Point;
    typedef typename Kernel::Segment_2 Segment;
    
private:
    double m_dx;
    double m_dy;
    std::string filename;
    cimg_library::CImg<unsigned char> m_image;
    bool m_invert;
    
public:
    CDomain()
    { 
        m_dx = 0.0;
        m_dy = 0.0;
        m_invert = false;
    }
    
    void toggle_invert() { m_invert = !m_invert; }
    
    double compute_area() const { return 4.0*std::abs(m_dx)*std::abs(m_dy); }
    
    double get_dx() const { return m_dx; }
    
    double get_dy() const { return m_dy; }
    
    double get_px() const { return 2.0*get_dx()/get_width(); } // pixel width
    
    double get_py() const { return 2.0*get_dy()/get_height(); } // pixel height

    void tonemap(double /*key*/) { /*m_image.tonemap(key);*/ }

    bool load(const std::string& filename, const int width);

    double compute_x(unsigned i) const
    {
        unsigned w = get_width();
        return m_dx*(2.0*double(i)/double(w) - 1.0);
    }
    
    double compute_y(unsigned j) const
    {
        unsigned h = get_height();
        return m_dy*(1.0 - 2.0*double(j)/double(h));
    }
    
    void locate(const Point& p, 
                unsigned& i, 
                unsigned& j,
                bool pixelid = true) const
    {
        double x = 0.5 * (p.x() + m_dx) / m_dx;
        double y = 0.5 * (m_dy - p.y()) / m_dy;
        x = std::max(0.0, x);
        x = std::min(1.0, x);
        y = std::max(0.0, y);
        y = std::min(1.0, y);
        unsigned w = get_width() ;
        unsigned h = get_height();
        i = (unsigned) std::floor( x * double(w) + EPS );
        j = (unsigned) std::floor( y * double(h) + EPS );
        if (pixelid && i == w) i--;
        if (pixelid && j == h) j--;
    }
    
    bool is_grid_line(const Segment& segment) const
    {
        const Point& p0 = segment.source();
        const Point& p1 = segment.target();
        if (p0 == p1) return true;
        
        unsigned i0, j0;
        locate(p0, i0, j0, false);
        double x0 = compute_x(i0);
        double y0 = compute_y(j0);        
        if (std::abs(x0 - p0.x()) > EPS || 
            std::abs(y0 - p0.y()) > EPS) return false;

        unsigned i1, j1;
        locate(p1, i1, j1, false);        
        double x1 = compute_x(i1);
        double y1 = compute_y(j1);
        if (std::abs(x1 - p1.x()) > EPS || 
            std::abs(y1 - p1.y()) > EPS) return false;

        if (i0 == i1) return true;
        if (j0 == j1) return true;
        return false;
    }
    
    bool is_valid() const
    {
        return true;
    }
    
    unsigned get_width() const
    {
        return m_image.width();
    }
    
    unsigned get_height() const
    {
        return m_image.height();
    }
    
    double get_max_value() const
    {
        //return m_image.max_value();
        return 1.0;
    }

    std::string get_filename()
    {
        return filename;
    }
    
    double get_value(unsigned i, unsigned j, const bool normalized = false) const
    {
        double r = m_image(i, get_height() - j, 0) / 255.0; // Normalize R component
        double g = m_image(i, get_height() - j, 1) / 255.0; // Normalize G component
        double b = m_image(i, get_height() - j, 2) / 255.0; // Normalize B component
        double value = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method

        if (m_invert) value = get_max_value() - value;
        if (normalized) value /= get_max_value();
        return value + PIXEL_EPS;
    }
    
    double get_value(const Point& p, const bool normalized = false) const
    {
        unsigned i, j;
        locate(p, i, j);
        return get_value(i, j, normalized);
    }
    
    double integrate_intensity() const
    {
        double sum = 0.0;
        double pixel_area = get_px()*get_py();
        for (unsigned i = 0; i < get_width(); ++i)
        {
            for (unsigned j = 0; j < get_height(); ++j)
            {
                double value = get_value(i, j);
                sum += value*pixel_area;
            }
        }
        return sum;
    }
    
    bool is_outside(const Point& p) const
    {
        if (std::abs(p.x()) > m_dx) return true;
        if (std::abs(p.y()) > m_dy) return true;
        return false;
    }
    
    Point clamp(const Point& p) const
    {
        double x = p.x();
        if (x < -m_dx) x = -m_dx;
        if (x >  m_dx) x =  m_dx;
        
        double y = p.y();
        if (y < -m_dy) y = -m_dy;
        if (y >  m_dy) y =  m_dy;
        
        return Point(x, y);
    }
};

template <class Kernel>
bool CDomain<Kernel>::load(const std::string& filename, const int width)
{
    m_image.load(filename.c_str());

    std::cout << m_image.width() << std::endl;

    if(width != 0)
    {
        double scale = (double)width / double(m_image.width());
        int height = (int)(((double)m_image.height()) * scale);
        m_image.resize(width, height);
    }

    this->filename = filename;

    unsigned w = get_width();
    unsigned h = get_height();
    m_dx = DOMAIN_WIDTH;
    m_dy = m_dx * double(h) / double(w);
    return true;
}

#endif
