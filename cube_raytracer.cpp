#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>

// Implementación personalizada de clamp ya que std::clamp
// requiere C++17 y no tengo esa versión en este compu ni en los de Turing

template <typename T>
inline T clamp(T x, T min, T max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

// Constantes y utilidades basicas del programa
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Clase Vector3 para representar puntos, vectores y colores de la lectura
class Vec3
{
public:
    double e[3];

    Vec3() : e{0, 0, 0} {}
    Vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double &operator[](int i) { return e[i]; }

    Vec3 &operator+=(const Vec3 &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vec3 &operator*=(const double t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vec3 &operator/=(const double t)
    {
        return *this *= 1 / t;
    }

    double length() const
    {
        return sqrt(length_squared());
    }

    double length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }
};

// Alias para Vec3
using Point3 = Vec3; // Punto 3D
using Color = Vec3;  // Color RGB

// Funciones de utilidad para Vec3
inline std::ostream &operator<<(std::ostream &out, const Vec3 &v)
{
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline Vec3 operator+(const Vec3 &u, const Vec3 &v)
{
    return Vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline Vec3 operator-(const Vec3 &u, const Vec3 &v)
{
    return Vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline Vec3 operator*(const Vec3 &u, const Vec3 &v)
{
    return Vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline Vec3 operator*(double t, const Vec3 &v)
{
    return Vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline Vec3 operator*(const Vec3 &v, double t)
{
    return t * v;
}

inline Vec3 operator/(Vec3 v, double t)
{
    return (1 / t) * v;
}

inline double dot(const Vec3 &u, const Vec3 &v)
{
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

inline Vec3 cross(const Vec3 &u, const Vec3 &v)
{
    return Vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline Vec3 unit_vector(Vec3 v)
{
    return v / v.length();
}

// Clase Ray para representar rayos
class Ray
{
public:
    Point3 orig;
    Vec3 dir;

    Ray() {}
    Ray(const Point3 &origin, const Vec3 &direction) : orig(origin), dir(direction) {}

    Point3 origin() const { return orig; }
    Vec3 direction() const { return dir; }

    Point3 at(double t) const
    {
        return orig + t * dir;
    }
};

// Estructura para almacenar información de intersección
struct HitRecord
{
    Point3 p;
    Vec3 normal;
    double t;
    bool front_face;

    inline void set_face_normal(const Ray &r, const Vec3 &outward_normal)
    {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

// Clase abstracta para objetos que pueden ser golpeados por un rayo
class Hittable
{
public:
    virtual bool hit(const Ray &r, double t_min, double t_max, HitRecord &rec) const = 0;
};

// Clase para representar un cubo
class Cube : public Hittable
{
public:
    Point3 center;
    double side_length;
    Color color;

    Cube(Point3 cen, double side, Color col) : center(cen), side_length(side), color(col) {}

    virtual bool hit(const Ray &r, double t_min, double t_max, HitRecord &rec) const override
    {
        // Calculamos los límites del cubo
        double half_side = side_length / 2.0;
        Point3 min_point = Point3(center.x() - half_side, center.y() - half_side, center.z() - half_side);
        Point3 max_point = Point3(center.x() + half_side, center.y() + half_side, center.z() + half_side);

        // Algoritmo de intersección rayo-caja (AABB)
        double tmin, tmax, tymin, tymax, tzmin, tzmax;

        double invD = 1.0 / r.direction().x();
        if (invD >= 0)
        {
            tmin = (min_point.x() - r.origin().x()) * invD;
            tmax = (max_point.x() - r.origin().x()) * invD;
        }
        else
        {
            tmin = (max_point.x() - r.origin().x()) * invD;
            tmax = (min_point.x() - r.origin().x()) * invD;
        }

        invD = 1.0 / r.direction().y();
        if (invD >= 0)
        {
            tymin = (min_point.y() - r.origin().y()) * invD;
            tymax = (max_point.y() - r.origin().y()) * invD;
        }
        else
        {
            tymin = (max_point.y() - r.origin().y()) * invD;
            tymax = (min_point.y() - r.origin().y()) * invD;
        }

        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        invD = 1.0 / r.direction().z();
        if (invD >= 0)
        {
            tzmin = (min_point.z() - r.origin().z()) * invD;
            tzmax = (max_point.z() - r.origin().z()) * invD;
        }
        else
        {
            tzmin = (max_point.z() - r.origin().z()) * invD;
            tzmax = (min_point.z() - r.origin().z()) * invD;
        }

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        // Verificamos si la intersección está dentro del rango válido
        if (tmin < t_max && tmin > t_min)
        {
            rec.t = tmin;
            rec.p = r.at(rec.t);

            // Calculamos la normal en el punto de intersección
            Vec3 outward_normal;
            Point3 relative_p = rec.p - center;
            double epsilon = 1e-8;

            if (std::abs(relative_p.x() - half_side) < epsilon)
                outward_normal = Vec3(1, 0, 0);
            else if (std::abs(relative_p.x() + half_side) < epsilon)
                outward_normal = Vec3(-1, 0, 0);
            else if (std::abs(relative_p.y() - half_side) < epsilon)
                outward_normal = Vec3(0, 1, 0);
            else if (std::abs(relative_p.y() + half_side) < epsilon)
                outward_normal = Vec3(0, -1, 0);
            else if (std::abs(relative_p.z() - half_side) < epsilon)
                outward_normal = Vec3(0, 0, 1);
            else
                outward_normal = Vec3(0, 0, -1);

            rec.set_face_normal(r, outward_normal);
            return true;
        }

        return false;
    }
};

// Clase para representar un plano (suelo)
class Plane : public Hittable
{
public:
    Point3 point;
    Vec3 normal;
    Color color;

    Plane(Point3 p, Vec3 n, Color col) : point(p), normal(unit_vector(n)), color(col) {}

    virtual bool hit(const Ray &r, double t_min, double t_max, HitRecord &rec) const override
    {
        double denom = dot(normal, r.direction());

        // Si el rayo es paralelo al plano
        if (std::abs(denom) < 1e-8)
            return false;

        double t = dot(point - r.origin(), normal) / denom;

        if (t < t_min || t > t_max)
            return false;

        rec.t = t;
        rec.p = r.at(t);
        rec.set_face_normal(r, normal);

        return true;
    }
};

// Lista de objetos que pueden ser golpeados por un rayo
class HittableList : public Hittable
{
public:
    std::vector<Hittable *> objects;

    HittableList() {}
    HittableList(Hittable *object) { add(object); }

    void clear() { objects.clear(); }
    void add(Hittable *object) { objects.push_back(object); }

    virtual bool hit(const Ray &r, double t_min, double t_max, HitRecord &rec) const override
    {
        HitRecord temp_rec;
        bool hit_anything = false;
        auto closest_so_far = t_max;

        for (const auto &object : objects)
        {
            if (object->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        return hit_anything;
    }
};

// Clase para la cámara
class Camera
{
public:
    Point3 origin;
    Point3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;

    Camera()
    {
        auto aspect_ratio = 1;
        auto viewport_height = 2.0;
        auto viewport_width = aspect_ratio * viewport_height;
        auto focal_length = 1.0;

        origin = Point3(0, 0, 0);
        horizontal = Vec3(viewport_width, 0, 0);
        vertical = Vec3(0, viewport_height, 0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3(0, 0, focal_length);
    }

    Ray get_ray(double u, double v) const
    {
        return Ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
    }
};

// Función para escribir el color en el archivo de imagen
void write_color(std::ostream &out, Color pixel_color, int samples_per_pixel)
{
    // Dividimos el color acumulado por el número de muestras y aplicamos corrección gamma
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Dividimos el color por el número de muestras
    auto scale = 1.0 / samples_per_pixel;
    r *= scale;
    g *= scale;
    b *= scale;

    // Aplicamos corrección gamma (gamma 2) con una raíz cuadrada simple
    r = sqrt(r);
    g = sqrt(g);
    b = sqrt(b);

    // Escribimos los valores [0,255] para RGB con clamp
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

// Función para generar un vector aleatorio en una esfera unitaria
Vec3 random_in_unit_sphere()
{
    static std::mt19937 generator;
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);

    while (true)
    {
        auto p = Vec3(distribution(generator), distribution(generator), distribution(generator)) * 2.0 - Vec3(1, 1, 1);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

// Función para generar un vector unitario aleatorio
Vec3 random_unit_vector()
{
    return unit_vector(random_in_unit_sphere());
}

// Función para calcular el color de un rayo con profundidad de rebote
Color ray_color(const Ray &r, int depth, const Hittable &world)
{
    // Si hemos excedido el límite de rebotes, no se recoge más luz
    if (depth <= 0)
        return Color(0, 0, 0);

    HitRecord rec;

    if (world.hit(r, 0.001, infinity, rec))
    {
        // Calculamos la dirección de rebote difuso
        Vec3 target = rec.p + rec.normal + random_unit_vector();
        Ray scattered(rec.p, target - rec.p);

        // Buscamos el objeto que fue golpeado para obtener su color
        Color object_color(0.8, 0.8, 0.8); // Color por defecto gris claro

        if (const HittableList *list = dynamic_cast<const HittableList *>(&world))
        {
            // Encontramos qué objeto fue golpeado para obtener su color
            for (const auto &object : list->objects)
            {
                // Evitamos la recursión infinita
                if (object == const_cast<Hittable *>(&world))
                    continue;

                // Verificamos si el punto de intersección está en este objeto
                HitRecord temp_rec;
                if (object->hit(Ray(rec.p - 0.001 * rec.normal, rec.normal), 0.001, 0.1, temp_rec))
                {
                    // Obtenemos el color según el tipo de objeto
                    if (Cube *cube = dynamic_cast<Cube *>(object))
                    {
                        object_color = cube->color;
                        break;
                    }
                    else if (Plane *plane = dynamic_cast<Plane *>(object))
                    {
                        object_color = plane->color;
                        break;
                    }
                }
            }
        }

        // Aplicamos el color del objeto al rebote difuso
        return 0.5 * ray_color(scattered, depth - 1, world) * object_color;
    }

    // Si no golpeamos nada, devolvemos el color del cielo (gradiente)
    Vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * Color(1.0, 1.0, 1.0) + t * Color(0.5, 0.7, 1.0);
}

int main()
{
    // Configuración de la imagen
    const auto aspect_ratio = 1;
    const int image_width = 600;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 200;
    const int max_depth = 200;

    // Configuración de la cámara - Posicionada para ver mejor el cubo rotado
    Camera cam;
    cam.origin = Point3(2, 2, 4); // Posición diagonal para ver mejor las caras del cubo
    cam.horizontal = Vec3(4, 0, 0);
    cam.vertical = Vec3(0, 4 / aspect_ratio, 0);
    // Ajustamos la dirección para mirar hacia el centro del cubo
    Vec3 look_at = Point3(0.5, 0, -0.5); // Centro del cubo
    Vec3 view_direction = unit_vector(look_at - cam.origin);
    cam.lower_left_corner = cam.origin - cam.horizontal / 2 - cam.vertical / 2 + view_direction;

    // Creamos la escena
    HittableList world;

    world.add(new Cube(Point3(0.5, 0, -0.5), 2.0, Color(0.3, 0.3, 0.5)));

    // Añadimos un plano como suelo
    world.add(new Plane(Point3(0, -1, 0), Vec3(0, 1, 0), Color(0.3, 0.3, 0.5)));

    // Generador de números aleatorios para el anti-aliasing
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Renderizamos la escena
    std::ofstream output("cube_render.ppm");
    output << "P3\n"
           << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            Color pixel_color(0, 0, 0);
            // Realizamos múltiples muestras por píxel para anti-aliasing
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + distribution(generator)) / (image_width - 1);
                auto v = (j + distribution(generator)) / (image_height - 1);
                Ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, max_depth, world);
            }
            write_color(output, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
    return 0;
}