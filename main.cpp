#include <iostream>
#include <cmath>
#include <vector>
#include <random>

struct Vector3
{
    double x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3(double s) : x(s), y(s), z(s) {}

    friend Vector3 operator*(double t, const Vector3 &v)
    {
        return Vector3(v.x * t, v.y * t, v.z * t);
    }

    Vector3 operator*(const Vector3 &v) const
    {
        return Vector3(x * v.x, y * v.y, z * v.z);
    }

    Vector3 operator+(const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
    Vector3 operator-(const Vector3 &v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
    Vector3 operator*(double t) const { return Vector3(x * t, y * t, z * t); }
    Vector3 operator/(double t) const { return Vector3(x / t, y / t, z / t); }

    double dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }
    Vector3 cross(const Vector3 &v) const
    {
        return Vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    double length() const { return sqrt(x * x + y * y + z * z); }
    Vector3 normalized() const { return *this / length(); }
};

struct Ray
{
    Vector3 origin;
    Vector3 direction;
    Ray(Vector3 o, Vector3 d) : origin(o), direction(d) {}
};

struct Material
{
    Vector3 albedo;
    double specular;
    Material(Vector3 a, double s) : albedo(a), specular(s) {}
};

struct Hittable
{
    virtual bool hit(const Ray &ray, double t_min, double t_max, double &t, Vector3 &normal) const = 0;
    Material material;
    Hittable(Material m) : material(m) {}
    virtual ~Hittable() = default;
};

struct Sphere : Hittable
{
    Vector3 center;
    double radius;
    Sphere(Vector3 c, double r, Material m) : Hittable(m), center(c), radius(r) {}

    bool hit(const Ray &ray, double t_min, double t_max, double &t, Vector3 &normal) const override
    {
        Vector3 oc = ray.origin - center;
        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;

        if (discriminant > 0)
        {
            double sqrtd = sqrt(discriminant);
            double root = (-b - sqrtd) / (2.0 * a);
            if (root < t_max && root > t_min)
            {
                t = root;
                return true;
            }
            root = (-b + sqrtd) / (2.0 * a);
            if (root < t_max && root > t_min)
            {
                t = root;
                return true;
            }
        }
        return false;
    }
};

struct Cube : Hittable
{
    Vector3 min;
    Vector3 max;
    Cube(Vector3 mn, Vector3 mx, Material m) : Hittable(m), min(mn), max(mx) {}

    bool hit(const Ray &ray, double t_min, double t_max, double &t, Vector3 &normal) const override
    {
        Vector3 inv_dir(1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z);
        Vector3 t0 = (min - ray.origin) * inv_dir;
        Vector3 t1 = (max - ray.origin) * inv_dir;

        Vector3 tmin = Vector3(std::min(t0.x, t1.x), std::min(t0.y, t1.y), std::min(t0.z, t1.z));
        Vector3 tmax = Vector3(std::max(t0.x, t1.x), std::max(t0.y, t1.y), std::max(t0.z, t1.z));

        double t_near = std::max(std::max(tmin.x, tmin.y), tmin.z);
        double t_far = std::min(std::min(tmax.x, tmax.y), tmax.z);

        if (t_near > t_far || t_far < 0)
            return false;

        t = t_near;
        Vector3 hit_point = ray.origin + ray.direction * t;

        // Determine which face was hit based on t_near
        Vector3 normal_dir;
        if (t_near == tmin.x)
        {
            normal_dir = Vector3(ray.direction.x > 0 ? -1 : 1, 0, 0);
        }
        else if (t_near == tmin.y)
        {
            normal_dir = Vector3(0, ray.direction.y > 0 ? -1 : 1, 0);
        }
        else
        {
            normal_dir = Vector3(0, 0, ray.direction.z > 0 ? -1 : 1);
        }
        normal = normal_dir.normalized();
        Vector3 delta = hit_point - (min + max) * 0.5;
        Vector3 abs_delta(fabs(delta.x), fabs(delta.y), fabs(delta.z));

        if (abs_delta.x > abs_delta.y && abs_delta.x > abs_delta.z)
            normal = Vector3(delta.x > 0 ? 1 : -1, 0, 0);
        else if (abs_delta.y > abs_delta.z)
            normal = Vector3(0, delta.y > 0 ? 1 : -1, 0);
        else
            normal = Vector3(0, 0, delta.z > 0 ? 1 : -1);
        return true;
    }
};

Vector3 color(const Ray &ray, const std::vector<Hittable *> &objects)
{
    double t;
    double closest = std::numeric_limits<double>::max();
    bool hit_anything = false;
    Vector3 final_normal;
    const Hittable *hit_obj = nullptr;

    for (const auto &obj : objects)
    {
        Vector3 tmp_normal;
        if (obj->hit(ray, 0.001, closest, t, tmp_normal))
        {
            hit_anything = true;
            closest = t;
            final_normal = tmp_normal;
            hit_obj = obj;
            hit_obj = obj;
        }
    }

    if (hit_anything)
    {
        Vector3 hit_point = ray.origin + ray.direction * closest;
        Vector3 light_dir = Vector3(0.5, 1.0, 0.3).normalized();
        double diffuse = std::max(0.0, final_normal.dot(light_dir));
        Vector3 view_dir = ray.direction.normalized();
        Vector3 reflect_dir = light_dir - final_normal * 2.0 * final_normal.dot(light_dir);
        double specular = pow(std::max(0.0, view_dir.dot(reflect_dir)), hit_obj->material.specular);
        return hit_obj->material.albedo * (diffuse + 0.1) + Vector3(1.0) * specular;
    }

    // Sky gradient
    Vector3 unit_dir = ray.direction.normalized();
    t = 0.5 * (unit_dir.y + 1.0);
    return Vector3(1.0, 1.0, 1.0) * (1.0 - t) + Vector3(0.2, 0.5, 1.0) * t;
}

int main()
{
    const int width = 800;
    const int height = 600;
    const int samples = 50;

    std::vector<Hittable *> objects;
    Material sphere_mat(Vector3(0.8, 0.3, 0.2), 128.0);
    Material ground_mat(Vector3(0.2, 0.6, 0.2), 16.0);
    Material cube_mat(Vector3(0.2, 0.3, 0.8), 256.0);

    objects.push_back(new Sphere(Vector3(0, 0, -1), 0.5, sphere_mat));
    objects.push_back(new Sphere(Vector3(0, -100.5, -1), 100, ground_mat));
    objects.push_back(new Cube(Vector3(-1, -0.5, -2), Vector3(1, 1.5, -1.5), cube_mat));

    Vector3 lower_left(-1.5, -1.0, -1.0);
    Vector3 horizontal(3.0, 0.0, 0.0);
    Vector3 vertical(0.0, 2.0, 0.0);
    Vector3 origin(0.0, 0.0, 2.0);

    std::cout << "P3\n"
              << width << " " << height << "\n255\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int y = height - 1; y >= 0; --y)
    {
        for (int x = 0; x < width; ++x)
        {
            Vector3 col;
            for (int s = 0; s < samples; ++s)
            {
                double u = double(x + dis(gen)) / double(width);
                double v = double(y + dis(gen)) / double(height);
                Ray ray(origin, lower_left + horizontal * u + vertical * v);
                col = col + color(ray, objects);
            }
            col = col / samples;

            int ir = int(255.99 * col.x);
            int ig = int(255.99 * col.y);
            int ib = int(255.99 * col.z);

            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }
    return 0;
}