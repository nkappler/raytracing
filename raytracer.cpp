#include "rtweekend.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

#include <iostream>
#include <chrono>
#include <ctime>
#include <thread>
#include <future>

color ray_color(const ray &r, const hittable &world, int depth)
{
	hit_record rec;

	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth < 0)
	{
		return color(0, 0, 0);
	}

	if (world.hit(r, 0.001, infinity, rec))
	{
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(0, 0, 0);
	}
	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

hittable_list random_scene()
{
	hittable_list world;

	auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
	world.add(make_shared<sphere>(point3(0, -10000, 0), 10000, ground_material));

	for (int a = -11; a < 11; a++)
	{
		for (int b = -11; b < 11; b++)
		{
			auto choose_mat = random_double();
			point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

			if ((center - point3(4, 0.2, 0)).length() > 0.9)
			{
				shared_ptr<material> sphere_material;

				if (choose_mat < 0.8)
				{
					// diffuse
					auto albedo = color::random() * color::random();
					sphere_material = make_shared<lambertian>(albedo);
					world.add(make_shared<sphere>(center, 0.2, sphere_material));
				}
				else if (choose_mat < 0.95)
				{
					// metal
					auto albedo = color::random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<sphere>(center, 0.2, sphere_material));
				}
				else
				{
					// glass
					sphere_material = make_shared<dielectric>(1.5);
					world.add(make_shared<sphere>(center, 0.2, sphere_material));
				}
			}
		}
	}

	auto material1 = make_shared<dielectric>(1.5);
	world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

	auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
	world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

	auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
	world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

	return world;
}

void getPixelColor(const int samples_per_pixel, int i, const int image_width, int j, const int image_height, camera &cam, color &pixel_color, hittable_list &world, const int max_bounces)
{
	for (int s = 0; s < samples_per_pixel; s++)
	{
		auto u = (i + random_double()) / (image_width - 1);
		auto v = (j + random_double()) / (image_height - 1);
		ray r = cam.get_ray(u, v);
		pixel_color += ray_color(r, world, max_bounces);
	}
}

auto renderLines(const int startLine, const int endLine, const int image_height, const int image_width, const int samples_per_pixel, const int max_bounces, camera cam, hittable_list world)
{
	std::vector<color> image(image_width * (endLine - startLine), color(0, 0, 0));

	for (int y = endLine - 1; y >= startLine; --y)
	{
		for (int x = image_width - 1; x >= 0; --x)
		{
			getPixelColor(samples_per_pixel, x, image_width, y, image_height, cam, image[(y - startLine) * image_width + (image_width - x - 1)], world, max_bounces);
		}
	}

	return image;
}

int main()
{
	auto start = std::chrono::system_clock::now();

	// Image
	const auto aspect_ratio = 3.0 / 2.0;
	const int image_width = 1200;
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	const int threadCount = 16;
	const int samples_per_pixel = 100;
	const int max_bounces = 6;

	// World
	auto world = random_scene();

	// Camera
	point3 lookfrom(13, 2, 3);
	point3 lookat(0, 0, 0);
	vec3 vup(0, 1, 0);
	auto dist_to_focus = 10.0;
	auto aperture = 0.05;

	camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

	// Render multithreaded
	const int linesPerThread = image_height / threadCount;
	std::vector<std::future<std::vector<color>>> threads(threadCount);
	for (int t = 0; t < threadCount - 1; t++)
	{
		threads[t] = std::async(renderLines, t * linesPerThread, (t + 1) * linesPerThread, image_height, image_width, samples_per_pixel, max_bounces, cam, world);
	}
	threads[threadCount - 1] = std::async(renderLines, (threadCount - 1) * linesPerThread, image_height, image_height, image_width, samples_per_pixel, max_bounces, cam, world);

	// Output
	std::cout << "P3\n"
			  << image_width << " " << image_height << "\n255\n";
	for (int t = threadCount - 1; t >= 0; --t)
	{
		auto imagePart = threads[t].get();

		for (int i = imagePart.size() - 1; i >= 0; --i)
		{
			write_color(std::cout, imagePart[i], samples_per_pixel);
		}
	}

	std::cerr << "\nDone.\n";

	// Statistics
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cerr << "finished computation at " << std::ctime(&end_time)
			  << "elapsed time: " << elapsed_seconds.count() << "s"
			  << std::endl;
}
