#include <stdio.h>
#if ((__STDC_VERSION__ >= 199901L) || (_MSC_VER))
	#include <stdint.h>
#endif
#include <string.h>
#include <errno.h>

#include <png.h>

#include "OpenSimplex2F.h"

#define WIDTH 512
#define HEIGHT 512
#define FEATURE_SIZE 24

static int write_png_image(const char *filename, unsigned char *pixels, int w, int h, int has_alpha)
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_byte **row;
	int x, y, rc, colordepth = 8;
	int bytes_per_pixel = has_alpha ? 4 : 3;
	FILE *f;

	f = fopen(filename, "w");
	if (!f) {
		fprintf(stderr, "fopen: %s:%s\n", filename, strerror(errno));
		return -1;
	}
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr)
		goto cleanup1;
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		goto cleanup2;
	if (setjmp(png_jmpbuf(png_ptr))) /* oh libpng, you're old as dirt, aren't you. */
		goto cleanup2;

	png_set_IHDR(png_ptr, info_ptr, (size_t) w, (size_t) h, colordepth,
			has_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
			PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
			PNG_FILTER_TYPE_DEFAULT);

	row = png_malloc(png_ptr, h * sizeof(*row));
	for (y = 0; y < h; y++) {
		row[y] = png_malloc(png_ptr, w * bytes_per_pixel);
		for (x = 0; x < w; x++) {
			unsigned char *r = (unsigned char *) row[y];
			unsigned char *src = (unsigned char *)
				&pixels[y * w * bytes_per_pixel + x * bytes_per_pixel];
			unsigned char *dest = &r[x * bytes_per_pixel];
			memcpy(dest, src, bytes_per_pixel);
		}
	}

	png_init_io(png_ptr, f);
	png_set_rows(png_ptr, info_ptr, row);
	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_PACKING, NULL);

	for (y = 0; y < h; y++)
		png_free(png_ptr, row[y]);
	png_free(png_ptr, row);
	rc = 0;
cleanup2:
	png_destroy_write_struct(&png_ptr, &info_ptr);
cleanup1:
	fclose(f);
	return rc;
}

int main(__attribute__((unused)) int argc, __attribute__((unused)) char *argv[])
{
	int x, y;
	double value;
	double v0, v1, v2; /* values from different octaves. */
	uint32_t rgb;
	uint32_t image2d[HEIGHT][WIDTH];
	uint32_t image3d[HEIGHT][WIDTH];
	uint32_t image4d[HEIGHT][WIDTH];
	struct OpenSimplex2F_context *ctx;

	OpenSimplex2F(77374, &ctx);

	for (y = 0; y < HEIGHT; y++) {
		for (x = 0; x < WIDTH; x++) {
#if defined(SINGLE_OCTAVE)
			value = open_simplex_noise4(ctx, (double) x / FEATURE_SIZE, (double) y / FEATURE_SIZE, 0.0, 0.0);
#else
			/* Use three octaves: frequency N, N/2 and N/4 with relative amplitudes 4:2:1. */
			v0 = OpenSimplex2F_noise4_Classic(ctx, (double) x / FEATURE_SIZE / 4,
						(double) y / FEATURE_SIZE / 4, 0.0, 0.0);
			v1 = OpenSimplex2F_noise4_Classic(ctx, (double) x / FEATURE_SIZE / 2,
						(double) y / FEATURE_SIZE / 2, 0.0, 0.0);
			v2 = OpenSimplex2F_noise4_Classic(ctx, (double) x / FEATURE_SIZE / 1,
						(double) y / FEATURE_SIZE / 1, 0.0, 0.0);
			value = v0 * 4 / 7.0 + v1 * 2 / 7.0 + v2 * 1 / 7.0;
#endif
			rgb = 0x010101 * (uint32_t) ((value + 1) * 127.5);
			image2d[y][x] = ((uint32_t) 0x0ff << 24) | (rgb);

			value = OpenSimplex2F_noise2(ctx, (double) x / FEATURE_SIZE, (double) y / FEATURE_SIZE);
			rgb = 0x010101 * (uint32_t) ((value + 1) * 127.5);
			image3d[y][x] = ((uint32_t) 0x0ff << 24) | (rgb);

			value = OpenSimplex2F_noise3_Classic(ctx, (double) x / FEATURE_SIZE, (double) y / FEATURE_SIZE, 0.0);
			rgb = 0x010101 * (uint32_t) ((value + 1) * 127.5);
			image4d[y][x] = ((uint32_t) 0x0ff << 24) | (rgb);
		}
	}
	write_png_image("test2d.png", (unsigned char *) image2d, WIDTH, HEIGHT, 1);
	write_png_image("test3d.png", (unsigned char *) image3d, WIDTH, HEIGHT, 1);
	write_png_image("test4d.png", (unsigned char *) image4d, WIDTH, HEIGHT, 1);
	OpenSimplex2F_free(ctx);
	OpenSimplex2F_shutdown();
	return 0;
}

