#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "utils.h"
#include "frame.h"
#include "common.h"
#include "dwt.h"

struct view {
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
};

int is_max(int *this, size_t stride, size_t half_size)
{
	int *other;
	int *begin = this - half_size*stride;
	int *end = this + half_size*stride;

#if 1
	/* suppose *this > others */
	for (other = begin; other <= end; other += stride) {
		if (*other > *this)
			goto next;
	}

	return +1;

	next:

	/* suppose *this < others */
	for (other = begin; other <= end; other += stride) {
		if (*other < *this)
			return 0;
	}

	return -1;
#else
	/* suppose |*this| > |others| */
	for (other = begin; other <= end; other += stride) {
		if (abs_(*other) > abs_(*this))
			return 0;
	}

	return 1;
#endif
}

int is_double_max(int *this, size_t stride, size_t half_size, size_t block_distance)
{
	int max0 = is_max(this, stride, half_size);
	int max1 = is_max(this+block_distance*stride, stride, half_size);

	return max0 && max1;
}

void filter1(int *this, size_t stride, size_t half_size, size_t block_distance, int threshold)
{
#if 1
	if (abs_(*this) > threshold)
		return;
#endif

	if (is_double_max(this, stride, half_size, block_distance))
		*this = 0;
}

void filter2(int *this, size_t stride_x, size_t stride_y, size_t half_size, size_t block_distance, int threshold)
{
#if 1
	if (abs_(*this) > threshold)
		return;
#endif

	if (is_double_max(this, stride_x, half_size, block_distance) &&
	    is_double_max(this, stride_y, half_size, block_distance)) {
		*this = 0;
	}
}

void denoise_hl2band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0; y < height; ++y) {
		for (x = 0+1+block_distance; x < width-1-block_distance; ++x) {
			filter1(data + y*stride_y + x*stride_x, stride_x, 1, block_distance, threshold);
		}
	}
}

void denoise_hh2band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0; y < height; ++y) {
		for (x = 0+1+block_distance; x < width-1-block_distance; ++x) {
			filter2(data + y*stride_y + x*stride_x, stride_x, stride_x, 1, block_distance, threshold);
		}
	}
}

void denoise_hh1band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0; y < height; ++y) {
		for (x = 0+3+block_distance; x < width-3-block_distance; ++x) {
			filter2(data + y*stride_y + x*stride_x, stride_x, stride_x, 3, block_distance, threshold);
		}
	}
}

void denoise_lh2band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0+1+block_distance; y < height-1-block_distance; ++y) {
		for (x = 0; x < width; ++x) {
			filter1(data + y*stride_y + x*stride_x, stride_y, 1, block_distance, threshold);
		}
	}
}

void denoise_hl1band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0; y < height; ++y) {
		for (x = 0+3+block_distance; x < width-3-block_distance; ++x) {
			filter1(data + y*stride_y + x*stride_x, stride_x, 3, block_distance, threshold);
		}
	}
}

void denoise_lh1band(struct view *view, size_t block_distance, int threshold)
{
	int *data;
	size_t stride_x, stride_y;
	size_t width, height;
	size_t x, y;

	assert(view);

	stride_x = view->stride_x;
	stride_y = view->stride_y;
	width = view->width;
	height = view->height;
	data = view->data;

	for (y = 0+3+block_distance; y < height-3-block_distance; ++y) {
		for (x = 0; x < width; ++x) {
			filter1(data + y*stride_y + x*stride_x, stride_y, 3, block_distance, threshold);
		}
	}
}

void denoise(struct frame *frame)
{
	size_t width, height;
	int *data;
	int i;

	assert(frame);

	width = ceil_multiple8(frame->width);
	height = ceil_multiple8(frame->height);

	data = frame->data;

	for (i = 0; i < 1; ++i) {
		struct view view;

		printf("pass %i\n", i);

		/* HL2 */
		view.data = data + (0)*width + (2);
		view.stride_x = 4;
		view.stride_y = 4*width;
		view.width = width/4;
		view.height = height/4;
		denoise_hl2band(&view, 2, 32);

		/* LH2 */
		view.data = data + (2)*width + (0);
		view.stride_x = 4;
		view.stride_y = 4*width;
		view.width = width/4;
		view.height = height/4;
		denoise_lh2band(&view, 2, 32);

		/* HH2 */
		view.data = data + (2)*width + (2);
		view.stride_x = 4;
		view.stride_y = 4*width;
		view.width = width/4;
		view.height = height/4;
		denoise_hh2band(&view, 2, 16);

		/* HL1 */
		view.data = data + (0)*width + (1);
		view.stride_x = 2;
		view.stride_y = 2*width;
		view.width = width/2;
		view.height = height/2;
		denoise_hl1band(&view, 4, 16);

		/* LH1 */
		view.data = data + (1)*width + (0);
		view.stride_x = 2;
		view.stride_y = 2*width;
		view.width = width/2;
		view.height = height/2;
		denoise_lh1band(&view, 4, 16);

		/* TODO HH1 */
		view.data = data + (1)*width + (1);
		view.stride_x = 2;
		view.stride_y = 2*width;
		view.width = width/2;
		view.height = height/2;
		denoise_hh1band(&view, 4, 16);
	}
}

int main(int argc, char *argv[])
{
	struct frame frame, original;
	struct parameters parameters;
	float mse;
	float psnr;

	assert( ~-1 == 0 );
	assert( -1 >> 1 == -1 );

	if (argc < 3) {
		fprintf(stderr, "[ERROR] argument expected\n");
		return EXIT_FAILURE;
	}

	/** (1) load input image */

	dprint (("[DEBUG] loading...\n"));

	if ( frame_load_pgm(&frame, argv[1]) ) {
		fprintf(stderr, "[ERROR] unable to load image\n");
		return EXIT_FAILURE;
	}

	if ( frame_load_pgm(&original, argv[2]) ) {
		fprintf(stderr, "[ERROR] unable to load image\n");
		return EXIT_FAILURE;
	}

	if ( frame.width > (1<<20) || frame.width < 17 ) {
		fprintf(stderr, "[ERROR] unsupported image width\n");
		return EXIT_FAILURE;
	}

	if ( frame.height < 17 ) {
		fprintf(stderr, "[ERROR] unsupported image height\n");
		return EXIT_FAILURE;
	}

	frame_save_pgm(&frame, "input.pgm");

	parameters.DWTtype = 0;

	dprint (("[DEBUG] transform...\n"));

	/** (2) forward DWT */
	if (dwt_encode(&frame, &parameters)) {
		fprintf(stderr, "[ERROR] transform failed\n");
		return EXIT_FAILURE;
	}

	dprint (("[DEBUG] dump...\n"));

	frame_dump_chunked_as_semiplanar(&frame, "dwt3.pgm", 1);
#if 1
	denoise(&frame);
#endif
	frame_dump_chunked_as_semiplanar(&frame, "dwt3-denoised.pgm", 1);

	dprint (("[DEBUG] inverse transform...\n"));

	/** (2) inverse DWT */
	if (dwt_decode(&frame, &parameters)) {
		fprintf(stderr, "[ERROR] inverse transform failed\n");
		return EXIT_FAILURE;
	}

	frame_dump(&frame, "decoded.pgm", 1);

	mse = frame_get_mse(&frame, &original, frame.bpp);
	psnr = convert_mse_to_psnr(mse, frame.bpp);

	fprintf(stderr, "[DEBUG] psnr = %f\n", psnr);

	/** (1) save output image */

	dprint (("[DEBUG] saving...\n"));

	if ( frame_save_pgm(&frame, "output.pgm") ) {
		fprintf(stderr, "[ERROR] unable to save an output raster\n");
		return EXIT_FAILURE;
	}

	/** (1) release resources */

	frame_destroy(&frame);

	frame_destroy(&original);

	return EXIT_SUCCESS;
}
