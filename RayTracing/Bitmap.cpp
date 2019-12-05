#include "Bitmap.h"

/* This class is based of an article from Solarian Programmer */
// Source: https://solarianprogrammer.com/2018/11/19/cpp-reading-writing-bmp-images/

Bitmap::Bitmap(int32_t width, int32_t height, bool has_alpha) 
{
	if (width <= 0 || height <= 0) {
		throw std::runtime_error("Error::Bitmap::Constructor : Width and height must be positive integers!");
	}
	// Set dimensions in the info header
	bmp_info_header.width = width;
	bmp_info_header.height = height;

	// If has_alpha is true, we make space for 32 bit color image
	if (has_alpha) {
		bmp_info_header.size = sizeof(BitmapInfoHeader) + sizeof(BitmapColorHeader);
		file_header.offset_data = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) + sizeof(BitmapColorHeader);

		bmp_info_header.bit_count = 32;
		bmp_info_header.compression = 3;
		_row_stride = width * 4;

		auto capacity = _row_stride * height;
		image_data.resize(capacity);
		file_header.file_size = file_header.offset_data + capacity;
	}
	else { // 24 bit color
		bmp_info_header.size = sizeof(BitmapInfoHeader);
		file_header.offset_data = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);

		bmp_info_header.bit_count = 24;
		bmp_info_header.compression = 0;
		_row_stride = width * 3;

		auto capacity = _row_stride * height;
		image_data.resize(capacity);
		uint32_t new_stride = alignStride(4);
		file_header.file_size = file_header.offset_data + capacity + bmp_info_header.height * (new_stride - _row_stride);
	}
}

void Bitmap::writeToFile(const char *file_name) 
{
	std::ofstream out_file{ file_name, std::ios_base::binary };
	if (out_file) {
		if (bmp_info_header.bit_count == 32) {
			writeHeadersAndData(out_file);
		}
		else if (bmp_info_header.bit_count == 24) {
			if (bmp_info_header.width % 4 == 0) {
				writeHeadersAndData(out_file);
			}
			else {
				uint32_t new_stride = alignStride(4);
				std::vector<uint8_t> padding_row(new_stride - _row_stride);

				writeHeaders(out_file);

				for (int y = 0; y < bmp_info_header.height; ++y) {
					out_file.write((const char*)(image_data.data() + _row_stride * y), _row_stride);
					out_file.write((const char*)padding_row.data(), padding_row.size());
				}
			}
		}
	}
	else {
		throw std::runtime_error("Error::Bitmap::writeToFile : Failed to open file!");
	}
}

void Bitmap::setPixelColor(uint32_t x, uint32_t y, uint8_t R, uint8_t G, uint8_t B, uint8_t A) 
{
	if (x > (uint32_t)bmp_info_header.width || y > (uint32_t)bmp_info_header.height) {
		throw std::runtime_error("Error::Bitmap::setPixelColor : Coordinate point exceeds image dimensions!");
	}

	uint32_t channels = bmp_info_header.bit_count / 8;
	int32_t w = bmp_info_header.width;
	int32_t h = bmp_info_header.height;

	auto index = channels * (w * (h - y - 1) + x);

	image_data[index + 0] = B;
	image_data[index + 1] = G;
	image_data[index + 2] = R;
	if (channels == 4) {
		image_data[index + 3] = A;
	}
}

//***************** PRIVATE ****************//

void Bitmap::writeHeaders(std::ofstream &out_file)
{
	out_file.write((const char*)&file_header, sizeof(file_header));
	out_file.write((const char*)&bmp_info_header, sizeof(bmp_info_header));
	if (bmp_info_header.bit_count == 32) {
		out_file.write((const char*)&bmp_color_header, sizeof(bmp_color_header));
	}
}

void Bitmap::writeHeadersAndData(std::ofstream &out_file)
{
	writeHeaders(out_file);
	out_file.write((const char*)image_data.data(), image_data.size());
}

uint32_t Bitmap::alignStride(uint32_t align_stride) 
{
	uint32_t new_stride = _row_stride;

	// Add 1 to stride until it is evenly divisible with given alignment
	while (new_stride % align_stride != 0)
		new_stride++;

	return new_stride;
}