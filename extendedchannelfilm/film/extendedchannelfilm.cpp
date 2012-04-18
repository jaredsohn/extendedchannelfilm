/*
   pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
  
   All Rights Reserved.
   For educational use only; commercial use expressly forbidden.
   NO WARRANTY, express or implied, for this software.
   (See file License.txt for complete license)
  
    ExtendedChannelFilm v0.23
    A PBRT film plug-in
    Created by Jared Sohn (sohn@alumni.cs.wisc.edu)
    Based on PBRT's ImageFilm plug-in
    July 3, 2005
 
   See http://www.cs.wisc.edu/~sohn/extendedchannelfilm/index.html for instructions and version history.
 */

// Extendedchannelfilm.cpp
#include <sstream>

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include "pbrt.h"
#include "color.h"
using namespace Imath;

#include "film.h"
#include "paramset.h"
#include "tonemap.h"
#include "sampling.h"
#include "primitive.h"
// ExtendedChannelFilm Declarations
class ExtendedChannelFilm : public Film {
public:
	// ExtendedChannelFilm Public Methods
	ExtendedChannelFilm::ExtendedChannelFilm(int xres, int yres,
	                     Filter *filt, const float crop[4],
			             const string &filename, bool premult,
			             int wf, const string &channelnames, const string &renameas);
	~ExtendedChannelFilm() {
		delete pixels;
		delete filter;
		delete channelNames;
		delete[] channelIndex;
		delete renameAs;
		delete nonRgbaChannelData;		
		delete[] nonRgbaChannelNames;
		delete[] nonRgbaShouldAve;
		delete[] nonRgbaChannelIndex;
		delete[] filterTable;
	}

	void AddSample(const Sample &sample, const Ray &ray,
	               const Spectrum &L, float alpha);
	void GetSampleExtent(int *xstart, int *xend,
	                     int *ystart, int *yend) const;
	void WriteImage();
private:
	// ExtendedChannelFilm Private Data
	Filter *filter;
	int writeFrequency, sampleCount;
	string filename;
	bool premultiplyAlpha;
	float cropWindow[4];
	int xPixelStart, yPixelStart, xPixelCount, yPixelCount;
	struct Pixel {
		Pixel() : L(0.f) {
			alpha = 0.f;
			weightSum = 0.f;
		}
		Spectrum L;
		float alpha, weightSum;
	};
	BlockedArray<Pixel> *pixels;
	BlockedArray<float> *nonRgbaChannelData;

	float *filterTable;
	
	std::vector<string> *channelNames;
	int *channelIndex;
	unsigned int numChannels;
	std::vector<Imf::PixelType> types;
	std::vector<float> minBounds;
	std::vector<float> maxBounds;

	std::string *nonRgbaChannelNames;	
	int *nonRgbaChannelIndex;
	int numNonRgbaChannels;
	bool *nonRgbaShouldAve;

	std::vector<string> *renameAs;

	void UpdateChannel(const string &channelname, int index, bool shouldave, const Ray &ray, int x, int y, float weight);
	float ComputeVal(const string &channelname, const Ray &ray) const;
	static std::vector<string> *Tokenize(const string &str, const string &delims);

	static void WriteExtendedChannelImage(const string &name, const std::vector<string> &channelnames,
		const std::vector<Imf::PixelType> &types, float *pixels,
		int xRes, int yRes,
		int totalXRes, int totalYRes,
		int xOffset, int yOffset);

	inline static float GetMinBound(const std::string &channelname);
	inline static float GetMaxBound(const std::string &channelname);
	inline static Imf::PixelType GetVarType(const std::string &channelname);
	inline static bool ShouldAverage(const std::string &channelname);
};

// ExtendedChannelFilm Method Definitions
ExtendedChannelFilm::ExtendedChannelFilm(int xres, int yres,
                     Filter *filt, const float crop[4],
		             const string &fn, bool premult, int wf, const string &channelnames, const string &renameas)
	: Film(xres, yres) {

	filter = filt;
	memcpy(cropWindow, crop, 4 * sizeof(float));
	filename = fn;
	premultiplyAlpha = premult;
	writeFrequency = sampleCount = wf;
	// Compute film image extent
	xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
	xPixelCount =
		max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
	yPixelStart =
		Ceil2Int(yResolution * cropWindow[2]);
	yPixelCount =
		max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);
	// Allocate film image storage
	pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
	// Precompute filter weight table
	#define FILTER_TABLE_SIZE 16
	filterTable =
		new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
	float *ftp = filterTable;
	for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
		float fy = ((float)y + .5f) * filter->yWidth /
			FILTER_TABLE_SIZE;
		for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
			float fx = ((float)x + .5f) * filter->xWidth /
				FILTER_TABLE_SIZE;
			*ftp++ = filter->Evaluate(fx, fy);
		} 
	}

	// Tokenize channelnames and place unique non-RGBA names into a map.
	channelNames = Tokenize(channelnames, ", ");
	int max = channelNames->size();
	int mapindex = 0;
	numNonRgbaChannels = max;
	numChannels = max;
	channelIndex = new int[numChannels];	
	map<string, int> channelmap;

	for (int i = 0; i < max; i++)
	{
		if ((*channelNames)[i] == "R") { numNonRgbaChannels--; continue; }
		if ((*channelNames)[i] == "G") { numNonRgbaChannels--; continue; }
		if ((*channelNames)[i] == "B") { numNonRgbaChannels--; continue; }
		if ((*channelNames)[i] == "A") { numNonRgbaChannels--; continue; }

		if (channelmap[(*channelNames)[i]] == 0)
		{
			channelmap[(*channelNames)[i]] = mapindex + 1; // We add one since we assume zero = no entry in map 
			channelIndex[i] = mapindex;
			mapindex++;
		} else
		{
			// We allow duplicate channels but we only have to store one set
			numNonRgbaChannels--;
			channelIndex[i] = -1;
		}
	}

	// Allocate and clear a place to store data found in addSample().
	if (numNonRgbaChannels != 0)
	{
		nonRgbaChannelData = new BlockedArray<float>(xPixelCount*numNonRgbaChannels, yPixelCount);
		for (int i = 0; i < xPixelCount * numNonRgbaChannels; i++)
		{
			for (int j = 0; j < yPixelCount; j++)
			{
				(*nonRgbaChannelData)(i,j) = 0.0f;
			}
		}

		// Iterate through channel names and place them into nonrgbachannels.
		// (We do this because it is much faster to iterate through an array than a map.)
		// 
		// Also, determine and cache indices and which should be averaged.
		nonRgbaChannelNames = new std::string[numNonRgbaChannels];
		nonRgbaShouldAve = new bool[numNonRgbaChannels];
		nonRgbaChannelIndex = new int[numNonRgbaChannels];

		mapindex = 0;
		for (map<string, int>::iterator it = channelmap.begin(); it != channelmap.end(); it++)
		{
			nonRgbaChannelNames[mapindex] = (*it).first;
			nonRgbaChannelIndex[mapindex] = (*it).second - 1; 
			nonRgbaShouldAve[mapindex] = ShouldAverage((*it).first);
			mapindex++;
		}
	} else
	{
		nonRgbaChannelData = NULL;
		nonRgbaChannelNames = NULL;
		nonRgbaChannelIndex = NULL;
		nonRgbaShouldAve = NULL;
	}

	// We must ensure that all channel names as written to the output EXR file are unique.
	//
	// To do this, we append [pbrt<indexno>] to the beginning of non-unique channel names.
	//
	// This strategy fails, however, when we have the following situation:
	//
	// renameas channels = {R, pbrt3-R, R}.  Here, the third channel will initially be called 
	// pbrt3-R, but since that name is not unique, it will then be called pbrt3-pbrt3-R. 
	// (If that name also exists, then we continue appending pbrt3- to the name until it is unique.
	//
	// Of course, having a variable with such a long nonuseful name is not ideal, but
	// one shouldn't name their variables with pbrt[x] anyway.
	renameAs = Tokenize(renameas, ", ");
	unsigned int numwords = renameAs->size();
	map<string, int> renameashash;

	std::ostringstream stream;
	std::string temp2;
	std::string temp;

	for (unsigned int i = 0; i < numwords; i++)
	{
		temp = (*renameAs)[i];
		
		// If name begins with 'pbrt_' then show a warning and quit
		if (temp.find("pbrt",0) != string::npos)
		{
			Warning("Channel names should not begin with the string 'pbrt'");
			exit(1);
		}

		if (renameashash[temp] != 0)	
		{
			stream << i;
			temp2 = stream.str();
			stream.str("");	

			temp = "pbrt" + temp2 + "_" + temp;
			Warning("Channel '%s' was renamed to '%s' since output channel names must be unique.", (*renameAs)[i].c_str(), temp.c_str());
			(*renameAs)[i] = temp;
		}

		renameashash[temp]++;
	}
	
	// We set the number of channels to min(# channels, # renameas)
	if (numChannels < (*renameAs).size())
	{
		Warning("Not going to rename all channels since only %i of them were initially specified.", numChannels);		
		(*renameAs).erase(renameAs->begin() + numChannels - 1, renameAs->begin()+ (*renameAs).size() - 1);
	} else if (numChannels > (*renameAs).size())
	{		
		Warning("Only saving the first %i channels because others were not renamed.", numChannels);
		numChannels = (*renameAs).size();
	}

	// Here we precompute the minimum, maximum, and type of each channel.
	// We also determine if samples should be averaged.
	for (unsigned int i = 0; i < numChannels; i++)
	{
		std::string c = (*channelNames)[i];

		types.push_back(GetVarType(c));
		minBounds.push_back(GetMinBound(c));
		maxBounds.push_back(GetMaxBound(c));
	}
}

// This method is a threadsafe version of strtok which returns all of the results
// immediately within a vector.
//
// str = The string to tokenize.
// delims = Each character in this string represents a delimiter that will be used for tokenizing.
// Courtesy of John Danks
std::vector<string> *ExtendedChannelFilm::Tokenize(const string &str, const string &delims)
{
        std::vector<string> *vec = new std::vector<string>;
        if (str.size() == 0 || delims.size() == 0) {
                Warning("Channel string or delimiter string was zero-length!");
                exit(1);
        }
        string::size_type pos = 0;
        while (pos != string::npos) {
                string::size_type pos2 = str.find_first_of(delims, pos);
                if (pos2 == string::npos) {
                        // no more delimiters
                        if (pos != str.size()) {
                                // no trailing delimiters, pick up the last word
                                vec->push_back(str.substr(pos));
                        }
                        break;
                }
                if (pos2 > pos) {
                        // not multiple delimiters, pick up current word
                        vec->push_back(str.substr(pos, pos2-pos));
                }
                // advance pos to next delimiter
                pos = pos2+1;
        }
        return vec;
}

void ExtendedChannelFilm::AddSample(const Sample &sample,
		const Ray &ray, const Spectrum &L, float alpha) {
						
	// Compute sample's raster extent
	float dImageX = sample.imageX - 0.5f;
	float dImageY = sample.imageY - 0.5f;
	int x0 = Ceil2Int (dImageX - filter->xWidth);
	int x1 = Floor2Int(dImageX + filter->xWidth);
	int y0 = Ceil2Int (dImageY - filter->yWidth);
	int y1 = Floor2Int(dImageY + filter->yWidth);
	x0 = max(x0, xPixelStart);
	x1 = min(x1, xPixelStart + xPixelCount - 1);
	y0 = max(y0, yPixelStart);
	y1 = min(y1, yPixelStart + yPixelCount - 1);
	// Loop over filter support and add sample to pixel arrays
	// Precompute $x$ and $y$ filter table offsets
	int *ifx = (int *)alloca((x1-x0+1) * sizeof(int));
	for (int x = x0; x <= x1; ++x) {
		float fx = fabsf((x - dImageX) *
		                 filter->invXWidth * FILTER_TABLE_SIZE);
		ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
	}
	int *ify = (int *)alloca((y1-y0+1) * sizeof(int));
	for (int y = y0; y <= y1; ++y) {
		float fy = fabsf((y - dImageY) *
		                 filter->invYWidth * FILTER_TABLE_SIZE);
		ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
	}
	for (int y = y0; y <= y1; ++y)
		for (int x = x0; x <= x1; ++x) {
			// Evaluate filter value at $(x,y)$ pixel
			int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
			float filterWt = filterTable[offset];
			// Update pixel values with filtered sample contribution
			Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
			pixel.L.AddWeighted(filterWt, L);
			pixel.alpha += alpha * filterWt;
			pixel.weightSum += filterWt;

			// Also write non-RGBA channels
			for (int i = 0; i < numNonRgbaChannels; i++)
			{
				UpdateChannel(nonRgbaChannelNames[i], nonRgbaChannelIndex[i], nonRgbaShouldAve[i], ray, x, y, filterWt);
			}
		}

	// Possibly write out in-progress image
	if (--sampleCount == 0) {
		WriteImage();
		sampleCount = writeFrequency;
	}
}

void ExtendedChannelFilm::UpdateChannel(const string &channelname, int index, bool shouldave, const Ray &ray, int x, int y, float weight)
{
	float val = ComputeVal(channelname, ray);

	if (shouldave) 
	{ 		
		(*nonRgbaChannelData)((x * numNonRgbaChannels) + (index), y) += val * weight;
	} else
	{
		// We will arbitrarily choose the last sample that we looked at.  (Ideally we would
		// keep a histogram of all values and choose the mode...)
		(*nonRgbaChannelData)((x * numNonRgbaChannels) + (index), y) = val;
	}
}

// Look up the channel name and return the appropriate value.  If the channel name 
// is not found, we show a warning and quit.
float ExtendedChannelFilm::ComputeVal(const string &channelname, const Ray &ray) const
{
	std::string s = channelname;

	transform (s.begin(), s.end(), s.begin(), tolower);
	if (s == "z") return (float) ray.maxt;
	if (s == "wx") return ray.p.x;
	if (s == "wy") return ray.p.y;
	if (s == "wz") return ray.p.z;
	if (s == "nx") return ray.nn.x;
	if (s == "ny") return ray.nn.y;
	if (s == "nz") return ray.nn.z;
	if (s == "ou") return ray.u;
	if (s == "ov") return ray.v;
	if (s == "id") return ray.id;

	// The specified channel does not exist.  Exit and show a warning.
	std::cout << "\n";
	Warning("Channel '%s' does not exist.  Exiting...", channelname.c_str());

	exit(1);
}
		
void ExtendedChannelFilm::GetSampleExtent(int *xstart,
		int *xend, int *ystart, int *yend) const {
	*xstart = Floor2Int(xPixelStart + .5f - filter->xWidth);
	*xend   = Floor2Int(xPixelStart + .5f + xPixelCount  +
		filter->xWidth);
	*ystart = Floor2Int(yPixelStart + .5f - filter->yWidth);
	*yend   = Floor2Int(yPixelStart + .5f + yPixelCount +
		filter->yWidth);
}

float ExtendedChannelFilm::GetMinBound(const std::string &channelname)
{
	if (channelname == "R") return 0.0f;
	if (channelname == "G") return 0.0f;
	if (channelname == "B") return 0.0f;
	if (channelname == "A") return 0.0f;

	return -(FLT_MAX + 1);
}

float ExtendedChannelFilm::GetMaxBound(const std::string &channelname)
{
	if (channelname == "A") return 1.0f;

	return INFINITY;
}

Imf::PixelType ExtendedChannelFilm::GetVarType(const std::string &channelname)
{
	if (channelname == "Z") return Imf::FLOAT;

	return Imf::HALF;
}

bool ExtendedChannelFilm::ShouldAverage(const std::string &channelname)
{
	if (channelname == "id") 
	{
		return false;
	}
	return true;
}

// Compute final pixel values and write to an EXR file
void ExtendedChannelFilm::WriteImage()
{
	int nPix = xPixelCount * yPixelCount;

	float *channeldata = new float[numChannels*nPix];
	float *alpha = new float[1*nPix];

	int offset = 0;
	for (int y = 0; y < yPixelCount; ++y) {
		for (int x = 0; x < xPixelCount; ++x) {
			// Convert pixel spectral radiance to RGB
			float xyz[3];
			(*pixels)(x, y).L.XYZ(xyz);
			const float
				rWeight[3] = { 3.240479f, -1.537150f, -0.498535f };
			const float
				gWeight[3] = {-0.969256f,  1.875991f,  0.041556f };
			const float
				bWeight[3] = { 0.055648f, -0.204043f,  1.057311f };					

			for (unsigned int i = 0; i < numChannels; i++)
			{
				if ((*channelNames)[i] == "R")
				{
					channeldata[numChannels * offset + i] = rWeight[0]*xyz[0] +
						rWeight[1]*xyz[1] +
						rWeight[2]*xyz[2];
				} else if ((*channelNames)[i] == "G")
				{
					channeldata[numChannels * offset + i] = gWeight[0]*xyz[0] +
									gWeight[1]*xyz[1] +
									gWeight[2]*xyz[2];
				} else if ((*channelNames)[i] == "B")
				{
					channeldata[numChannels * offset + i] = bWeight[0]*xyz[0] +
							bWeight[1]*xyz[1] +
							bWeight[2]*xyz[2];
				} else if ((*channelNames)[i] == "A")
				{
					channeldata[numChannels * offset + i] = (*pixels)(x, y).alpha;
				} else
				{
					int index = channelIndex[i];
					assert(index != -1);					
					channeldata[numChannels * offset + i] = (*nonRgbaChannelData)(x * numNonRgbaChannels + index, y);
				}
			}			
		
			// Now we have copied all of our data into channels

			// Normalize pixel with weight sum
			float weightSum = (*pixels)(x, y).weightSum;
			float invWt;
			if (weightSum != 0)
			{
				invWt = 1 / weightSum;
			} else
			{
				invWt = 1.0f;
			}

			// Clamp the data within our bounds and divide by the sum of the weights (if any...)
			for (unsigned int i = 0; i < numChannels; i++)
			{										
				std::string c = (*channelNames)[i];
				
				float temp;
				if (ShouldAverage(c))
				{
					temp = invWt;
				} else
				{
					temp = 1.0f;
				}

				channeldata[numChannels * offset + i] = Clamp(channeldata[numChannels*offset + i] * temp, minBounds[i], maxBounds[i]);

				if ((premultiplyAlpha) && (c == "A"))
				{
					// We make a potentially unnecessary copy of alpha (since we must ensure that we have it...)
					alpha[offset] = Clamp((*pixels)(x, y).alpha, 0.f, 1.f);
				}
			}

			// Compute premultiplied alpha color
			if (premultiplyAlpha) 
			{				
				for (unsigned int i = 0; i < numChannels; i++)
				{
					if ((*channelNames)[i] == "R") channeldata[numChannels * offset + i] *= alpha[offset];
					if ((*channelNames)[i] == "G") channeldata[numChannels * offset + i] *= alpha[offset];
					if ((*channelNames)[i] == "B") channeldata[numChannels * offset + i] *= alpha[offset];
				}			
			}
			++offset;
		}
	}

	// Write image
	WriteExtendedChannelImage(filename, *renameAs, types, channeldata, 
		xPixelCount, yPixelCount,
		xResolution, yResolution,
		xPixelStart, yPixelStart);
	
	// Release temporary image memory
	delete[] channeldata;
	delete[] alpha;
}

// Note that this method would fit better in exrio.cpp.  However, we place it here so that we can reduce by
// two the number of main pbrt source files that must be altered for this plug-in (exrio.cpp and pbrt.h).
//
// If this function does get added to exrio.cpp, then one must be careful to either remove Imf::PixelType from
// the method signature or add the EXR directory to the include path list for every plug-in that #includes pbrt.h.
void ExtendedChannelFilm::WriteExtendedChannelImage(const string &name, const std::vector<string> &channelnames,
		const std::vector<Imf::PixelType> &types, float *pixels,
		int xRes, int yRes,
		int totalXRes, int totalYRes,
		int xOffset, int yOffset)
{
	Imf::Header header(totalXRes, totalYRes);
	Box2i dataWindow(V2i(xOffset, yOffset), V2i(xOffset + xRes - 1, yOffset + yRes - 1));
	header.dataWindow() = dataWindow;
	
	int numchannels = channelnames.size();
	for (int i = 0; i < numchannels; i++)
	{
		header.channels().insert(channelnames[i].c_str(), Imf::Channel (types[i]));
	}

	// Copy all data for each channel into a single array and then
	// add that into the framebuffer
	Imf::FrameBuffer fb;

	half **hptrs = new half*[numchannels];
	unsigned int **iptrs = new unsigned int *[numchannels];
	float **fptrs = new float*[numchannels];	

	for (int i = 0; i < numchannels; i++)
	{
		hptrs[i] = NULL;
		fptrs[i] = NULL;
		iptrs[i] = NULL;
	}

	for (int i = 0; i < numchannels; i++)
	{
		if (types[i] == Imf::HALF)
		{
			half *h = new half[xRes * yRes];
			hptrs[i] = h;
			
			for (int j = 0; j < xRes * yRes; j++)
			{
				h[j] = (half) pixels[j * numchannels + i];
			}

			h -= (xOffset + yOffset * xRes);

			fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)h, sizeof(half),
				xRes*sizeof(half)));

		} else if ((types[i]) == Imf::FLOAT)
		{
			float *f = new float[xRes * yRes];
			fptrs[i] = f;

			for (int j = 0; j < xRes * yRes; j++)
			{
				f[j] = pixels[j * numchannels + i];
			}

			f -= (xOffset + yOffset * xRes);

			fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)f, sizeof(float),
				xRes*sizeof(float)));

		} else
		{
			unsigned int *ui = new unsigned int[xRes * yRes];
			iptrs[i] = ui;

			for (int j = 0; j < xRes * yRes; j++)
			{
				ui[j] = (unsigned int) pixels[j * numchannels + i];
			}

			ui -= (xOffset + yOffset * xRes);

			fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)ui, sizeof(unsigned int),
				xRes*sizeof(unsigned int)));

		} 
	}

	Imf::OutputFile file(name.c_str(), header);

	file.setFrameBuffer(fb);
	try {
		file.writePixels(yRes);
	}
	catch (const std::exception &e) {
		Error("Unable to write image file \"%s\": %s", name.c_str(),
			e.what());
	}
	// Free the memory that was allocated for the slices
	for (int i = 0; i < numchannels; i++)
	{
		delete[] fptrs[i];
		delete[] iptrs[i];
		delete[] hptrs[i];
	}
	delete[] fptrs;
	delete[] iptrs;
	delete[] hptrs;
}

extern "C" DLLEXPORT Film *CreateFilm(const ParamSet &params, Filter *filter)
{
	string filename = params.FindOneString("filename", "pbrt.exr");
	bool premultiplyAlpha = params.FindOneBool("premultiplyalpha", true);
	string channelnames = params.FindOneString("channels", "R,G,B,A");
	string renameas = params.FindOneString("renameas", channelnames);

	int xres = params.FindOneInt("xresolution", 640);
	int yres = params.FindOneInt("yresolution", 480);
	float crop[4] = { 0, 1, 0, 1 };
	int cwi;
	const float *cr = params.FindFloat("cropwindow", &cwi);
	if (cr && cwi == 4) {
		crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
		crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
		crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
		crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
	}
	int writeFrequency = params.FindOneInt("writefrequency", -1);

	return new ExtendedChannelFilm(xres, yres, filter, crop,
		filename, premultiplyAlpha, writeFrequency, channelnames, renameas);
}
