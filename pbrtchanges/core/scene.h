
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_SCENE_H
#define PBRT_SCENE_H
// scene.h*
#include "pbrt.h"
#include "primitive.h"
#include "transport.h"
// Scene Declarations
class COREDLL Scene {
public:
	// Scene Public Methods
	void Render();
	Scene(Camera *c, SurfaceIntegrator *in,
		VolumeIntegrator *vi, Sampler *s,
		Primitive *accel, const vector<Light *> &lts,
		VolumeRegion *vr);
	~Scene();
	bool Intersect(const Ray &ray, Intersection *isect) const {
		bool results = aggregate->Intersect(ray, isect);

		ray.nn = isect->dg.nn;
		ray.p = isect->dg.p;
		ray.u = isect->dg.u;
		ray.v = isect->dg.v;
		if (isect->primitive != NULL)
		{
			ray.id = isect->primitive->id;
		} else
		{
			ray.id = 0;
		}

		return results;
	}
	bool IntersectP(const Ray &ray) const {
		return aggregate->IntersectP(ray);
	}
	const BBox &WorldBound() const;
	Spectrum Li(const RayDifferential &ray, const Sample *sample,
		float *alpha = NULL) const;
	Spectrum Transmittance(const Ray &ray) const;
	// Scene Data
	Primitive *aggregate;
	vector<Light *> lights;
	Camera *camera;
	VolumeRegion *volumeRegion;
	SurfaceIntegrator *surfaceIntegrator;
	VolumeIntegrator *volumeIntegrator;
	Sampler *sampler;
	BBox bound;
};
#endif // PBRT_SCENE_H
