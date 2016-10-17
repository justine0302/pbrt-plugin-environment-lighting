
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H
#define PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H

// lights/MedianCutEnvironmentLight.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
    #include <iostream>

struct float2
{
    float x, y;
};

// Median cut environment light Declarations
class MedianCutEnvironmentLight : public Light {
public:
    // MedianCutEvironmentLight Public Methods
    MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
        const string &texmap);
    ~MedianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; } //if  false -> importance sampling
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;
private:
    // MedianCutEvironmentLight Private Data
    MIPMap<RGBSpectrum> *radianceMap;
    Distribution2D *distribution;
    std::vector<float2> lights;
    std::vector<RGBSpectrum> spectrums;
};


class summed_area_table{

public:
    int width_, height_;
    std::vector<float> sat_;

    float I(int x, int y) const
        {
                if (x < 0 || y < 0) return 0;
                int i = y*width_+ x;
                return sat_[i];
        }

        void create_sat(float *img, int width, int height){
            width_ = width; height_ = height;
            int i=0;
            for(int y=0;y<height_;y++){
                for( int x=0; x<width_; x++){
                    float ixy=img[i];
                    float value=ixy+I(x-1,y)+I(x,y-1)-I(x-1,y-1);
                    sat_.push_back(value);
                    i++;
                }
            }
        }

        float sum(int ax, int ay, int bx, int by, int cx, int cy, int dx, int dy) const
        {
                return I(cx, cy) + I(ax, ay) - I(bx, by) - I(dx, dy);
        }

};

struct sat_region
{
    int x_, y_, w_, h_;
    float sum_;
    const summed_area_table* sat_;

    void create(int x, int y, int w, int h, const summed_area_table* sat, float init_sum = -1)
    {
        x_ = x; y_ = y; w_ = w; h_ = h; sum_ = init_sum; sat_ = sat;

        if (sum_ < 0)
            sum_ = sat_->sum(x,       y, 
                    x+(w-1), y, 
                    x+(w-1), y+(h-1),
                    x,       y+(h-1));
    }


    void split_w(sat_region& A, sat_region& B) const
    {
        for (size_t w = 1; w <= w_; ++w)
        {
            A.create(x_, y_, w, h_, sat_);
            if (A.sum_*2.f >= sum_)
                break;
        }
        B.create(x_ + (A.w_-1), y_, w_ - A.w_, h_, sat_, sum_ - A.sum_);
    }


    void split_h(sat_region& A, sat_region& B) const
    {
        for (size_t h = 1; h <= h_; ++h)
        {
            A.create(x_, y_, w_, h, sat_);
            if (A.sum_*2.f >= sum_)
                break;
        }
        B.create(x_, y_ + (A.h_-1), w_, h_ - A.h_, sat_, sum_ - A.sum_);
    }

    float2 centroid() const
    {
        float2 c;

        sat_region A;
            
        for (size_t w = 1; w <= w_; ++w)
        {
            A.create(x_, y_, w, h_, sat_);
            if (A.sum_*2.f >= sum_)
                break;
        }
        c.x = A.x_ + (A.w_-1);

        for (size_t h = 1; h <= h_; ++h)
        {
            A.create(x_, y_, w_, h, sat_);
            if (A.sum_*2.f >= sum_)
                break;
        }
        c.y = A.y_ + (A.h_-1);
            
        return c;
    }

};



MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H
