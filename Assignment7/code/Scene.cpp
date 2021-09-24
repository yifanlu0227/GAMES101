//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){ // 在这个光源里均匀采样
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}


// Implementation of Path Tracing
// shade in fact
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    // the ray is definite, include origin, direction, direction_inv

    Intersection inter_p = intersect(ray);
    if(inter_p.happened == false)
        return Vector3f{0,0,0};
    // 打到光源
    if(inter_p.m->hasEmission())
        return inter_p.m->getEmission();

    // 直接光照
    // Get x, ws, NN, emit from inter. 
    // Shoot a ray from p to x
    // If the ray is not blocked in the middle
    // L_dir = emit * eval(wo, ws, N) * dot(ws, N) * dot(ws, NN) / |x-p|^2 / pdf_light
    // 光线和场景求交
    Vector3f p = inter_p.coords;
    Vector3f N = inter_p.normal; // point normal
    Material* m = inter_p.m;

    // 光源采样
    Intersection inter;
    float pdf_light;
    sampleLight(inter, pdf_light);
    Vector3f x = inter.coords;
    Vector3f NN = inter.normal; // light normal
    Vector3f emit = inter.emit; // light intensity
    Vector3f ws = normalize(x - p); // p to x , point to light
    Vector3f wo = - normalize(ray.direction); // point to eye

    Ray ray_p2x = Ray(p,ws); 
    Intersection inter_direct = intersect(ray_p2x);
    float epsilon = 0.00005;
    Vector3f L_dir = {0.f,0.f,0.f};
    float light_distance2 = dotProduct(x - p, x - p);
    if(inter_direct.distance - std::sqrt(light_distance2) > - epsilon)
        L_dir = emit * m->eval(wo, ws, N) * dotProduct(ws, N) * dotProduct(-ws, NN) / light_distance2 / pdf_light;

    // 间接光照
    // Test Russian Roulette with probability RussianRoulette wi = sample(wo, N)
    // Trace a ray r(p, wi)
    // If ray r hit a non-emitting object at q
    // L_indir = shade(q, wi) * eval(wo, wi, N) * dot(wi, N) / pdf(wo, wi, N) / RussianRoulette
    // Return L_dir + L_indir
    Vector3f L_indir = {0.0f,0.0f,0.0f};
    float Prob = ((float) rand() / (RAND_MAX));
    if(Prob > RussianRoulette){
        L_indir = {0.0f,0.0f,0.0f};
    } 
    else{
        Vector3f wi = normalize(m->sample(wo,N));
        Ray ray_q = Ray(p,wi);
        Intersection inter_q = intersect(ray_q);
        if(inter_q.happened && !inter_q.m->hasEmission()) // 注意，不可以是射向光源的光线
            L_indir = castRay(ray_q, depth+1) * m->eval(wo, wi, N) * dotProduct(wi, N) / m->pdf(wo, wi, N) / RussianRoulette;
    }
    return L_dir + L_indir;
}