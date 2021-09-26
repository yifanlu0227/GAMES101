#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        
        // create particle
        for(int i=0;i<num_nodes;i++){
            Vector2D node = start * i / num_nodes + end * (num_nodes-i) / num_nodes;
            Mass * tmp = new Mass(node, node_mass, false);
            masses.push_back(tmp);
        }
        // create spring
        for(int i=0;i<num_nodes-1;i++){
            Spring *string = new Spring(masses[i], masses[i+1], k);
            springs.push_back(string);
        }
        // Comment-in this part when you implement the constructor
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D m1_m2 = (s->m1->position - s->m2->position);
            s->m1->forces += - s->k *  m1_m2 / m1_m2.norm() * (m1_m2.norm() - s->rest_length);
            Vector2D m2_m1 = (s->m2->position - s->m1->position); 
            s->m2->forces += - s->k *  m2_m1 / m2_m1.norm() * (m2_m1.norm() - s->rest_length);
        }

        for (auto &m : masses)
        {   
            float k_d = 0.03;
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += m->mass * gravity;
                // TODO (Part 2): Add global damping
                m->forces += - k_d * m->velocity;

                m->velocity += delta_t * m->forces / m->mass;
                m->position += delta_t * m->velocity;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D m1_m2 = (s->m1->position - s->m2->position);
            s->m1->forces += - s->k *  m1_m2 / m1_m2.norm() * (m1_m2.norm() - s->rest_length);
            Vector2D m2_m1 = (s->m2->position - s->m1->position); 
            s->m2->forces += - s->k *  m2_m1 / m2_m1.norm() * (m2_m1.norm() - s->rest_length);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                Vector2D last_position = m->last_position;
                // TODO (Part 3.1): Set the new position of the rope mass
                m->forces += m->mass * gravity;

                // TODO (Part 4): Add global Verlet damping
                float damping_factor = 0.00005;
                m->position = temp_position + (1-damping_factor) * (temp_position - last_position) + 
                                        m->forces / m->mass * delta_t * delta_t;
                m->last_position = temp_position;

            }
            m->forces =  Vector2D(0,0);
        }
    }
}
