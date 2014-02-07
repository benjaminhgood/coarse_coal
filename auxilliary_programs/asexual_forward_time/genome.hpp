#ifndef GENOME_HPP
#define GENOME_HPP

typedef unsigned long Label;

class LabelGenerator{
    private:
        Label current_label;
    public:
        LabelGenerator(): current_label(0) {};
        Label get_next_label() {return current_label++;};
};

class Mutation{
    public:
        Label label;
        double s;
        double W_background;
};

//typedef Label Mutation;

#include<vector>
typedef std::vector<Mutation> MutationList;
#endif
