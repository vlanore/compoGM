#include "graphicalModel.hpp"

int main() {
    Assembly model;

    // graphical model part
    model.component<Exponential>("Sigma");
    model.property("Sigma", "paramConst", 1.0);

    model.component<Exponential>("Theta");
    model.property("Theta", "paramConst", 1.0);

    model.component<Array<Gamma, 5>>("Omega");
    model.connect<MultiProvide<Real>>("Omega", "paramPtr", "Theta");

    model.component<Array<Product, 5>>("rate");
    model.connect<ArrayOneToOne<Real>>("rate", "aPtr", "Omega");
    model.connect<MultiProvide<Real>>("rate", "bPtr", "Sigma");

    model.component<Array<Poisson, 5>>("X");
    model.connect<ArrayOneToOne<Real>>("X", "paramPtr", "rate");

    // moves part
    model.component<MultiSample>("Sampler");
    model.connect<UseProvide<RandomNode>>("Sampler", "register", "Sigma");
    model.connect<UseProvide<RandomNode>>("Sampler", "register", "Theta");
    model.connect<MultiUse<RandomNode>>("Sampler", "register", "Omega");
    model.connect<MultiUse<RandomNode>>("Sampler", "register", "X");

    model.component<RejectionSampling>("RS", 10000);
    model.connect<UseProvide<Go>>("RS", "sampler", "Sampler");
    model.connect<MultiUse<RandomNode>>("RS", "data", "X");

    // model.component<SimpleMove>("Move1");
    // model.connect<UseProvide<RandomNode>>("Move1", "target", "Theta");

    // model.component<SimpleMove>("Move2");
    // model.connect<UseProvide<RandomNode>>("Move2", "target", "Sigma");

    // model.component<Array<SimpleMove, 5>>("MoveArray");
    // model.connect<ArrayOneToOne<RandomNode>>("MoveArray", "target", "Omega");

    // model.component<Scheduler>("Scheduler");
    // model.connect<UseProvide<SimpleMove>>("Scheduler", "register", "Move1");
    // model.connect<UseProvide<SimpleMove>>("Scheduler", "register", "Move2");
    // model.connect<MultiUse<SimpleMove>>("Scheduler", "register", "MoveArray");

    // instantiate everything!
    model.instantiate();

    // hacking the model to clamp observed data
    auto &Xref = model.get_ref<ComponentArray>("X");
    Xref.get_ref_at<RandomNode>(0).clamp(1);
    Xref.get_ref_at<RandomNode>(1).clamp(0);
    Xref.get_ref_at<RandomNode>(2).clamp(1);
    Xref.get_ref_at<RandomNode>(3).clamp(0);
    Xref.get_ref_at<RandomNode>(4).clamp(0);

    // do some things
    model.call("RS", "go");

    model.print_all();
}
