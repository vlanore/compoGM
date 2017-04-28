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
    model.component<SimpleMove>("Move1");
    model.connect<UseProvide<Real>>("Move1", "target", "Theta");

    model.component<SimpleMove>("Move2");
    model.connect<UseProvide<Real>>("Move2", "target", "Sigma");

    model.component<Array<SimpleMove, 5>>("MoveArray");
    model.connect<ArrayOneToOne<Real>>("MoveArray", "target", "Omega");

    model.component<MCMCScheduler>("Scheduler");
    model.connect<UseProvide<SimpleMove>>("Scheduler", "register", "Move1");
    model.connect<UseProvide<SimpleMove>>("Scheduler", "register", "Move2");
    model.connect<MultiUse<SimpleMove>>("Scheduler", "register", "MoveArray");

    // instantiate everything!
    model.instantiate();

    // do some things
    model.call("Theta", "sample");
    model.call("Scheduler", "go");


    model.print_all();
}
