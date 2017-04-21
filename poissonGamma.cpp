#include "model.hpp"

class Real_I {
  public:
    virtual double getValue() = 0;
};

class Constant_C : public Instance, public Real_I {
  public:
    double value;
    Constant_C(double value) : value(value) {}
    double getValue() override { return value; }
    std::string _sdebug() override { return "Constant"; }
};

class Exponential_C : public Instance, public Real_I {
  public:
    Real_I *lambda;
    double value;
    std::string _sdebug() override { return "Exponential"; }
    double getValue() override { return value; }
};

class Gamma_C : public Instance {
public:
    Real_I *k, *theta;
    std::string _sdebug() override { return "Gamma"; }
};

class PoissonGamma_A : public Assembly {
  public:
    PoissonGamma_A() {
        node<Constant_C>("One", 1.0);
        node<Exponential_C>("Theta");
        node<Exponential_C>("Sigma");
        connection<UseProvide<Exponential_C,Real_I>>("Theta", "One", &Exponential_C::lambda);
        connection<UseProvide<Exponential_C,Real_I>>("Sigma", "One", &Exponential_C::lambda);

        node<Array<Gamma_C>>("Omega", 5, [](int){ return Gamma_C();});
        connection<MultiProvideArray<Gamma_C,Real_I>>("Omega", "Theta", &Gamma_C::theta);
    }

    void go() {
        instantiate();
        print_all();
        to_dot();
    }
};

int main() {
    PoissonGamma_A model{};
    model.go();
}
