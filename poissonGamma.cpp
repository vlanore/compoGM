#include "model.hpp"

// interface to get a Real
class Real_I {
  public:
    virtual double getValue() const = 0;
};

// little object to encapsulate having a constant OR a pointer to Real_I
class RealProp {
    int mode{0};  // 0:unset, 1:constant, 2:pointer
    Real_I *ptr{nullptr};
    double value{0};

  public:
    RealProp() = default;
    explicit RealProp(double value) : mode(1), value(value) {}
    explicit RealProp(Real_I *ptr) : mode(2), ptr(ptr) {}
    double getValue() const {
        if (mode == 1) {
            return value;
        } else if (mode == 2) {
            return ptr->getValue();
        } else {
            std::cerr << "Property is no set!\n";
            exit(1);
        }
    }
};

class Exponential_C : public Component, public Real_I {
    RealProp lambda{};
    double value{0.0};

  public:
    explicit Exponential_C(double value = 0.0) : value(value) {
        port("lambdaConst", &Exponential_C::setLambda<double>);
        port("lambdaPtr", &Exponential_C::setLambda<Real_I *>);
    };

    std::string _debug() const override {
        std::stringstream ss;
        ss << "Exponential(" << lambda.getValue() << "):" << value;
        return ss.str();
    }
    double getValue() const override { return value; }

    template <class... Args>
    void setLambda(Args... args) {
        lambda = RealProp(std::forward<Args>(args)...);
    }
};

class Gamma_C : public Component {
  public:
    Real_I *k, *theta;
    std::string _debug() const override { return "Gamma"; }
};

int main() {
    Assembly model;

    model.component<Exponential_C>("Exp", 2.0);
    model.component<Exponential_C>("Exp2");
    model.property("Exp", "lambdaConst", 1.0);
    model.connect<UseProvide<Real_I>>("Exp2", "lambdaPtr", "Exp");

    model.instantiate();
    model.print_all();
}
