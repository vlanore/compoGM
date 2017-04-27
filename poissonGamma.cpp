#include "model.hpp"

// interface to get a Real
class Real {
  public:
    virtual double getValue() const = 0;
};

// little object to encapsulate having a constant OR a pointer to Real
class RealProp {
    int mode{0};  // 0:unset, 1:constant, 2:pointer
    Real *ptr{nullptr};
    double value{0};

  public:
    RealProp() = default;
    explicit RealProp(double value) : mode(1), value(value) {}
    explicit RealProp(Real *ptr) : mode(2), ptr(ptr) {}
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

class UnaryReal : public Component, public Real {
    RealProp param{};
    double value{0.0};

  public:
    std::string name{};

    UnaryReal() = delete;
    explicit UnaryReal(const std::string &name) : name(name) {
        port("paramConst", &UnaryReal::setLambda<double>);
        port("paramPtr", &UnaryReal::setLambda<Real *>);
    };

    std::string _debug() const override {
        std::stringstream ss;
        ss << name << "(" << param.getValue() << "):" << value;
        return ss.str();
    }
    double getValue() const override { return value; }
    void setValue(double valuein) { value = valuein; }

    template <class... Args>
    void setLambda(Args... args) {
        param = RealProp(std::forward<Args>(args)...);
    }
};

class Exponential : public UnaryReal {
  public:
    explicit Exponential(double value = 0.0) : UnaryReal("Exponential") { setValue(value); }
};

class Gamma : public UnaryReal {
  public:
    explicit Gamma() : UnaryReal("Gamma") {}
};

class Poisson : public UnaryReal {
  public:
    explicit Poisson() : UnaryReal("Poisson") {}
};

class Product : public Component, public Real {
    RealProp a{};
    RealProp b{};

  public:
    Product() {
        port("aPtr", &Product::setA<Real *>);
        port("bPtr", &Product::setB<Real *>);
        port("bConst", &Product::setB<double>);
    }

    template <class... Args>
    void setA(Args... args) {
        a = RealProp(std::forward<Args>(args)...);
    }

    template <class... Args>
    void setB(Args... args) {
        b = RealProp(std::forward<Args>(args)...);
    }

    void setA(Real *ptr) { a = RealProp(ptr); }
    double getValue() const override { return a.getValue() * b.getValue(); }
    std::string _debug() const override {
        std::stringstream ss;
        ss << "Product(" << a.getValue() << "," << b.getValue() << "):" << getValue();
        return ss.str();
    }
};

int main() {
    Assembly model;

    model.component<Exponential>("Sigma");
    model.property("Sigma", "paramConst", 1.0);

    model.component<Exponential>("Theta");
    model.property("Theta", "paramConst", 1.0);

    model.component<Array<Gamma, 5>>("Omega");
    model.connect<MultiUse<Real>>("Omega", "paramPtr", "Theta");

    model.component<Array<Product, 5>>("rate");
    model.connect<ArrayOneToOne<Real>>("rate", "aPtr", "Omega");
    model.connect<MultiUse<Real>>("rate", "bPtr", "Sigma");

    model.component<Array<Poisson, 5>>("X");
    model.connect<ArrayOneToOne<Real>>("X", "paramPtr", "rate");

    model.instantiate();
    model.print_all();
}
