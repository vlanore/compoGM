#include <random>
#include "model.hpp"

class Go : public Component {
  public:
    Go() { port("go", &Go::go); }
    virtual void go() = 0;
};

// interface to get a Real
class Real {
  public:
    virtual double getValue() const = 0;
    virtual void setValue(double value) = 0;
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
  protected:
    RealProp param{};
    double value{0.0};

  public:
    std::string name{};

    // constructors
    UnaryReal() = delete;
    explicit UnaryReal(const std::string &name) : name(name) {
        port("paramConst", &UnaryReal::setParam<double>);
        port("paramPtr", &UnaryReal::setParam<Real *>);
        port("sample", &UnaryReal::sample);
    };

    // methods required from parent classes
    std::string _debug() const override {
        std::stringstream ss;
        ss << name << "(" << param.getValue() << "):" << value;
        return ss.str();
    }
    double getValue() const override { return value; }
    void setValue(double valuein) override { value = valuein; }

    // class-specific methods
    template <class... Args>
    void setParam(Args... args) {
        param = RealProp(std::forward<Args>(args)...);
    }

    virtual void sample() = 0;
};

class Exponential : public UnaryReal {
  public:
    explicit Exponential(double value = 0.0) : UnaryReal("Exponential") { setValue(value); }

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::exponential_distribution<> d(param.getValue());
        setValue(d(gen));
    }
};

class Gamma : public UnaryReal {
  public:
    explicit Gamma() : UnaryReal("Gamma") {}

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::gamma_distribution<> d(param.getValue(), param.getValue());
        setValue(d(gen));
    }
};

class Poisson : public UnaryReal {
  public:
    explicit Poisson() : UnaryReal("Poisson") {}

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::poisson_distribution<> d(param.getValue());
        setValue(d(gen));
    }
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
    void setValue(double) override {
        std::cerr << "-- Warning! Trying to set a deterministic node!\n";
    }
    std::string _debug() const override {
        std::stringstream ss;
        ss << "Product(" << a.getValue() << "," << b.getValue() << "):" << getValue();
        return ss.str();
    }
};

class SimpleMove : public Go {
    Real *target{nullptr};

  public:
    SimpleMove() { port("target", &SimpleMove::setTarget); }

    std::string _debug() const override { return "SimpleMove"; }
    void go() override {}

    void setTarget(Real *targetin) { target = targetin; }
};

class MCMCScheduler : public Go {
    std::vector<SimpleMove *> moves;

  public:
    MCMCScheduler() { port("register", &MCMCScheduler::registerMove); }

    void go() override { std::cout << "MCMCScheduler started!\n"; }
    std::string _debug() const override { return "MCMCScheduler"; }

    void registerMove(SimpleMove *ptr) { moves.push_back(ptr); }
};

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
