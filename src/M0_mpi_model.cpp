
struct M0 : public Composite {
    static void contents(
        Model& m, IndexSet& experiments, IndexSet& samples, map<string, map<string, int>>& data) {
        m.component<OrphanExp>("alpha", 1, 10);
        m.component<OrphanExp>("mu", 1, 1);

        m.component<Array<Gamma>>("lambda", experiments, 10)
            .connect<ArrayToValue>("a", "alpha")
            .connect<ArrayToValue>("b", "mu");

        if (p.rank != 0) {  // slave only
            m.component<Matrix<Poisson>>("K", experiments, samples, 0)
                .connect<MatrixLinesToValueArray>("a", "lambda")
                .connect<SetMatrix<int>>("x", data);
        }
    }
};
