
struct M3 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
        map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
        map<string, double>& size_factors) {
        // global variables
        m.component<OrphanNormal>("a0", 1, -2, 2);
        m.component<OrphanNormal>("a1", 1, 0, 2);
        m.component<OrphanExp>("sigma_alpha", 1, 1);

        m.component<Matrix<OrphanNormal>>("log10(q)", genes, conditions, 1, 2, 2);  // changed
        m.connect<MapPower10>("log10(q)", "q");

        m.component<Array<Mean>>("q_bar", genes).connect<ArrayToValueMatrixLines>("parent", "q");

        m.component<Array<DeterministicTernaryNode<double>>>("log10(alpha_bar)", genes,
             [](double a0, double a1, double q_bar) { return log10(a0 + a1 / q_bar); })
            .connect<ArrayToValue>("a", "a0")
            .connect<ArrayToValue>("b", "a1")
            .connect<ArrayToValueArray>("c", "q_bar");

        m.component<Array<Normal>>("log10(alpha)", genes, 1)
            .connect<ArrayToValueArray>("a", "log10(alpha_bar)")
            .connect<ArrayToValue>("b", "sigma_alpha");
        m.connect<MapInversePower10>("log10(alpha)", "1/alpha");

        m.component<Matrix<GammaSR>>("tau", genes, samples, 1)
            .connect<MatrixLinesToValueArray>("a", "1/alpha")
            .connect<MatrixLinesToValueArray>("b", "1/alpha");

        m.component<Array<Constant<double>>>("sf", samples, 0)
            .connect<SetArray<double>>("x", size_factors);

        m.component<Matrix<DeterministicTernaryNode<double>>>(
             "lambda", genes, samples, [](double a, double b, double c) { return a * b * c; })
            .connect<MatrixColumnsToValueArray>("a", "sf")
            .connect<ManyToMany<ArraysMap<UseValue>>>("b", "q", condition_mapping)
            .connect<MatrixToValueMatrix>("c", "tau");

        m.component<Matrix<Poisson>>("K", genes, samples, 0)
            .connect<SetMatrix<int>>("x", counts)
            .connect<MatrixToValueMatrix>("a", "lambda");
    }
};
