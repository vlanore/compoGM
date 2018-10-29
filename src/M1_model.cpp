
struct M1 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
        map<string, map<string, int>>& counts, map<string, string>& condition_mapping) {
        m.component<Matrix<OrphanNormal>>("log10(lambda)", genes, conditions, 1, 3, pow(1.5, 2));
        m.connect<MapPower10>("log10(lambda)", "lambda");

        m.component<Matrix<Poisson>>("K", genes, samples, 0)
            .connect<SetMatrix<int>>("x", counts)
            .connect<ManyToMany<ArraysMap<UseValue>>>("a", "lambda", condition_mapping);
    }
};